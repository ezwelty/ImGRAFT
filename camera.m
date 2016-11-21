% TODO: Use degrees for yaw, pitch, and roll angles.
% TODO: Split out DEM intersection function into DEM class method.
% TODO: Limit voxelviewshed to camera view (not just position)
% TODO: Remove imgsz from fullmodel (since not needed in optimization)
% TODO: Clean up optimizecam. Add support for named arguments.
% TODO: Extend optimizecam to work with camera bundle.

classdef camera
  % camera Distorted camera model
  %
  % This class is an implementation of the distorted camera model used by OpenCV:
  % http://docs.opencv.org/modules/calib3d/doc/camera_calibration_and_3d_reconstruction.html
  %
  % Note: Pixel coordinates are defined such that [1, 1] is the center of the
  % upper left pixel of the image and [nx, ny] is the center of the lower right
  % pixel, where nx and ny are the width and height of the image. As a result,
  % the image extends from [0.5, 0.5] to [nx + 0.5, ny + 0.5] and the center is
  % at [(nx + 1) / 2, (ny + 1) / 2].
  %
  % camera Properties:
  % imgsz    - Size of image in pixels [nx|ncols|width, ny|nrows|height]
  % f        - Focal length in pixels [fx, fy]
  % c        - Camera center coordinates in pixels [cx, cy]
  % k        - Radial distortion coefficients [k1, ..., k6]
  % p        - Tangential distortion coefficients [p1, p2]
  % xyz      - Camera position in world coordinates [x, y, z]
  % viewdir  - Camera view direction in degrees [yaw, pitch, roll]
  %            yaw: clockwise rotation about z-axis (0 = look north)
  %            pitch: rotation from horizon (+ look up, - look down)
  %            roll: rotation about optical axis (+ down right, - down left, from behind)
  % sensorsz - Size of camera sensor in mm [width, height] (optional)
  %
  % camera Properties (dependent):
  % fullmodel - Vector containing all 20 camera parameters
  %             [xyz(1:3), imgsz(1:2), viewdir(1:3), f(1:2), c(1:2), k(1:6), p(1:2)]
  % fmm       - Focal length in mm [fx, fy] (not set unless sensorsz is defined)
  % R         - Rotation matrix corresponding to camera view direction (read-only)
  % K         - Camera matrix [fx 0 cx; 0 fy cy; 0 0 1] (read-only)
  %
  % camera Methods:
  % camera      - Construct a new camera object
  % project     - Project world coordinates to image coordinates (3D -> 2D)
  % invproject  - Project image coordinates to world coordinates (2D -> 3D)
  % optimizecam - Optimize camera parameters to mimimize the distance between
  %               projected world coordinates and their expected image coordinates
  %
  % ImGRAFT - An image georectification and feature tracking toolbox for MATLAB
  % Copyright (C) 2014 Aslak Grinsted (<www.glaciology.net glaciology.net>)

  % Permission is hereby granted, free of charge, to any person obtaining a copy
  % of this software and associated documentation files (the "Software"), to deal
  % in the Software without restriction, including without limitation the rights
  % to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  % copies of the Software, and to permit persons to whom the Software is
  % furnished to do so, subject to the following conditions:
  %
  % The above copyright notice and this permission notice shall be included in
  % all copies or substantial portions of the Software.
  %
  % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  % IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  % FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  % AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  % LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  % OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  % THE SOFTWARE.

  properties
    xyz
    viewdir
    f
    c
    k
    p
    imgsz
    sensorsz
  end

  properties (Dependent)
    R
    K
    fmm
    fullmodel
    framebox
    framepoly
  end

  methods

    function cam = camera(varargin)
      % CAMERA  Construct a new camera object.
      %
      % There are several ways to call this method -
      %
      % 1. Initialize with default parameters, then edit individual parameters:
      %
      %   cam = camera()
      %   cam.viewdir = [pi 0 0] % look west
      %
      % 2. Specify camera parameters as a list:
      %
      %   cam = camera(xyz, imgsz, viewdir, f, c, k, p, sensorsz)
      %
      % 3. Specify camera parameters as named arguments:
      %
      %   cam = camera('xyz', [0 0 10], 'viewdir', [pi 0 0])
      %
      % 4. Specify camera parameters in a structure with matching field names:
      %
      %   S.xyz = [0 0 10]; S.viewdir = [pi 0 0];
      %   cam = camera(S)
      %
      % 5. Specify camera parameters as a 20-element or shorter (fullmodel) vector:
      %
      %   cam = camera([1 1 0 1024 768 pi 0 0 1000 1000 512 384 0 0 0 0 0 0 0 0])

      % Single vector argument: build camera from fullmodel vector
      if nargin == 1 && isnumeric(varargin{1})
        cam.fullmodel = varargin{1};
        return
      end

      % All other cases: set and validate camera parameters individually
      p = inputParser;
      p.CaseSensitive = false;
      p.StructExpand = true;
      p.addOptional('xyz', [0 0 0], @(x) isnumeric(x) && length(x) <= 3);
      p.addOptional('imgsz', [], @(x) isnumeric(x) && length(x) <= 2);
      p.addOptional('viewdir', [0 0 0], @(x) isnumeric(x) && length(x) <= 3);
      p.addOptional('f', [], @(x) isnumeric(x) && length(x) <= 2);
      p.addOptional('c', [], @(x) isnumeric(x) && length(x) == 2);
      p.addOptional('k', [0 0 0 0 0 0], @(x) isnumeric(x) && length(x) <= 6);
      p.addOptional('p', [0 0], @(x) isnumeric(x) && length(x) <= 2);
      p.addOptional('sensorsz', [], @(x) isnumeric(x) && length(x) <= 2);
      % HACK: Removes non-matching field names from structure
      if nargin == 1 && isstruct(varargin{1})
        fields = fieldnames(varargin{1});
        values = struct2cell(varargin{1});
        include = ismember(fields, p.Parameters);
        varargin{1} = cell2struct(values(include), fields(include));
      end
      p.parse(varargin{:});

      % Set parameters
      for field = fieldnames(p.Results)'
        cam.(field{1}) = p.Results.(field{1});
      end

      % Validate parameters
      % TODO: Move to camera.set functions?
      % xyz: 3-element vector, default = 0
      %cam.xyz(end + 1:3) = 0;
      % imgsz: 2-element vector, no default, expand [#] to [#, #]
      if length(cam.imgsz) == 1, cam.imgsz(end + 1) = cam.imgsz(end); end
      % viewdir: 3-element vector, default = 0
      cam.viewdir(end + 1:3) = 0;
      % f: 2-element vector, no default, expand [f] to [f, f]
      if length(cam.f) == 1, cam.f(end + 1) = cam.f(end); end
      % c: 2-element vector, default = (imgsz[2 1] + 1) / 2
      if isempty(cam.c) && ~isempty(cam.imgsz)
        cam.c = (cam.imgsz + 1) / 2;
        cam.c(cam.imgsz == 0) = 0;
      end
      % k: 6-element vector, default = 0
      cam.k(end + 1:6) = 0;
      % p: 2-element vector, default = 0
      cam.p(end + 1:2) = 0;
      % sensorsz: 2-element vector, no default, expand [#] to [#, #]
      if length(cam.sensorsz) == 1, cam.sensorsz(end + 1) = cam.sensorsz(end); end
    end

    function value = get.R(cam)
      % Initial rotations of camera reference frame
      % (camera +z pointing up, with +x east and +y north)
      % Point camera north: -90 deg counterclockwise rotation about x-axis
      %   ri = [1 0 0; 0 cosd(-90) sind(-90); 0 -sind(-90) cosd(-90)];
      % (camera +z now pointing north, with +x east and +y down)

      % View direction rotations
      C = cosd(cam.viewdir); S = sind(cam.viewdir);
      % syms C1 C2 C3 S1 S2 S3
      % yaw: counterclockwise rotation about y-axis (relative to north, from above: +cw, - ccw)
      %   ry = [C1 0 -S1; 0 1 0; S1 0 C1];
      % pitch: counterclockwise rotation about x-axis (relative to horizon: + up, - down)
      %   rp = [1 0 0; 0 C2 S2; 0 -S2 C2];
      % roll: counterclockwise rotation about z-axis (from behind camera: + ccw, - cw)
      %   rr = [C3 S3 0; -S3 C3 0; 0 0 1];

      % Apply all rotations in order
      %   R = rr * rp * ry * ri;
      value = [ C(1) * C(3) + S(1) * S(2) * S(3),  C(1) * S(2) * S(3) - C(3) * S(1), -C(2) * S(3); ...
                C(3) * S(1) * S(2) - C(1) * S(3),  S(1) * S(3) + C(1) * C(3) * S(2), -C(2) * C(3); ...
                C(2) * S(1)                     ,  C(1) * C(2)                     ,  S(2)       ];
    end

    function cam = set.R(cam, value)
      cam.viewdir = rot2oblang(value);
    end

    function value = get.fullmodel(cam)
      fullmodel = [cam.xyz, cam.imgsz, cam.viewdir, cam.f, cam.c, cam.k, cam.p];
      if length(fullmodel) == 20
        value = [cam.xyz, cam.imgsz, cam.viewdir, cam.f, cam.c, cam.k, cam.p];
      else
        error('Camera parameters missing or not of the correct length.')
      end
    end

    function cam = set.fullmodel(cam, value)
      if length(value) < 20, error('Camera.fullmodel must have 20 elements.'), end
      cam.xyz = value(1:3);
      cam.imgsz = value(4:5); cam.viewdir = value(6:8); cam.f = value(9:10);
      cam.c = value(11:12); cam.k = value(13:18); cam.p = value(19:20);
    end

    function value = get.K(cam)
      value = [cam.f(1) 0 cam.c(1); 0 cam.f(2) cam.c(2); 0 0 1];
    end

    function cam = set.fmm(cam, value)
      if (isempty(cam.sensorsz))
        error('Camera sensor size not set.');
      else
        cam.f = value .* cam.imgsz ./ cam.sensorsz;
      end
    end

    function value = get.fmm(cam)
      if (isempty(cam.sensorsz))
        error('Camera sensor size not set.');
      else
        value = cam.f .* cam.sensorsz ./ cam.imgsz;
      end
    end

    function value = get.framebox(cam)
      value = [0 cam.imgsz(1) 0 cam.imgsz(2)] + 0.5;
    end

    function value = get.framepoly(cam)
      value = flipud([0 0 ; cam.imgsz(1) 0 ; cam.imgsz(1) cam.imgsz(2) ; 0 cam.imgsz(2) ; 0 0] + 0.5); % flip so cw in typical xy-space
    end

    function plot(cam, radius, color)
      plot3(cam.xyz(1), cam.xyz(2), cam.xyz(3), 'r.'); hold on;
      corner_uv = [0 0; cam.imgsz(1) 0; cam.imgsz; 0 cam.imgsz(2)] + 0.5;
      corner_xyz = cam.invproject(corner_uv);
      directions = bsxfun(@minus, corner_xyz, cam.xyz);
      if nargin < 2
        radius = 1;
      end
      scaled_directions = bsxfun(@times, directions, radius ./ sum(directions.^2, 2));
      corner_xyz = bsxfun(@plus, scaled_directions, cam.xyz);
      plot3(corner_xyz([1:4 1], 1), corner_xyz([1:4 1], 2), corner_xyz([1:4 1], 3), 'k-');
      plot3([repmat(cam.xyz(1), 4, 1) corner_xyz(:, 1)]', [repmat(cam.xyz(2), 4, 1) corner_xyz(:, 2)]', [repmat(cam.xyz(3), 4, 1) corner_xyz(:, 3)]', 'k-');
      if nargin < 3
        color = rand(1, 3);
      end
      patch(corner_xyz(:, 1), corner_xyz(:, 2), corner_xyz(:, 3), color, 'FaceAlpha', 0.5);
    end

    function value = viewpyramid(cam, radius)
      u = [0:cam.imgsz(1)]'; v = [0:cam.imgsz(2)]';
      edge_uv = [
        [u repmat(0, length(u), 1)];
        [repmat(cam.imgsz(1), length(v), 1) v];
        [flip(u) repmat(cam.imgsz(2), length(u), 1)];
        [repmat(0, length(v), 1) flip(v)]
      ];
      edge_xyz = cam.invproject(edge_uv);
      directions = bsxfun(@minus, edge_xyz, cam.xyz);
      if nargin < 2
        radius = 1;
      end
      scaled_directions = bsxfun(@times, directions, radius ./ sum(directions.^2, 2));
      edge_xyz = bsxfun(@plus, scaled_directions, cam.xyz);
      value = [edge_xyz; cam.xyz];
    end

    function value = viewbox(cam, radius)
      if nargin < 2
        pyramid = cam.viewpyramid();
      else
        pyramid = cam.viewpyramid(radius);
      end
      value = [min(pyramid); max(pyramid)];
    end

    function cam = resize(cam, scale)
      % RESIZE  Resize a camera image.
      %
      %   cam = cam.resize(scale)
      %
      % The focal length (f) and principal point (c) are adjusted as needed to
      % account for the change in image size (imgsz).
      %
      % NOTE: Since an image can only contain whole pixels, the resulting image size
      % may not equal the target image size. Because of this rounding, calls to
      % this function may be nonreversible.
      %
      % Inputs:
      %   scale - Scale factor [s] or desired image size [ny|nrows, nx|ncols]

      if length(scale) > 1
        imgsz1 = round(scale(1:2));
        scale = imgsz1 ./ cam.imgsz;
        if any(abs(diff(scale)) * cam.imgsz > 1)
          error(['Aspect ratio of target image size (' num2str(imgsz1(1) / imgsz1(2), 2) ') too different from original (' num2str(cam.imgsz(1) / cam.imgsz(2), 2) ') to be accounted for by rounding.']);
        end
      else
        imgsz1 = round(scale * cam.imgsz);
        scale = imgsz1 ./ cam.imgsz;
      end
      f = cam.f .* flip(scale);
      c = ((imgsz1 + 1) / 2) + (cam.c - ((cam.imgsz + 1) / 2)) .* flip(scale);
      % Save to new object
      cam = cam;
      cam.f = f;
      cam.c = c;
      cam.imgsz = imgsz1;
    end

    function [uv, depth, inframe] = project(cam, xyz)
      % PROJECT  Project 3D world coordinates to 2D image coordinates.
      %
      %   [uv, depth, inframe] = cam.project(xyz)
      %
      % Inputs:
      %   xyz - World coordinates [x1 y1 z1; x2 y2 z2; ...]
      %
      % Outputs:
      %   uv      - Image coordinates [u1 v1; u2 v2; ...]
      %   depth   - Distance of each point from the camera
      %   inframe - Boolean whether each point is in the image
      %
      % See also: camera.invproject

      % Convert to camera coordinates
      xyz = bsxfun(@minus, xyz, cam.xyz);
      xyz = xyz * cam.R';

      % Normalize by perspective division
      xy = bsxfun(@rdivide, xyz(:, 1:2), xyz(:, 3));

      % Apply distortion
      infront = xyz(:, 3) > 0;
      xy(infront, :) = cam.distort(xy(infront, :));

      % Convert to image coordinates
      % [x', y'] = [fx * x + cx, fy * y + cy]
      % (set points behind camera to NaN)
      uv = nan(size(xy));
      uv(infront, :) = [cam.f(1) * xy(infront, 1) + cam.c(1), cam.f(2) * xy(infront, 2) + cam.c(2)];

      if nargout > 1
        depth = xyz(:, 3);
      end
      if nargout > 2
        % TODO: additional constraint for negative k1 and r2 > 1. (See orthorectification example) NOTE: Why? A point is either in or out of the image frame.
        inframe = (depth > 0) & (uv(:, 1) >= 0.5) & (uv(:, 1) <= cam.imgsz(1) + 0.5) & (uv(:, 2) >= 0.5) & (uv(:, 2) <= cam.imgsz(2) + 0.5);
      end
    end

    function xyz = invproject(cam, uv, S, xy0)
        % INVPROJECT  Project 2D image coordinates into 3D world coordinates.
        %
        %   xyz = cam.invproject(uv)
        %   xyz = cam.invproject(uv, X, Y, Z)
        %   xyz = cam.invproject(uv, X, Y, Z, xy0)
        %
        % Image coordinates are projected out the camera as rays. If a surface is specified, the intersections of the rays with the surface are returned (and NaN if intersection failed). If a surface is not specified, points of unit distance along the rays are returned.
        % There are three ways to call this method:
        %
        % 1. Return a set of 3 coordinates which when projected
        %    results in the pixel coordinates (uv).
        %
        %   xyz = cam.invproject(uv)
        %
        % 2. Inverse projection constrained to DEM surface.
        %    (Project DEM to camera-view and use interpolation to find xyz.)
        %
        %   xyz=cam.invproject(uv,X,Y,Z)
        %
        %
        % 3. Minimize misfit between projected DEM point and UV. (least squares)
        %
        %   xyz=cam.invproject(uv,X,Y,Z,xy0)
        %
        %
        %
        % Inputs:
        %   uv: 2 column matrix with pixel coordinates.
        %   [X,Y,Z]: DEM. (X,Y expected to be consistent with output from meshgrid.)
        %   [xy0]: initial guess of some xy points on the DEM which are
        %          consistent with pixel coordinates in uv.
        %
        % Outputs:
        %   xyz: 3-column matrix with world coordinates consistent with
        %        pixel coordinates in uv.
        %
        %
        nanix = any(isnan(uv), 2);
        anynans = any(nanix);
        if anynans
          uv(nanix, :) = [];
        end

        npts = size(uv, 1);

        % Convert to normalized camera coordinates
        % [x, y] = [(x' - cx) / fx, (y' - cy) / fy]
        xy = [(uv(:, 1) - cam.c(1)) / cam.f(1), (uv(:, 2) - cam.c(2)) / cam.f(2)];
        xy = cam.undistort(xy);
        % Convert to ray direction vectors [dx dy dz]
        v = [xy ones(npts, 1)] * cam.R;

        if nargin == 2
          % No surface: Return points along rays ~ unit distance from camera
          % v = bsxfun(@rdivide, v, sqrt(sum(v.^2, 2))); % force unit distance (slower)
          xyz = bsxfun(@plus, v, cam.xyz);

        elseif nargin > 2 && isnumeric(S) && length(S) == 4
          % Plane: Return intersection of rays with plane
          xyz = intersectRayPlane(cam.xyz, v, S);

        elseif nargin > 2 && strcmp(class(S), 'DEM')
          % DEM: Return intersection of rays with DEM
          rays = [repmat(cam.xyz, npts, 1) v];
          xyz = nan(size(rays, 1), 3);
          for i = 1:size(rays, 1)
            xyz(i, :) = intersectRayDEM(rays(i, :), S);
          end

        %   % DEM: Return intersection of rays with DEM
        %   [X, Y, Z] = deal(S.X, S.Y, S.Z);
        %   visible = voxelviewshed(X, Y, Z, cam.xyz);
        %   Z = Z ./ visible;
        %   xyz = nan(npts, 3);
        %   if nargin < 4
        %       [uv0, ~ , inframe] = cam.project([X(visible(:)), Y(visible(:)), Z(visible(:))]);
        %       uv0(:,3) = X(visible(:));
        %       uv0(:,4) = Y(visible(:));
        %       uv0(:,5) = Z(visible(:));
        %       uv0 = uv0(inframe, :);
        %       if exist('scatteredInterpolant','file') > 1
        %           Xscat = scatteredInterpolant(uv0(:,3), uv0(:,4), uv0(:,3));
        %           Xscat.Points = uv0(:,1:2);
        %           Yscat = Xscat; Yscat.Values = uv0(:,4);
        %           Zscat = Xscat; Zscat.Values = uv0(:,5);
        %       else
        %           %fallback for older versions of matlab.
        %           Xscat = TriScatteredInterp(uv0(:,3), uv0(:,4), uv0(:,3));  %#ok<REMFF1>
        %           Xscat.X = uv0(:,1:2);
        %           Yscat = Xscat; Yscat.V = uv0(:,4);
        %           Zscat = Xscat; Zscat.V = uv0(:,5);
        %       end
        %       xy0 = [Xscat(uv(:,1), uv(:,2)) Yscat(uv(:,1), uv(:,2)) Zscat(uv(:,1), uv(:,2))];
        %       xyz = xy0;
        %
        %       if anynans
        %           xyz(find(~nanix), :) = xyz; %find necessary because it ensures that xyz can grow.
        %           xyz(find(nanix), :) = nan;
        %       end
        %       return
        %   end
        %
        %   if Y(2, 2) < Y(1, 1)
        %       X = flipud(X); Y = flipud(Y); Z = flipud(Z);
        %   end
        %   if X(2, 2) < X(1, 1)
        %       X = fliplr(X); Y = fliplr(Y); Z = fliplr(Z);
        %   end
        %
        %   if exist('griddedInterpolant', 'file') > 1
        %       zfun = griddedInterpolant(X', Y', Z'); %TODO: improve robustness.
        %   else
        %       %fallback for older versions of matlab. slower
        %       zfun = @(x, y) interp2(X, Y, Z, x, y);
        %   end
        %   for ii = 1:length(uv)
        %       %misfit=@(xy)sum((cam.project([xy zfun(xy(1),xy(2))])-uv(ii,1:2)).^2);
        %       misfitlm = @(xy) (cam.project([xy(:)' zfun(xy(1), xy(2))]) - uv(ii, 1:2))'.^2;
        %       try
        %           %[xyz(ii,1:2),err]=fminunc(misfit,xy0(ii,1:2),optimset('LargeScale','off','Display','off','TolFun',0.001)); %TODO: remove dependency. can i use LMFnlsq?
        %           xyz(ii, 1:2) = LMFnlsq(misfitlm, xy0(ii, 1:2));
        %           xyz(ii, 3) = zfun(xyz(ii, 1), xyz(ii, 2));
        %           if sum(misfitlm(xyz(ii, 1:2))) > 2^2
        %               xyz(ii, :) = nan; %do not accept greater than 2 pixel error.
        %           end
        %       catch
        %       end
        %   end
        end

        if anynans
            xyz(find(~nanix), :) = xyz; % find() necessary because it ensures that xyz can grow.
            xyz(find(nanix), :) = NaN;
        end
    end

    function [newcam, rmse, aic] = optimizecam(cam, xyz, uv, freeparams)
        % OPTIMIZECAM  Calibrate a camera from paired image-world coordinates.
        %
        %   [newcam, rmse, aic] = cam.optimizecam(xyz, uv, freeparams)
        %
        % Uses an optimization routine to minize the root-mean-square reprojection
        % error of image-world point correspondences (xyz, uv) by adjusting the
        % specified camera parameters.
        %
        % If uv has three columns, the third column is interpreted as a weight
        % in the misfit function.
        %
        % Inputs:
        %   xyz        - World coordinates [x1 y1 z1; x2 y2 z2; ...]
        %   uv         - Image coordinates [u1 v1; u2 v2; ...]
        %                (optional 3rd column may specify weights)
        %   freeparams - 20-element vector describing which camera parameters
        %                should be optimized. Follows same order as cam.fullmodel.
        %
        % Outputs:
        %   newcam - The optimized camera
        %   rmse   - Root-mean-square reprojection error
        %   aic    - Akaike information criterion for reprojection errors, which
        %            can help determine an appropriate degree of complexity for
        %            the camera model (i.e. avoid overfitting).
        %            NOTE: Only strictly applicable for unweighted fitting.
        %
        % Example:
        %   % Optimize the 3 viewdir parameters (elements 6-8 in fullmodel).
        %   [newcam, rmse, aic] = cam.optimizecam(xyz, uv, '00000111000000000000')

        % Drop any points with NaN coordinates
        nanrows = any(isnan(xyz), 2) | any(isnan(uv), 2);
        xyz(nanrows, :) = [];
        uv(nanrows, :) = [];

        fullmodel0 = cam.fullmodel; % Describes the initial camera being perturbed.

        % Convert fullmodel to boolean
        freeparams = ~(freeparams(:) == 0 | freeparams(:) == '0')';
        paramix = find(freeparams);
        Nfree = length(paramix);
        mbest = zeros(1, Nfree);

        newcam = @(m) camera(fullmodel0 + sparse(ones(1, Nfree), paramix, m, 1, length(fullmodel0)));

        if size(uv, 2) == 3
            % Weighted least squares
            misfit = @(m) reshape((project(newcam(m), xyz) - uv(:, 1:2)) .* uv(:, [3 3]), [], 1);
        else
            % Unweighted least squares
            misfit = @(m) reshape(project(newcam(m), xyz) - uv, [], 1);
        end
        if isnan(misfit(mbest))
            error('All GCPs must be in front of the initial camera location for optimizecam to work.'); %TODO: write better explanation. and remove requirement.
        end
        [mbest, RSS] = LMFnlsq(misfit, mbest); %WORKS SUPER FAST

        Nuv = size(uv, 1);
        newcam = newcam(mbest);
        rmse = sqrt(RSS / Nuv);
        aic = numel(uv) * log(RSS / numel(uv)) + 2 * Nfree;
    end

    % Sky: Horizon detection (works great!)
    % function X = horizon(cam, dem, ddeg)
    %   if nargin < 3
    %     ddeg = 1;
    %   end
    %   viewedges = cam.viewpyramid();
    %   dxy = bsxfun(@minus, viewedges(1:(end - 1), 1:2), viewedges(end, 1:2));
    %   angles = atan2d(dxy(:, 2), dxy(:, 1));
    %   ray_angles = [min(angles):ddeg:max(angles)]';
    %   X = dem.horizon(cam.xyz, ray_angles);
    % end

    function in = inframe(cam, uv)
      box = cam.framebox;
      in = uv(:, 1) >= box(1) & uv(:, 1) <= box(2) & uv(:, 2) >= box(3) & uv(:, 2) <= box(4);
    end

    % function uv = idealize(cam, uv)
    %   % Convert to normalized camera coordinates
    %   % [x, y] = [(x' - cx) / fx, (y' - cy) / fy]
    %   xy = cam.normalize(uv);
    %   % Project back undistorted
    %   uv = [cam.f(1) * xy(:, 1) + cam.c(1), cam.f(2) * xy(:, 2) + cam.c(2)];
    % end

    function xy = normalize(cam, pts)
      if size(pts, 2) == 2
        % Convert image to camera coordinates
        % [x, y] = [(x' - cx) / fx, (y' - cy) / fy]
        xy = [(pts(:, 1) - cam.c(1)) / cam.f(1), (pts(:, 2) - cam.c(2)) / cam.f(2)];
        % Remove distortion
        xy = cam.undistort(xy);
      elseif size(pts, 2) == 3
        % Convert world to camera coordinates
        xyz = bsxfun(@minus, pts, cam.xyz);
        xyz = xyz * cam.R';
        % Normalize by perspective division
        xy = bsxfun(@rdivide, xyz(:, 1:2), xyz(:, 3));
        % Convert points behind camera to NaN
        infront = xyz(:, 3) > 0;
        xy(infront, :) = NaN;
      end
    end

    function uv = normalized2image(cam, xy, distort)
      if nargin < 3
        distort = true;
      end
      if distort
        xy = cam.distort(xy);
      end
      uv = [cam.f(1) * xy(:, 1) + cam.c(1), cam.f(2) * xy(:, 2) + cam.c(2)];
    end

    function dxyz = normalized2world(cam, xy)
      dxyz = [xy ones(size(xy, 1), 1)] * cam.R;
    end

    function h = plotDistortion(cam, scale)
      if nargin < 2
        scale = 1;
      end
      du = 100; nu = round(cam.imgsz(1) / du); du = cam.imgsz(2) / nu;
      u = (0:du:cam.imgsz(1)) + 0.5;
      v = (0:du:cam.imgsz(2)) + 0.5;
      [pu pv] = meshgrid(u, v);
      P0 = [pu(:) pv(:)]; P1 = cam.idealize(P0);
      h = plot(cam.framepoly(:, 1), cam.framepoly(:, 2), 'k:'); hold on;
      quiver(P0(:, 1), P0(:, 2), scale * (P1(:, 1) - P0(:, 1)), scale * (P1(:, 2) - P0(:, 2)), 0, 'r');
      set(gca, 'xlim', cam.imgsz(1) * [-0.1 1.1], 'ylim', cam.imgsz(2) * [-0.1 1.1], 'ydir', 'reverse');
      axis equal; hold off;
    end

  end % methods

  methods (Static)

    function R = optimizeR(xyA, xyB)
      ang = [0 0 0];
      misfit = @(ang) camera.Re(ang, xyA, xyB);
      [ang, RSS] = LMFnlsq(misfit, ang); %WORKS SUPER FAST
      R = ang2rot(ang);
    end

    function e = Re(ang, xyA, xyB)
      R = ang2rot(ang);
      pxyB = hnormalise(R * homog(xyA'))';
      e = sum((pxyB(:, 1:2) - xyB).^2, 2);
    end

  end

  methods (Access = private)

    function xy = distort(cam, xy)
      % DISTORT  Apply radial and tangential lens distortion to normalized camera coordinates.
      %
      %   xy = distort(xy)
      %
      % Inputs:
      %   xy - normalized camera coordinates [x1 y1; x2 y2; ...]
      %
      % Outputs:
      %   xy - distorted normalized camera coordinates [x1 y1; x2 y2; ...]

      if any([cam.k, cam.p] ~= 0)
        % r = sqrt(x^2 + y^2)
        r2 = sum(xy.^2, 2);
        if any(cam.k ~= 0)
          % Radial lens distortion
          % dr = (1 + k1 * r^2 + k2 * r^4 + k3 * r^6) / (1 + k4 * r^2 + k5 * r^4 + k6 * r^6)
          dr = 1 + cam.k(1) * r2 + cam.k(2) * r2.^2 + cam.k(3) * r2.^3;
          if any(cam.k(4:6) ~= 0)
            dr = dr ./ (1 + cam.k(4) * r2 + cam.k(5) * r2.^2 + cam.k(6) * r2.^3);
          end
        end
        if any(cam.p ~= 0)
          % Tangential lens distortion
          % dtx = 2xy * p1 + p2 * (r^2 + 2x^2)
          % dty = p1 * (r^2 + 2y^2) + 2xy * p2
          xty = xy(:, 1) .* xy(:, 2);
          dtx = 2 * xty * cam.p(1) + cam.p(2) * (r2 + 2 * xy(:, 1).^2);
          dty = cam.p(1) * (r2 + 2 * xy(:, 2).^2) + 2 * xty * cam.p(2);
        end
        % Compute distorted camera coordinates
        % x' = dr * x + dtx
        % y' = dr * y + dty
        if any(cam.k ~= 0)
          xy = bsxfun(@times, xy, dr);
        end
        if any(cam.p ~= 0)
          xy = xy + [dtx dty];
        end
      end
    end

    function xy = undistort(cam, xy)
      % UNDISTORT  Undo radial and tangential lens distortion on normalized camera coordinates.
      %
      %   xy = undistort(xy)
      %
      % Inputs:
      %   xy - distorted normalized camera coordinates [x1 y1; x2 y2; ...]
      %
      % Outputs:
      %   xy - normalized camera coordinates [x1 y1; x2 y2; ...]

      if any([cam.k, cam.p] ~= 0)

        % May fail for large negative k1.
        if cam.k(1) < -0.5
          warning(['Large, negative k1 (', num2str(cam.k(1), 3), '). Undistort may fail.'])
        end

        % If only k1 is nonzero, use closed form solution.
        % Cubic roots solution from Numerical Recipes in C 2nd Edition:
        % http://apps.nrbook.com/c/index.html (pages 183-185)
        if sum([cam.k cam.p] ~= 0) == 1 && cam.k(1) ~= 0
          phi = atan2(xy(:, 2), xy(:, 1));
          Q = -1 / (3 * cam.k(1));
          R = -xy(:, 1) ./ (2 * cam.k(1) * cos(phi));
          % For negative k1
          if cam.k(1) < 0
            th = acos(R ./ sqrt(Q^3));
            r = -2 * sqrt(Q) * cos((th - 2 * pi) / 3);
          % For positive k1
          else
            A = (sqrt(R.^2 - Q^3) - R).^(1/3);
            B = Q * (1 ./ A);
            r = (A + B);
          end
          xy = [r .* cos(phi), r .* sin(phi)];
          xy = real(xy);

        % Otherwise, use iterative solution.
        else
          xyi = xy; % initial guess
          for n = 1:20
            % r = sqrt(x^2 + y^2)
            r2 = sum(xy.^2, 2);
            if any(cam.k ~= 0)
              % Radial lens distortion
              % dr = (1 + k1 * r^2 + k2 * r^4 + k3 * r^6) / (1 + k4 * r^2 + k5 * r^4 + k6 * r^6)
              dr = 1 + cam.k(1) * r2 + cam.k(2) * r2.^2 + cam.k(3) * r2.^3;
              if any(cam.k(4:6) ~= 0)
                dr = dr ./ (1 + cam.k(4) * r2 + cam.k(5) * r2.^2 + cam.k(6) * r2.^3);
              end
            end
            if any(cam.p ~= 0)
              % Tangential lens distortion
              % dtx = 2xy * p1 + p2 * (r^2 + 2x^2)
              % dty = p1 * (r^2 + 2y^2) + 2xy * p2
              xty = xy(:, 1) .* xy(:, 2);
              dtx = 2 * xty * cam.p(1) + cam.p(2) * (r2 + 2 * xy(:, 1).^2);
              dty = cam.p(1) * (r2 + 2 * xy(:, 2).^2) + 2 * xty * cam.p(2);
            end
            % Remove distortion
            % x = (x' - dtx) / dr
            % y = (y' - dty) / dr
            if any(cam.p ~= 0)
              xy = xyi - [dtx dty];
              if any(cam.k ~= 0)
                xy = bsxfun(@rdivide, xy, dr);
              end
            else
              xy = bsxfun(@rdivide, xyi, dr);
            end
          end
        end
      end
    end

  end % methods (private)

end % classdef
