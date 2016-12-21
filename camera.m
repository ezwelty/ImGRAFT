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
    fmm
    fullmodel
    framebox
    framepoly
  end

  methods

    % Camera creation

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
    end

    function cam = set.xyz(cam, value)
      % xyz: 3-element vector, default = 0
      value(end + 1:3) = 0;
      cam.xyz = value;
    end

    function cam = set.imgsz(cam, value)
      % imgsz: 2-element vector, no default, expand [#] to [#, #]
      if length(value) == 1
        value(end + 1) = value(end);
      end
      cam.imgsz = value;
      % c: 2-element vector, default = (imgsz[2 1] + 1) / 2
      if isempty(cam.c) && ~isempty(cam.imgsz)
        cam.c = (cam.imgsz + 1) / 2;
        cam.c(cam.imgsz == 0) = 0;
      end
    end

    function cam = set.viewdir(cam, value)
      % viewdir: 3-element vector, default = 0
      value(end + 1:3) = 0;
      cam.viewdir = value;
    end

    function cam = set.f(cam, value)
      % f: 2-element vector, no default, expand [f] to [f, f]
      if length(value) == 1
        value(end + 1) = value(end);
      end
      cam.f = value;
    end

    function cam = set.k(cam, value)
      % k: 6-element vector, default = 0
      value(end + 1:6) = 0;
      cam.k = value;
    end

    function cam = set.p(cam, value)
      % p: 2-element vector, default = 0
      value(end + 1:2) = 0;
      cam.p = value;
    end

    function cam = set.sensorsz(cam, value)
      % sensorsz: 2-element vector, no default, expand [#] to [#, #]
      if length(value) == 1
        value(end + 1) = value(end);
      end
      cam.sensorsz = value;
    end

    function cam = set.fmm(cam, value)
      if isempty(cam.sensorsz)
        error('Camera sensor size not set.');
      else
        cam.f = value .* cam.imgsz ./ cam.sensorsz;
      end
    end

    function cam = set.R(cam, value)
      % cos(elevation) != 0
      if abs(value(3, 3)) ~= 1
        w = asin(value(3, 3));
        p = atan2(value(3, 1) / cos(w), value(3, 2) / cos(w));
        k = atan2(-value(1, 3) / cos(w), -value(2, 3) / cos(w));
      % cos(elevation) == 0
      else
        k = 0; % (unconstrained)
        if value(3, 3) == 1
            w = pi / 2;
            p = -k + atan2(-value(1, 2), value(1, 1));
        else
            w = -pi / 2;
            p = k + atan2(-value(1, 2), value(1, 1));
        end
      end
      cam.viewdir = rad2deg([p w k]);
    end

    function cam = set.fullmodel(cam, value)
      if length(value) < 20
        error('Camera.fullmodel must have 20 elements.')
      end
      cam.xyz = value(1:3);
      cam.imgsz = value(4:5); cam.viewdir = value(6:8); cam.f = value(9:10);
      cam.c = value(11:12); cam.k = value(13:18); cam.p = value(19:20);
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

    function value = get.fmm(cam)
      if isempty(cam.sensorsz)
        error('Camera sensor size not set.');
      else
        value = cam.f .* cam.sensorsz ./ cam.imgsz;
      end
    end

    function value = get.fullmodel(cam)
      fullmodel = [cam.xyz, cam.imgsz, cam.viewdir, cam.f, cam.c, cam.k, cam.p];
      if length(fullmodel) == 20
        value = [cam.xyz, cam.imgsz, cam.viewdir, cam.f, cam.c, cam.k, cam.p];
      else
        error('Camera parameters missing or not of the correct length.')
      end
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
      % this function may be non-reversible.
      %
      % Inputs:
      %   scale - Scale factor [s] or desired image size [nx|ncols, ny|nrows]

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

    function cam = idealize(cam)
      % Make a copy
      cam = cam;
      % Strip distortion
      cam.k = zeros(size(cam.k));
      cam.p = zeros(size(cam.p));
      cam.c = [];
      cam.imgsz = cam.imgsz;
    end

    % Camera geometry

    function value = get.framebox(cam)
      value = [0 cam.imgsz(1) 0 cam.imgsz(2)] + 0.5;
    end

    function value = get.framepoly(cam)
      value = flipud([0 0 ; cam.imgsz(1) 0 ; cam.imgsz(1) cam.imgsz(2) ; 0 cam.imgsz(2) ; 0 0] + 0.5); % flip so cw in typical xy-space
    end

    function value = viewbox(cam, radius)
      if nargin < 2
        pyramid = cam.viewpyramid();
      else
        pyramid = cam.viewpyramid(radius);
      end
      value = [min(pyramid); max(pyramid)];
    end

    function value = viewpyramid(cam, radius)
      u = [0:cam.imgsz(1)]'; v = [0:cam.imgsz(2)]';
      edge_uv = [
        [u repmat(0, length(u), 1)];
        [repmat(cam.imgsz(1), length(v), 1) v];
        [flip(u) repmat(cam.imgsz(2), length(u), 1)];
        [repmat(0, length(v), 1) flip(v)]
      ];
      directions = cam.invproject(edge_uv);
      if nargin < 2
        radius = 1;
      end
      scaled_directions = bsxfun(@times, directions, radius ./ sum(directions.^2, 2));
      edge_xyz = bsxfun(@plus, scaled_directions, cam.xyz);
      value = [edge_xyz; cam.xyz];
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

    function h = plotDistortion(cam, scale)
      if nargin < 2
        scale = 1;
      end
      du = 100; nu = round(cam.imgsz(1) / du); du = cam.imgsz(2) / nu;
      u = (0:du:cam.imgsz(1)) + 0.5;
      v = (0:du:cam.imgsz(2)) + 0.5;
      [pu pv] = meshgrid(u, v);
      P0 = [pu(:) pv(:)]; P1 = camera2image(cam.idealize, cam.image2camera(P0));
      h = plot(cam.framepoly(:, 1), cam.framepoly(:, 2), 'k:'); hold on;
      quiver(P0(:, 1), P0(:, 2), scale * (P1(:, 1) - P0(:, 1)), scale * (P1(:, 2) - P0(:, 2)), 0, 'r');
      set(gca, 'xlim', cam.imgsz(1) * [-0.1 1.1], 'ylim', cam.imgsz(2) * [-0.1 1.1], 'ydir', 'reverse');
      axis equal; hold off;
    end

    % Transformations

    function in = inframe(cam, uv)
      box = cam.framebox;
      in = uv(:, 1) >= box(1) & uv(:, 1) <= box(2) & uv(:, 2) >= box(3) & uv(:, 2) <= box(4);
    end

    function in = infront(cam, xyz)
      xyz = bsxfun(@minus, xyz, cam.xyz);
      xyz = xyz * cam.R';
      in = xyz(:, 3) > 0;
    end

    function xy = image2camera(cam, uv)
      % Convert image to camera coordinates
      % [x, y] = [(x' - cx) / fx, (y' - cy) / fy]
      xy = [(uv(:, 1) - cam.c(1)) / cam.f(1), (uv(:, 2) - cam.c(2)) / cam.f(2)];
      % Remove distortion
      xy = cam.undistort(xy);
    end

    function dxyz = camera2world(cam, xy)
      % Convert camera coordinates to world ray directions
      dxyz = [xy ones(size(xy, 1), 1)] * cam.R;
    end

    function [xy, infront] = world2camera(cam, xyz, ray_directions)
      if nargin < 3
        ray_directions = false;
      end
      if ~ray_directions
        % Convert world to camera coordinates
        xyz = bsxfun(@minus, xyz, cam.xyz);
      end
      xyz = xyz * cam.R';
      % Normalize by perspective division
      xy = bsxfun(@rdivide, xyz(:, 1:2), xyz(:, 3));
      % Convert points behind camera to NaN
      infront = xyz(:, 3) > 0;
      xy(~infront, :) = NaN;
    end

    function uv = camera2image(cam, xy)
      xy = cam.distort(xy);
      uv = [cam.f(1) * xy(:, 1) + cam.c(1), cam.f(2) * xy(:, 2) + cam.c(2)];
    end

    function [uv, infront] = project(cam, xyz, ray_directions)
      % PROJECT  Project 3D world coordinates to 2D image coordinates.
      %
      %   uv = cam.project(xyz)
      %
      % Inputs:
      %   xyz - World coordinates [x1 y1 z1; x2 y2 z2; ...]
      %
      % Outputs:
      %   uv      - Image coordinates [u1 v1; u2 v2; ...]
      %
      % See also: camera.invproject

      if nargin < 3
        ray_directions = false;
      end
      [xy, infront] = cam.world2camera(xyz, ray_directions);
      uv = cam.camera2image(xy);
    end

    function xyz = invproject(cam, uv, S)
      % INVPROJECT  Project 2D image coordinates out as 3D world coordinates.
      %
      %   xyz = cam.invproject(uv[, S])
      %
      % Image coordinates are projected out of the camera as rays.
      % If a surface (S) is specified, the intersections with the surface are returned (or NaN if none).
      % Otherwise, ray directions are returned.
      %
      % Inputs:
      %   uv - Image coordinates [u1 v1; u2 v2; ...]
      %   S  - Surface, either as a DEM object or an infinite plane
      %        (defined as [a b c d], where ax + by + cz + d = 0)
      %
      % Outputs:
      %   xyz - World coordinates of intersections [x1 y1 z1; x2 y2 z2; ...], or
      %         ray directions [dx1 dy1 dz1; dx2 dy2 dz2; ...]
      %
      % See also: camera.project

      % Convert to normalized camera coordinates
      xy = cam.image2camera(uv);
      % Convert to ray direction vectors
      xyz = cam.camera2world(xy);
      % Track valid coordinates
      is_valid = ~any(isnan(xyz), 2);

      if nargin == 2
        % No surface: Return ray directions

      elseif nargin > 2 && isnumeric(S) && length(S) == 4
        % Plane: Return intersection of rays with plane
        xyz(is_valid, :) = intersectRayPlane(cam.xyz, xyz(is_valid, :), S);

      elseif nargin > 2 && strcmp(class(S), 'DEM')
        % DEM: Return intersection of rays with DEM
        for i = find(is_valid)'
          % xyz(i, :) = intersectRayDEM([cam.xyz xyz(i, :)], S);
          xyz(i, :) = S.sample(cam.xyz, xyz(i, :));
        end
      end
    end

    % Calibration

    function [newcam, fit] = optimize(cam, xyz, uv, freeparams)
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
      %   freeparams - Either a string, array, or 20-element vector describing
      %                which parameters should be optimized (see Examples).
      %
      % Outputs:
      %   newcam - Optimized camera
      %   rmse   - Root-mean-square reprojection error
      %   aic    - Akaike information criterion for reprojection errors, which
      %            can help determine an appropriate degree of complexity for
      %            the camera model (i.e. avoid overfitting).
      %            NOTE: Only strictly applicable for unweighted fitting.
      %
      % Examples:
      %   % Optimize all elements of viewdir:
      %   cam.optimizecam(xyz, uv, '00000111000000000000')
      %   cam.optimizecam(xyz, uv, 'viewdir')
      %   cam.optimizecam(xyz, uv, {'viewdir'})
      %   % Also optimize the third (z) element of xyz:
      %   cam.optimizecam(xyz, uv, '00100111000000000000')
      %   cam.optimizecam(xyz, uv, {'viewdir', 'xyz', 1})

      [newcam, fit] = camera.optimizeCams(cam, xyz, uv, freeparams);
      newcam = newcam{1};
    end

    function [newcam, fit] = optimizeR(cam, uv, uv2, angles)
      % Enforce defaults
      if nargin < 4
        angles = 1:3;
      end
      % Convert to ray directions
      dxyz = cam.invproject(uv);
      % Set up error function
      flexparams = camera.parseFreeparams({'viewdir', angles});
      fixparams = camera.parseFreeparams({});
      ray_direction = true;
      ef = @(d) reshape(projerror(camera.updateCams(d, {cam}, {flexparams}, fixparams){1}, dxyz, uv2, ray_direction), [], 1);
      % Optimize
      [d, rss] = LMFnlsq(ef, zeros(length(angles), 1));
      % Compile results
      newcam = camera.updateCams(d, {cam}, {flexparams}, fixparams){1};
      k = length(angles);
      n = size(uv, 1);
      fit = struct();
      fit.rmse = sqrt(rss ./ n);
      % AIC: https://en.wikipedia.org/wiki/Akaike_information_criterion
      % AIC small sample correction: http://brianomeara.info/tutorials/aic/
      % Camera model selectio: http://imaging.utk.edu/publications/papers/2007/ICIP07_vo.pdf
      fit.aic = n .* log(rss ./ n) + 2 * k * (n / (n - k - 1));
      fit.bic = n .* log(rss ./ n) + 2 * k * log(n);
      % fit.mdl = n .* log(rss ./ n) + 1 / (2 * k * log(n));
    end

    function e = projerror(cam, xyz, uv, ray_directions, lxyz, luv)
      if nargin < 4
        ray_directions = false;
      end
      e = [];
      if ~isempty(xyz) && ~isempty(uv)
        if size(xyz, 2) == 3
          puv = cam.project(xyz, ray_directions);
        elseif size(xyz, 2) == 2
          puv = cam.camera2image(xyz);
        end
        e = puv - uv(:, 1:2);
        if size(uv, 2) == 3
          e = e .* uv(:, [3 3]);
        end
      end
      if nargin >= 6 && ~isempty(lxyz) && ~isempty(luv)
        pluv = cam.project(lxyz);
        pos = polylinePoint(pluv, projPointOnPolyline(luv, pluv));
        % e = [sqrt(sum(e.^2, 2)); distancePointPolyline(luv, pluv)];
        e = [e ; pos - luv];
      end
    end

    function [X, edge] = horizon(cam, dem, ddeg)
      if nargin < 3
        ddeg = 1;
      end
      viewedges = cam.viewpyramid();
      dxy = bsxfun(@minus, viewedges(1:(end - 1), 1:2), viewedges(end, 1:2));
      angles = atan2d(dxy(:, 2), dxy(:, 1));
      ray_angles = [min(angles):ddeg:max(angles)]';
      [X, edge] = dem.horizon(cam.xyz, ray_angles);
    end

  end % methods

  methods (Static)

    function freeparams = parseFreeparams(freeparams)
      % Convert freeparams to boolean fullmodel selector
      if ischar(freeparams)
        % Convert to single-cell array
        freeparams = {freeparams};
      end
      if iscell(freeparams)
        % Parse cell array (given as name-element pairs)
        % e.g. {'viewdir', 'xyz', 3} => All viewdir elements and 3rd xyz element
        params = {'xyz', 'imgsz', 'viewdir', 'f', 'c', 'k', 'p'};
        param_indices = {[1:3], [4:5], [6:8], [9:10], [11:12], [13:18], [19:20]};
        is_param = cellfun(@ischar, freeparams);
        is_pos = cellfun(@isnumeric, freeparams);
        is_free = false(max(cell2mat(param_indices)), 1);
        for i = find(is_param)
          is_match = strcmpi(freeparams{i}, params);
          if any(is_match)
            ind = param_indices{is_match};
            if i < length(freeparams) && is_pos(i + 1)
              ind = ind(freeparams{i + 1});
            end
            is_free(ind) = true;
          end
        end
        freeparams = is_free;
      else
        % Convert numeric or character string to logical
        freeparams = ~(freeparams(:) == 0 | freeparams(:) == '0')';
      end
    end

    function [newcams, fit] = optimizeCams(cams, xyz, uv, flexparams, fixparams, lxyz, luv)

      if nargin < 7
        lxyz = [];
        luv = [];
      end

      % Convert inputs to cell arrays
      if ~iscell(cams), cams = {cams}; end
      if ~iscell(xyz), xyz = {xyz}; end
      if ~iscell(uv), uv = {uv}; end
      if ~iscell(lxyz), lxyz = {lxyz}; end
      if ~iscell(luv), luv = {luv}; end
      % Expand inputs
      n_cams = length(cams);
      if ~rem(n_cams, length(xyz)), xyz = repmat(xyz, 1, n_cams / length(xyz)); end
      if ~rem(n_cams, length(uv)), uv = repmat(uv, 1, n_cams / length(uv)); end
      if ~rem(n_cams, length(lxyz)), lxyz = repmat(lxyz, 1, n_cams / length(lxyz)); end
      if ~rem(n_cams, length(luv)), luv = repmat(luv, 1, n_cams / length(luv)); end
      if ~rem(n_cams, length(flexparams)), flexparams = repmat(flexparams, 1, n_cams / length(flexparams)); end
      if any(n_cams ~= [length(xyz), length(uv), length(flexparams), length(lxyz), length(luv)])
        error('Input arrays cannot be coerced to equal length')
      end
      % Enforce defaults
      if nargin < 5
        fixparams = {};
      end

      % Discard invalid points
      for i = 1:length(xyz)
        is_valid = cams{i}.infront(xyz{i}) & ~(any(isnan(xyz{i}), 2) | any(isnan(uv{i}), 2));
        xyz{i} = xyz{i}(is_valid, :);
        uv{i} = uv{i}(is_valid, :);
        if size(xyz{i}, 1) == 0
          error(['No valid control points found for camera ' str2num(i)]);
        end
      end

      % Initialize
      flexparams = cellfun(@camera.parseFreeparams, flexparams, 'UniformOutput', false);
      fixparams = camera.parseFreeparams(fixparams);
      n_flex = cellfun(@sum, flexparams);
      n_fix = sum(fixparams);

      % Optimize
      ef = @(m) reshape(cell2mat(cellfun(@projerror, camera.updateCams(m, cams, flexparams, fixparams), xyz, uv, repmat({false}, length(cams), 1), lxyz, luv, 'UniformOutput', false)), [], 1);
      m0 = zeros(n_fix + sum(n_flex), 1);
      mbest = LMFnlsq(ef, m0);

      % Compile results
      newcams = camera.updateCams(mbest, cams, flexparams, fixparams);
      e = cellfun(@projerror, newcams, xyz, uv, repmat({false}, length(cams), 1), lxyz, luv, 'UniformOutput', false);
      rss = cellfun(@(x) sum(sum(x.^2, 2)), e);
      n = cellfun(@(x, y) size(x, 1) + size(y, 1), uv, luv);
      fit = struct();
      fit.rmse = sqrt(rss ./ n);
      n = sum(n);
      k = n_fix + sum(n_flex);
      % AIC: https://en.wikipedia.org/wiki/Akaike_information_criterion
      % AIC small sample correction: http://brianomeara.info/tutorials/aic/
      % Camera model selectio: http://imaging.utk.edu/publications/papers/2007/ICIP07_vo.pdf
      fit.aic = n .* log(rss ./ n) + 2 * k .* (n ./ (n - k - 1));
      fit.bic = n .* log(rss ./ n) + 2 * k .* log(n);
      % fit.mdl = n .* log(rss ./ n) + 1 ./ (2 * k .* log(n));
    end

    function cams = updateCams(m, cams, flexparams, fixparams)
      % Initialize
      n_cams = length(cams);
      n_flex = cellfun(@sum, flexparams);
      n_fix = sum(fixparams);
      % m -> fullmodel for each camera
      temp = zeros(1, 20);
      temp(fixparams) = m(1:n_fix);
      m(1:n_fix) = [];
      for i = 1:n_cams
        d = temp;
        d(flexparams{i}) = m(1:n_flex(i));
        cams{i}.fullmodel = cams{i}.fullmodel + d;
        m(1:n_flex(i)) = [];
      end
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

      if any([cam.k, cam.p])
        % r = sqrt(x^2 + y^2)
        r2 = sum(xy.^2, 2);
        if any(cam.k)
          % Radial lens distortion
          % dr = (1 + k1 * r^2 + k2 * r^4 + k3 * r^6) / (1 + k4 * r^2 + k5 * r^4 + k6 * r^6)
          dr = 1 + cam.k(1) * r2 + cam.k(2) * r2.^2 + cam.k(3) * r2.^3;
          if any(cam.k(4:6))
            dr = dr ./ (1 + cam.k(4) * r2 + cam.k(5) * r2.^2 + cam.k(6) * r2.^3);
          end
        end
        if any(cam.p)
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
        if any(cam.k)
          xy = bsxfun(@times, xy, dr);
        end
        if any(cam.p)
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

      if any([cam.k, cam.p])

        % May fail for large negative k1.
        if cam.k(1) < -0.5
          warning(['Large, negative k1 (', num2str(cam.k(1), 3), '). Undistort may fail.'])
        end

        % If only k1 is nonzero, use closed form solution.
        % Cubic roots solution from Numerical Recipes in C 2nd Edition:
        % http://apps.nrbook.com/c/index.html (pages 183-185)
        if sum([cam.k cam.p] ~= 0) == 1 && cam.k(1)
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
            if any(cam.k)
              % Radial lens distortion
              % dr = (1 + k1 * r^2 + k2 * r^4 + k3 * r^6) / (1 + k4 * r^2 + k5 * r^4 + k6 * r^6)
              dr = 1 + cam.k(1) * r2 + cam.k(2) * r2.^2 + cam.k(3) * r2.^3;
              if any(cam.k(4:6))
                dr = dr ./ (1 + cam.k(4) * r2 + cam.k(5) * r2.^2 + cam.k(6) * r2.^3);
              end
            end
            if any(cam.p)
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
            if any(cam.p)
              xy = xyi - [dtx dty];
              if any(cam.k)
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
