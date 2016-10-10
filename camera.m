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
  % imgsz    - Size of image in pixels [ny|nrows|height, nx|ncols|width]
  % f        - Focal length in pixels [fx, fy]
  % c        - Camera center coordinates in pixels [cx, cy]
  % k        - Radial distortion coefficients [k1, ..., k6]
  % p        - Tangential distortion coefficients [p1, p2]
  % xyz      - Camera position in world coordinates [x, y, z]
  % viewdir  - Camera view direction in radians [yaw, pitch, roll]
  %            yaw: counterclockwise rotation about z-axis (0 = look east)
  %            pitch: rotation from horizon (+ look down, - look up)
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

  % TODO: Reverse directions of yaw and pitch angles.
  % TODO: Use degrees for yaw, pitch, and roll angles.
  % TODO: Reverse imgsz [ny, nx] to [nx, ny].
  % TODO: Split out DEM intersection function into DEM class method.
  % TODO: Limit voxelviewshed to camera view (not just position)
  % TODO: Remove imgsz from fullmodel (since not needed in optimization)
  % TODO: Clean up optimizecam. Add support for named arguments.
  % TODO: Extend optimizecam to work with camera bundle.

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
      p.addOptional('imgsz', [100 100], @(x) isnumeric(x) && length(x) == 2);
      p.addOptional('viewdir', [0 0 0], @(x) isnumeric(x) && length(x) <= 3);
      p.addOptional('f', [5000 5000], @(x) isnumeric(x) && length(x) <= 2);
      p.addOptional('c', [], @(x) isnumeric(x) && length(x) == 2);
      p.addOptional('k', [0 0 0 0 0 0], @(x) isnumeric(x) && length(x) <= 6);
      p.addOptional('p', [0 0], @(x) isnumeric(x) && length(x) <= 2);
      p.addOptional('sensorsz', [], @(x) isnumeric(x) && length(x) == 2);
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
      % xyz: 3-element vector, default = 0
      cam.xyz(end+1:3) = 0;
      % imgsz: 2-element vector, no default
      % viewdir: 3-element vector, default = 0
      cam.viewdir(end+1:3) = 0;
      % f: 2-element vector, no default, expand [f] to [f, f]
      if length(cam.f) == 1, cam.f(end + 1) = cam.f(end); end
      % c: 2-element vector, default = (imgsz[2 1] + 1) / 2
      if isempty(cam.c), cam.c = (cam.imgsz([2 1]) + 1) / 2; end
      % k: 6-element vector, default = 0
      cam.k(end+1:6) = 0;
      % p: 2-element vector, default = 0
      cam.p(end+1:2) = 0;
      % sensorsz: 2-element vector, no default
    end

    function value = get.R(cam)
      % Initial rotations of camera reference frame
      % (camera +z pointing up, with +x east and +y north)
      % Point camera north: -90 deg counterclockwise rotation about x-axis
      %   ri_1 = [1 0 0; 0 cosd(-90) sind(-90); 0 -sind(-90) cosd(-90)];
      % Point camera east: 90 deg counterclockwise rotation about y-axis
      %   ri_2 = [cosd(90) 0 -sind(90); 0 1 0; sind(90) 0 cosd(90)];
      % (camera +z now pointing east, with +x south and +y down)

      % View direction rotations
      C = cos(cam.viewdir); S = sin(cam.viewdir);
      % yaw: clockwise rotation about y-axis (relative to east, from above: + ccw, - cw)
      %   ry = [C(1) 0 S(1); 0 1 0; -S(1) 0 C(1)];
      % pitch: clockwise rotation about x-axis (relative to horizon: - up, + down)
      %   rp = [1 0 0; 0 C(2) -S(2); 0 S(2) C(2)];
      % roll: counterclockwise rotation about z-axis (from behind camera: + ccw, - cw)
      %   rr = [C(3) S(3) 0; -S(3) C(3) 0; 0 0 1];

      % Apply all rotations in order
      %   R = rr * rp * ry * ri_2 * ri_1;
      value = [ C(3) * S(1) - C(1) * S(2) * S(3), -C(1) * C(3) - S(1) * S(2) * S(3), -C(2) * S(3); ...
               -S(1) * S(3) - C(1) * C(3) * S(2),  C(1) * S(3) - C(3) * S(1) * S(2), -C(2) * C(3); ...
                C(1) * C(2)                     ,  C(2) * S(1)                     , -S(2)       ];
    end

    function value = get.fullmodel(cam)
      value = [cam.xyz, cam.imgsz, cam.viewdir, cam.f, cam.c, cam.k, cam.p];
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
        cam.f = value .* cam.imgsz([2 1]) ./ cam.sensorsz;
      end
    end

    function value = get.fmm(cam)
      if (isempty(cam.sensorsz))
        error('Camera sensor size not set.');
      else
        value = cam.f .* cam.sensorsz ./ cam.imgsz([2 1]);
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
      % this function may be nonreversible.
      %
      % Inputs:
      %   scale - Scale factor [s] or desired image size [ny|nrows, nx|ncols]

      imgsz0 = cam.imgsz;
      c0 = cam.c;
      if length(scale) > 1
        cam.imgsz = round(scale(1:2));
        scale = cam.imgsz ./ imgsz0;
        if any(abs(diff(scale)) * imgsz0 > 1)
          error(['Aspect ratio of target image size (' num2str(cam.imgsz(1) / cam.imgsz(2), 2) ') too different from original (' num2str(imgsz0(1) / imgsz0(2)) ') to be accounted for by rounding.']);
        end
      else
        cam.imgsz = round(scale * cam.imgsz);
        scale = cam.imgsz ./ imgsz0;
      end
      cam.f = cam.f .* flip(scale);
      cam.c = ((cam.imgsz([2 1]) + 1) / 2) + (c0 - ((imgsz0([2 1]) + 1) / 2)) .* flip(scale);
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
        inframe = (depth > 0) & (uv(:, 1) >= 0.5) & (uv(:, 1) <= cam.imgsz(2) + 0.5) & (uv(:, 2) >= 0.5) & (uv(:, 2) <= cam.imgsz(1) + 0.5);
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

        elseif nargin > 2 && isstruct(S)
          % DEM: Return intersection of rays with DEM
          [X, Y, Z] = deal(S.cx, S.cy, S.z);
          visible = voxelviewshed(X, Y, Z, cam.xyz);
          Z = Z ./ visible;
          xyz = nan(npts, 3);
          if nargin < 4
              [uv0,~,inframe]=cam.project([X(visible(:)),Y(visible(:)),Z(visible(:))]);
              uv0(:,3)=X(visible(:));
              uv0(:,4)=Y(visible(:));
              uv0(:,5)=Z(visible(:));
              uv0=uv0(inframe,:);
              if exist('scatteredInterpolant','file')>1
                  Xscat=scatteredInterpolant(uv0(:,3),uv0(:,4),uv0(:,3));
                  Xscat.Points=uv0(:,1:2);
                  Yscat=Xscat; Yscat.Values=uv0(:,4);
                  Zscat=Xscat; Zscat.Values=uv0(:,5);
              else
                  %fallback for older versions of matlab.
                  Xscat=TriScatteredInterp(uv0(:,3),uv0(:,4),uv0(:,3));  %#ok<REMFF1>
                  Xscat.X=uv0(:,1:2);
                  Yscat=Xscat; Yscat.V=uv0(:,4);
                  Zscat=Xscat; Zscat.V=uv0(:,5);
              end
              xy0=[Xscat(uv(:,1),uv(:,2)) Yscat(uv(:,1),uv(:,2)) Zscat(uv(:,1),uv(:,2))];
              xyz=xy0;

              if anynans
                  xyz(find(~nanix),:)=xyz; %find necessary because it ensures that xyz can grow.
                  xyz(find(nanix),:)=nan;
              end
              return
          end

          if Y(2,2)<Y(1,1)
              X=flipud(X);Y=flipud(Y);Z=flipud(Z);
          end
          if X(2,2)<X(1,1)
              X=fliplr(X);Y=fliplr(Y);Z=fliplr(Z);
          end

          if exist('griddedInterpolant','file')>1
              zfun=griddedInterpolant(X',Y',Z'); %TODO: improve robustness.
          else
              %fallback for older versions of matlab. slower
              zfun=@(x,y)interp2(X,Y,Z,x,y);
          end
          for ii=1:length(uv)
              %misfit=@(xy)sum((cam.project([xy zfun(xy(1),xy(2))])-uv(ii,1:2)).^2);
              misfitlm=@(xy)(cam.project([xy(:)' zfun(xy(1),xy(2))])-uv(ii,1:2))'.^2;
              try
                  %[xyz(ii,1:2),err]=fminunc(misfit,xy0(ii,1:2),optimset('LargeScale','off','Display','off','TolFun',0.001)); %TODO: remove dependency. can i use LMFnlsq?
                  xyz(ii,1:2)=LMFnlsq(misfitlm,xy0(ii,1:2));
                  xyz(ii,3)=zfun(xyz(ii,1),xyz(ii,2));
                  if sum(misfitlm(xyz(ii,1:2)))>2^2
                      xyz(ii,:)=nan; %do not accept greater than 2 pixel error.
                  end
              catch
              end
          end
        end

        if anynans
            xyz(find(~nanix), :) = xyz; % find() necessary because it ensures that xyz can grow.
            xyz(find(nanix), :) = NaN;
        end
    end


    function [result, rmse, AIC] = optimizecam(cam, xyz, uv, freeparams)
        % OPTIMIZECAM  Tune the camera so that projecting xyz results in uv (least squares)
        %
        % [newcamera,rmse,AIC]=cam.optimizecam(xyz,uv,freeparams)
        %
        %
        % If uv has three columns then the third column is interpreted as a
        % weight in the misfit function.
        %
        %
        % INPUTS:
        %    xyz: world coordinates.
        %    uv: target pixel coordinates.
        %        [optional 3rd column may specify weights]
        %    freeparams: a 20-element vector describing which camera
        %         parameters should be optimized. Follows same order as
        %         cam.fullmodel.
        %
        % OUTPUTS:
        %   newcamera: the optimized camera.
        %   rmse: root-mean-square-error
        %   aic: Akaike information criterion which can be used to
        %       help in determining an appropriate degree of complexity
        %       for the camera model (i.e. avoiding overfitting). [only
        %       strictly applicable for unweighted fitting]
        %
        % EXAMPLE:
        % %optimize the three view direction parameters.
        % %viewdir are the 6-8 columns in the <a href="matlab:help camera.fullmodel">camera.fullmodel</a>.
        %   [newcamera,rmse,AIC] = cam.optimizecam(xyz,uv,'00000111000000000000')
        %
        %

        nanrows=any(isnan(xyz),2)|any(isnan(uv),2);
        xyz(nanrows,:)=[];
        uv(nanrows,:)=[];

        fullmodel0=cam.fullmodel; %this describes the initial camera that is being perturbed.

        freeparams=~(freeparams(:)==0|freeparams(:)=='0')'; %convert to bools
        paramix=find(freeparams);
        Nfree=length(paramix);
        mbest=zeros(1,Nfree);

        newcam=@(m)camera( fullmodel0 + sparse(ones(1,Nfree),paramix,m,1,length(fullmodel0)) );

        if size(uv,2)==3
            misfit=@(m)reshape((project(newcam(m),xyz)-uv(:,1:2)).*uv(:,[3 3]),[],1);%weighted least squares
        else
            misfit=@(m)reshape(project(newcam(m),xyz)-uv,[],1);
        end
        if isnan(misfit(mbest))
            error('All GCPs must be infront of the initial camera location for optimizecam to work.'); %TODO: write better explanation. and remove requirement.
        end

        [mbest,RSS]=LMFnlsq(misfit,mbest); %WORKS SUPER FAST


        Nuv=size(uv,1);
        rmse=sqrt(RSS/Nuv);
        AIC=numel(uv)*log(RSS/numel(uv)) + 2*Nfree;
        result=newcam(mbest);
    end

  end % methods

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
