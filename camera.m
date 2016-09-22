%% camera class - a distorted camera model.
%
% This class is an implementation of a distorted camera model.
% Note: Uses <http://docs.opencv.org/modules/calib3d/doc/camera_calibration_and_3d_reconstruction.html the same distorted camera model as OpenCV>.
%
% camera Intrinsic parameters:
%  imgsz   - size of image in pixels [rows|height, columns|width] NOTE: y, x
%  f       - focal length in pixel units [f] or [fx, fy]
%  c       - camera center in pixel coordinates [cx, cy]
%  k       - radial distortion coefficients [k1...k6]
%  p       - tangential distortion coefficients [p1, p2]
%  xyz     - camera position [x, y, z]
%  viewdir - camera view direction in radians [yaw, pitch, roll]
%            yaw: counterclockwise rotation about z-axis (0 = look east)
%            pitch: rotation from horizon (+ look down, - look up)
%            roll: rotation about optical axis (+ down right, - down left, from behind)
%
% camera Derived parameters:
%  R         - rotation matrix calculated from camera view direction (read only)
%  K         - camera matrix [fx 0 cx; 0 fy cy; 0 0 1] (read only)
%  fullmodel - 20-element vector containing all camera parameters
%              [xyz(1:3), imgsz(1:2), viewdir(1:3), f(1:2), c(1:2), k(1:6), p(1:2)]
%
% camera Methods:
%  camera      - constructor
%  optimizecam - optimize camera parameters to mimimize pixel distances between projected world coordinates and associated image coordinates
%  project     - project world coordinates to image coordinates (3D -> 2D)
%  invproject  - project image coordinates to world coordinates (2D -> 3D)
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
% TODO: Reverse imgsz [y, x] to [x, y].

classdef camera

    properties
        xyz = [0 0 0];
        imgsz = [100 100];
        viewdir = [0 0 0];
        f = [5000 5000];
        c = [50 50];
        k = [0 0 0 0 0 0];
        p = [0 0];
    end

    properties (Dependent)
        R
        K
        fullmodel
    end

    methods

        function cam = camera(varargin)
            % CAMERA Construct a new camera object.
            %
            % There are three ways of calling this function –
            %
            % 1. Specify all camera parameters as a list:
            %
            %   cam = camera(xyz, imgsz, viewdir, f, c, k, p)
            %
            % 2. Specify all camera parameters as a 20-element (fullmodel) vector:
            %
            %   cam = camera([1 1 0 1024 768 pi 0 0 1000 1000 512 384 0 0 0 0 0 0 0 0])
            %
            % 3. Initialize with default parameters, then edit individual parameters:
            %
            %   cam = camera()
            %   cam.viewdir = [pi 0 0]; % look west

            % No arguments: return default camera
            if nargin == 0, return, end

            % Single argument: build camera from fullmodel vector
            if nargin == 1, cam.fullmodel = varargin{1}; end

            % Multiple arguments: set and validate camera parameters individually
            if nargin > 1

              % Set parameters
              cam.xyz = varargin{1};
              cam.imgsz = varargin{2};
              if nargin >= 3, cam.viewdir = varargin{3}; end
              if nargin >= 4, cam.f = varargin{4}; end
              if nargin >= 5, cam.c = varargin{5}; end
              if nargin >= 6, cam.k = varargin{6}; end
              if nargin >= 7, cam.p = varargin{7}; end

              % Validate parameters
              % xyz: 3-element vector, default = 0
              cam.xyz(end + 1:3) = 0; cam.xyz(4:end) = [];
              % imgsz: 2-element vector, no default
              if length(cam.imgsz) < 2, error('Image size (imgsz) must have two elements.'); end
              cam.imgsz(3:end)=[];
              % viewdir: 3-element vector, default = 0
              cam.viewdir(end + 1:3) = 0; cam.viewdir(4:end) = [];
              % f: 2-element vector, no default, expand [f] to [f, f]
              cam.f(end + 1:2) = cam.f(end); cam.f(3:end) = [];
              % c: 2-element vector, default = (imgsz[2 1] + 1) / 2
              if length(cam.c) < 2, error('Camera center (c) must have two elements.'); end
              if isempty(cam.c), cam.c = (cam.imgsz([2 1]) + 1)/2; end
              cam.c(3:end) = [];
              % k: 6-element vector, default = 0
              cam.k(end + 1:6) = 0; cam.k(7:end) = [];
              % p: 2-element vector, default = 0
              cam.p(end + 1:2) = 0; cam.p(3:end) = [];
            end
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
          % yaw: clockwise rotation about y-axis (relative to east, from above: + ccw, - cw) NOTE: - ccw, +cw
          %   ry = [C(1) 0 S(1); 0 1 0; -S(1) 0 C(1)];
          % pitch: clockwise rotation about x-axis (relative to horizon: - up, + down) NOTE: - ccw, +cw
          %   rp = [1 0 0; 0 C(2) -S(2); 0 S(2) C(2)];
          % roll: counterclockwise rotation about z-axis, (from behind camera: + ccw, - cw)
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
            if length(value) < 20, error('Camera model vector (fullmodel) must have 20 elements.'), end
            cam.xyz = value(1:3);
            cam.imgsz = value(4:5); cam.viewdir = value(6:8); cam.f = value(9:10);
            cam.c = value(11:12); cam.k = value(13:18); cam.p = value(19:20);
        end

        function value = get.K(cam)
          value = [cam.f(1) 0 cam.c(1); 0 cam.f(2) cam.c(2); 0 0 1];
        end

        function [uv, depth, inframe] = project(cam, xyz)
            % PROJECT Project 3D world coordinates into 2D image coordinates.
            %
            %   [uv, depth, inframe] = cam.project(xyz)
            %
            % Inputs:
            %    xyz – world coordinates [x1 y1 z1; x2 y2 z2; ...]
            %
            % Outputs:
            %    uv – image coordinates [u1 v1; u2 v2; ...]
            %    depth – distance of each point from the camera
            %    inframe – boolean whether each point is in the image

            if size(xyz, 2) > 3
                xyz = xyz';
            end

            % Convert to camera coordinates
            xyz = bsxfun(@minus, xyz, cam.xyz);
            xyz = xyz * cam.R';

            % Normalize by perspective division
            xy = bsxfun(@rdivide, xyz(:, 1:2), xyz(:, 3));

            % Apply lens distortion
            % TODO: Skip for xyz(:, 3) <= 0?
            if any([cam.k, cam.p] ~= 0)
              % r = sqrt(x^2 + y^2)
              r2 = sum(xy.^2, 2);
              if any(cam.k ~= 0)
                % Radial lens distortion
                % dr = (1 + k1 * r^2 + k2 * r^4 + k3 * r^6) / (1 + k4 * r^2 + k5 * r^4 + k6 * r^6)
                dr = 1 + (cam.k(1) * r2 + cam.k(2) * r2.^2 + cam.k(3) * r2.^3);
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

            % Convert to image coordinates
            % [x', y'] = [fx * x + cx, fy * y + cy]
            uv = [cam.f(1) * xy(:, 1) + cam.c(1), cam.f(2) * xy(:, 2) + cam.c(2)];

            % Set points behind camera to NaN
            uv(xyz(:, 3) <= 0, :) = NaN;

            if nargout > 1
              depth = xyz(:, 3);
            end
            if nargout > 2
              inframe = (depth > 0) & (uv(:, 1) >= 0.5) & (uv(:, 1) <= cam.imgsz(2) + 0.5) & (uv(:, 2) >= 0.5) & (uv(:, 2) <= cam.imgsz(1) + 0.5); % TODO: additional constraint for negative k1 and r2>1. (See orthorectification example)
            end
        end

        % TODO: Add plane intersection
        function xyz=invproject(cam,uv,X,Y,Z,xy0)
            % Inverse projection from 2d to 3d coordinates. (pixel->world)
            %
            % There are 3 ways to call this method:
            %
            % 1. Return a set of world coordinates which when projected
            %    results in the pixel coordinates (uv).
            %
            %   xyz=cam.invproject(uv)
            %
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
            nanix=any(isnan(uv),2);
            anynans=any(nanix);
            if anynans
                uv(nanix,:)=[];
            end
            if nargin==2

                %first an exact calculation based on non-distorted model...
                depth=1000;
                xyz=[(uv(:,1)-cam.c(1))/cam.f(1) (uv(:,2)-cam.c(2))/cam.f(2)]*depth;
                xyz(:,3)=depth;
                xyz=xyz*cam.R;
                xyz=bsxfun(@plus,xyz,cam.xyz);
                %then tune that so that it is consistent with distorted model
                if any(cam.k~=0)||any(cam.p~=0)
                    E=[1 0 0;0 1 0]*cam.R; % perturb in these directions.
                    for ii=1:size(uv,1)
                        nxyz=@(m)xyz(ii,:)+m'*E;
                        %misfit=@(m)sum((cam.project(nxyz(m))-uv(ii,:)).^2);
                        %m=fminsearch(misfit,[0 0]); %THIS REQUIRES A TOOLBOX!!!
                        misfitlm=@(m)(cam.project(nxyz(m))-uv(ii,:))'.^2;
                        m=LMFnlsq(misfitlm,[0;0]);
                        xyz(ii,:)=nxyz(m);
                    end
                end
            else
                visible=voxelviewshed(X,Y,Z,cam.xyz);
                Z=Z./visible;
                xyz=nan(size(uv,1),3);
                if nargin<6
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
                xyz(find(~nanix),:)=xyz; %find necessary because it ensures that xyz can grow.
                xyz(find(nanix),:)=nan;
            end
        end


        function [result,rmse,AIC]=optimizecam(cam,xyz,uv,freeparams)
            % Tune the camera so that projecting xyz results in uv (least squares)
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
end % classdef
