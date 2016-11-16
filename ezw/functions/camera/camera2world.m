% CAMERA2WORLD  Transform camera coordinates to either world rays or
% points.
%
%   X = camera2world(Xc, camera)
%   X = camera2world(Xc, camera, surface)
%
% Camera coordinates Xc are transformed out the specified camera as rays.
% If a surface is specified, the intersections of the rays with the surface
% are returned (and NaN if intersection failed).
%
%
% Input:    Xc          3XN camera coordinates
%           camera      camera structure
%           surface     3D plane [x0 y0 z0 dx1 dy1 dz1 dx2 d2 dz2] or DEM
%
% Output:   X           3XN outgoing ray directions, or
%                       3XN points
%
% See also image2world, image2camera, camera2image, intersectRayDEM.

function X = camera2world(Xc, camera, surface)

% Normalize (just in case)
Xc = Xc(1:3,:) ./ (ones(3,1) * Xc(3,:));

% Ray direction
direction = camera.R' * Xc;

% Build rays
rays = [camera.P' * ones(1,size(direction,2)) ; direction]';

% Intersect with surface
if (nargin < 3)

    % no surface, return ray directions
    X = direction;

elseif isnumeric(surface) && length(surface) == 9

    % intersect with plane
    X = intersectLinePlane(rays, surface)';

elseif isstruct(surface)

    % intersect with DEM
    X = nan(3,size(Xc,2));
    for i = 1:size(Xc,2)
        X(:,i) = intersectRayDEM(rays(i,:), surface)';
    end
end
