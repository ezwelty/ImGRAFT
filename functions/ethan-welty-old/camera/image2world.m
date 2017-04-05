% IMAGE2WORLD  Transform image coordinates to either world rays or
% points.
%
%   X = image2world(Xi, camera)
%   X = image2world(Xi, camera, surface)
%
% Image coordinates Xi are projected out the specified camera as rays. If a
% surface is specified, the intersections of the rays with the surface are
% returned (and NaN if intersection failed).
%
%
% Input:    Xi          2XN image coordinates
%           camera      camera structure
%           surface     3D plane [x0 y0 z0 dx1 dy1 dz1 dx2 d2 dz2] or DEM
%           
% Output:   X           3XN outgoing ray directions, or
%                       3XN points
%                       
% See also world2image, image2camera, camera2world, intersectRayDEM.

function X = image2world(Xi, camera, surface)

% Image to camera coordinates
Xc = image2camera(Xi, camera);

% Camera to world coordinates
if nargin < 3
    X = camera2world(Xc, camera);
else
    X = camera2world(Xc, camera, surface);
end
  
end