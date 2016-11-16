% WORLD2CAMERA  Transform world into camera coordinates.
%
%   Xc = world2camera(X, camera)
%
% Compute the camera coordinates Xc of world points X projected into the
% specified camera.
%
%
% Input:    X           3xN world coordinates
%           camera      camera structure
%
% Output:   Xc          3xN (normalized) camera coordinates
%
% See also camera2world, camera2image.

function Xc = world2camera(X, camera)

% Load camera
t = -camera.R * camera.P';

% World to camera coordinates
% p = R*X + t
X = [camera.R t] * homog(X);

% Normalize by perspective division
% p = [px/pz ; py/pz]
X = X(1:3,:) ./ (ones(3,1) * X(3,:));

% normalized camera coordinates
Xc = X;
