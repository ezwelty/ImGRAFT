% DISTORT  Apply radial and tangential lens distortion to normalized
% image coordinates.
%
%   Xd = distort(X, k)
%
%
% Input:    X           2XN undistorted normalized image coordinates
%           k           distortion parameters [k1 k2 p1 p2 k3]
%           
% Output:   Xd          2XN distorted normalized image coordinates
%
% See also world2image, camera2image.

function Xd = distort(X, k)

% Radial lens distortion
% r = sqrt( x^2 + y^2 )
% dr = 1 + k1 * r^2 + k2 * r^4 + k3 * r^6
r = sqrt(X(1,:).^2 + X(2,:).^2);
dr = 1 + ones(2,1) * (k(1) * r.^2 + k(2) * r.^4 + k(5) * r.^6);

% Tangential lens distortion
% dx = 2xy * p1 + p2 * (r^2 + 2x^2)
% dy = p1 * (r^2 + 2y^2) + 2xy * p2
dtx = 2 * X(1,:) .* X(2,:) * k(3) + k(4) * (r.^2 + 2 * X(1,:).^2);
dty = k(3) * (r.^2 + 2 * X(2,:).^2) + 2 * X(1,:) .* X(2,:) * k(4);

% Calculate distorted coordinates
Xd = X .* dr + [dtx ; dty];