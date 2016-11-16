% WORLD2IMAGE  Transform world into image coordinates.
%
%   Xc = world2image(X, camera)
%
% Compute the image coordinates Xi of world points X projected into the
% specified camera.
%
%
% Input:    X           3xN world coordinates
%           camera      camera structure
%
% Output:   Xi          2xN image coordinates
%
% See also image2world, world2camera, camera2image.

function Xi = world2image(X, camera)

Xc = world2camera(X, camera);
Xi = camera2image(Xc, camera);
