% CAMERA2IMAGE  Transform camera into image coordinates.
%
%   Xi = camera2image(Xc, camera)
%
% Compute the image coordinates Xi of camera coordinates Xc after
% distortion and transformation by the specified camera.
%
%
% Input:    Xc          3xN camera coordinates
%           camera      camera structure                  
%           
% Output:   Xi          2xN image coordinates
%
% See also image2world, world2camera, camera2image.

function Xi = camera2image(Xc, camera)

% Load camera
K = cameraMatrix(camera.fpx, camera.c);

% Normalize (just in case)
Xc = Xc(1:2,:) ./ (ones(2,1) * Xc(3,:));

% Apply distortion
Xc = distort(Xc, camera.k);

% convert to image coordinates
% p = [f * px + cx ; f * py + cy]
Xi = inhomog(K * homog(Xc));