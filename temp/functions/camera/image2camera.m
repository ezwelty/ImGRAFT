% IMAGE2CAMERA  Transform image into camera coordinates.
%
%   Xc = image2camera(Xi, camera)
%
% Compute the camera coordinates Xc for image coordinates Xi following
% transformation and distortion removal by the specified camera.
%
%
% Input:    Xi          2XN image coordinates
%           camera      camera structure                
%           
% Output:   Xc          3xN (normalized) camera coordinates
%
% See also image2world, world2camera, camera2image.

function Xc = image2camera(Xi, camera)

% Load camera
K = cameraMatrix(camera.fpx, camera.c);

% Convert to normalized image coordinates
Xi = inhomog(K \ homog(Xi));

% Undo lens distortion
Xi = undistort(Xi, camera.k);

% Convert to normalized camera coordinates
Xc = homog(Xi);