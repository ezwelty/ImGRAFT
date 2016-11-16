% IMGCENTER  Calculate the center coordinates of images.
%
%   c = imgcenter(imgsize)
%
% Calculates the center coordinates of images as [nx/2, ny/2].
%
%
% Input:    imgsize     Nx2 matrix of image sizes in pixels [nx ny]
%
% Output:   c           Nx2 matrix of image center coordinates [cx cy]
%
% See also cameraMatrix.

function c = imgcenter(imgsize)

% Test size of imgsize
if (size(imgsize,2) ~= 2)
    error('Specify an NX2 matrix of image dimensions [nx ny].')
end

% Calculate coordinates of image center
c = [ imgsize(:,1) / 2 , imgsize(:,2) / 2 ];