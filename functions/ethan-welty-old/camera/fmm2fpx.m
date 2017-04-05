% FMM2FPX  Convert focal length from mm to pixels.
% 
%   fpx = fmm2fpx(fmm, wmm, imgsize)
%
% Calculates the focal length in pixels (fpx), given the focal length in mm
% (fmm), the longest dimension of the CCD or scanned film in mm (wmm) and
% the image resolution in pixels (imgsize).
%
% Requires that the image has not been cropped or warped beyond 
% simple resizing.
%
%
% Input:    fmm         focal length in mm [fx (fy)]
%           wmm         longest dimension of CCD or film in mm
%           imgsize     image size in pixels [nx ny]
%
% Output:   fpx         focal length in pixels 
%                       (same dimensions as fmm)
% 
% See also fpx2fmm, sensorWidth.

function fpx = fmm2fpx(fmm, wmm, imgsize)

% Check number of arguments
if (nargin < 3)
    error('Specify focal length, CCD or film width, and image size.')
end

% Calculate focal length in pixels
fpx = fmm * (max(imgsize) / wmm);