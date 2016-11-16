% FPX2FMM  Convert focal length from pixels to mm.
% 
%   fmm = fpx2fmm(fpx, wmm, imgsize)
%
% Calculates the focal length in mm (fmm), given the focal length in pixels
% (fpx), the longest dimension of the CCD or scanned film in mm (wmm) and
% the image resolution in pixels (imgsize):
%
% Requires that the image has not been cropped or warped beyond 
% simple resizing.
%
%
% Input:    fpx         focal length in pixels [fx (fy)]
%           wmm         longest dimension of CCD or film in mm
%           imgsize     image size in pixels [nx ny]
%
% Output:   fmm         focal length in mm
%                       (same dimensions as fmm)
%
% See also fmm2fpx, sensorWidth.

function fmm = fpx2fmm(fpx, wmm, imgsize)

% Check number of arguments
if (nargin < 3)
    error('Specify focal length, CCD or film width, and image size.')
end

% Calculate focal length in pixels
fmm = fpx * wmm / max(imgsize);