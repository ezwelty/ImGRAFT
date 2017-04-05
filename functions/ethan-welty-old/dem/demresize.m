% DEMRESIZE  Resize a DEM.
%
%   rdem = demresize(dem, scale, method = 'bicubic')
%
% Resizes a DEM by the specified scale. Bicubic interpolation is used to
% resample z on the new grid.
%
% 
% Input:    dem         dem structure
%           scale       decimal multiplier
%           method      any method allowed by imresize
%           
% Output:   rdem        resized DEM structure         
%
% See also imresize, dembuild, demcrop, demplot.

function rdem = demresize(dem, scale, method)

if (nargin < 3)
    method = 'bicubic';
end

% Scale grid
rgrid = imresize(dem.z, scale, method);

% Create resized DEM
xlim = [dem.min(1) dem.max(1)];
ylim = [dem.min(2) dem.max(2)];
rdem = dembuild(rgrid, xlim, ylim);