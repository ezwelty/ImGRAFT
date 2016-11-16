%--------------------------------------------------------------------------
% [rdem] = resize_dem(dem, scale)
%
% Resizes a DEM structure according to the specified scale.
%
% 
% Input:    dem         dem structure
%           scale       decimal multiplier
%           
% Output:   rdem        resized DEM structure         
%
%
% Author:   Ethan Welty
%
%--------------------------------------------------------------------------

function [rdem] = resize_dem(dem, scale)

% Scale grid
rgrid = imresize(dem.z, scale);

% Create resized DEM
xlim = [dem.minBound(1) dem.maxBound(1)];
ylim = [dem.minBound(2) dem.maxBound(2)];
rdem = build_dem(rgrid, xlim, ylim);