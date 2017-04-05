% DEMPLOT  Plot a DEM in 2 or 3 dimensions.
%
%   h = plot_dem(dem, dim, scale = 1)
%
% Plots a DEM at the specified scale in either 2 or 3 dimensions.
%
% 
% Input:    dem         dem structure
%           dim         dimension of plot (2 or 3)
%           scale       decimal multiplier to resize dem before plotting
%           
% Output:   h           figure handle
%
% See also dembuild, demresize, demcrop.

function h = demplot(dem, dim, scale)

% Resize DEM
if (nargin > 2)
    dem = demresize(dem, scale);
end

% Force scaling of large DEM for 3D plotting
if (dem.nx * dem.ny > 1e6 && dim == 3)
    dem = demscale(dem, sqrt(1e6/(dem.nx*dem.ny)));
    disp('DEM was automatically scaled down for stable 3D plotting.')
end

% Plot DEM
h = figure;
if (dim == 3)
    shade = hillshade(dem.z,dem.cx,dem.cy);
    surf(dem.cx,dem.cy,dem.z,shade)
    shading interp
    axis equal
    colormap gray
elseif (dim == 2)
    shade = hillshade(dem.z,dem.cx,dem.cy);
    imagesc(dem.cx,dem.cy,shade)
    axis image
    set(gca,'ydir','normal')
    colormap gray
end