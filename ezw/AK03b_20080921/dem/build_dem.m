%--------------------------------------------------------------------------
% [dem] = build_dem(gridZ, xlim, ylim)
%
% Constructs a DEM structure from a grid of elevations (gridZ) and the
% outer x,y boundaries of the grid (xlim,ylim), with the following fields:
%
% z         grid of elevations (gridZ)
% minBound  minimum corner of bounding box [x y z]
% maxBound  maximum corner of bounding box [x y z]
% nx,ny     grid dimensions
% dx,dy     cell size
% cx,cy     cell center coordinates
%
% NOTE: Assumes that the grid of elevations is supplied such that +x points
% right (increasing column index), and +y points up (decreasing row index).
%
% 
% Input:    gridZ       grid of elevations
%           xlim        x boundaries of grid
%           ylim        y boundaries of grid
%           
% Output:   dem         resulting DEM structure         
%
%
% Author:   Ethan Welty
%
%--------------------------------------------------------------------------

function [dem] = build_dem(gridZ, xlim, ylim)

% Check number of arguments
if (nargin < 3)
    error('Specify grid of elevations (gridZ) and boundaries of the grid (xlim & ylim).')
end

% Test size of xlim
if (length(xlim) ~= 2)
    error('X bounds must be of length 2: [xmin xmax].')
end

% Test size of ylim
if (length(ylim) ~= 2)
    error('Y bounds must be of length 2: [ymin ymax].')
end


% Construct dem structure

% load elevations
dem.z = gridZ;

% calculate xyz bounding box
dem.minBound = [min(xlim) min(ylim) min(min(dem.z))];
dem.maxBound = [max(xlim) max(ylim) max(max(dem.z))];
 
% store DEM dimensions
dem.nx = size(dem.z,2);
dem.ny = size(dem.z,1);

% calculate cell dimensions
dem.dx = (dem.maxBound(1) - dem.minBound(1)) / dem.nx;
dem.dy = (dem.maxBound(2) - dem.minBound(2)) / dem.ny;

% construct matrices of cell center coordinates (grid)
%dem.cx = repmat((dem.minBound(1) + dem.dx/2) : dem.dx : (dem.maxBound(1) - dem.dx/2), dem.ny, 1);
%dem.cy = repmat(fliplr((dem.minBound(2) + dem.dy/2) : dem.dy : (dem.maxBound(2) - dem.dy/2))', 1, dem.nx);

% construct matrices of cell center coordinates (vectors)
dem.cx = ((dem.minBound(1) + dem.dx/2) : dem.dx : (dem.maxBound(1) - dem.dx/2))';
dem.cy = fliplr((dem.minBound(2) + dem.dy/2) : dem.dy : (dem.maxBound(2) - dem.dy/2))';