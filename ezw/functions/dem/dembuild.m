% DEMBUILD  Build a DEM structure from gridded values and grid boundaries.
%
%   dem = dembuild(zgrid, xlim, ylim)
%
% Constructs a DEM structure from a grid of elevations (zgrid) and the
% outer x,y boundaries of the grid (xlim,ylim), with the following fields:
%
%   z         grid of elevations (gridZ)
%   min       minimum corner of bounding box [x y z]
%   max       maximum corner of bounding box [x y z]
%   nx,ny     grid dimensions
%   dx,dy     cell size
%   cx,cy     cell center coordinate vectors
%
% Assumes that the grid of elevations is supplied such that +x points
% right (increasing column index), and +y points up (decreasing row index).
%
%
% Input:    zgrid       NxM matrix
%           xlim        limits of grid in x [min max]
%           ylim        limits of grid in y [min max]
%
% Output:   dem         DEM structure
%
% See also arc2dem.

function dem = dembuild(zgrid, xlim, ylim)

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

% Construct DEM structure

% load elevations
dem.z = zgrid;

% calculate xyz bounding box
dem.min = [min(xlim) min(ylim) min(min(dem.z))];
dem.max = [max(xlim) max(ylim) max(max(dem.z))];

% store DEM dimensions
dem.nx = size(dem.z,2);
dem.ny = size(dem.z,1);

% calculate cell dimensions
dem.dx = (dem.max(1) - dem.min(1)) / dem.nx;
dem.dy = (dem.max(2) - dem.min(2)) / dem.ny;

% construct matrices of cell center coordinates (vectors)
dem.cx = ((dem.min(1) + dem.dx/2) : dem.dx : (dem.max(1) - dem.dx/2));
dem.cy = fliplr((dem.min(2) + dem.dy/2) : dem.dy : (dem.max(2) - dem.dy/2))';

% construct matrices of cell center coordinates (grid)
dem.cx = repmat(dem.cx, dem.ny, 1);
dem.cy = repmat(dem.cy, 1, dem.nx);
