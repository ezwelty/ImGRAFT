%--------------------------------------------------------------------------
% [cdem] = clip_dem(gridZ, xlim, ylim)
%
% Clips a DEM structure to the specified X and Y coordinate boundaries.
%
% 
% Input:    dem         dem structure
%           xlim        x boundaries of grid
%           ylim        y boundaries of grid
%           
% Output:   cdem        clipped DEM structure         
%
%
% Author:   Ethan Welty
%
%--------------------------------------------------------------------------

function [cdem] = clip_dem(dem, xlim, ylim)

% Test size of xlim
if (length(xlim) ~= 2)
    error('X bounds must be of length 2: [xmin xmax].')
end

% Test size of ylim
if (nargin >=2 && length(ylim) ~= 2)
    error('Y bounds must be of length 2: [ymin ymax].')
end

% Convert x,y clipping bounds to grid indices
mincol = floor((min(xlim) - dem.minBound(1)) / dem.dx);
maxcol = ceil((max(xlim) - dem.minBound(1)) / dem.dx);
maxrow = ceil((dem.maxBound(2) - min(ylim)) / dem.dy);
minrow = floor((dem.maxBound(2) - max(ylim)) / dem.dy);

% Test for no overlap
if (mincol > dem.nx || maxcol < 1 || maxrow < 1 || minrow > dem.ny)
    error('Clipping bounds do not intersect DEM.')
end

% Intersect clipping bounds with dem extent
if (mincol < 1) mincol = 1; end
if (maxcol > dem.nx) maxcol = dem.nx; end
if (maxrow > dem.ny) maxrow = dem.ny; end
if (minrow < 1) minrow = 1; end

% Perform clip
cgrid = dem.z(minrow:maxrow, mincol:maxcol);
cxlim = [dem.cx(mincol)-dem.dx/2 , dem.cx(maxcol)+dem.dx/2]; 
cylim = [dem.cy(maxrow)-dem.dy/2 , dem.cy(minrow)+dem.dy/2];

% Create clipped DEM
cdem = build_dem(cgrid, cxlim, cylim);