% DEMCROP  Crop a DEM.
%
%   cdem = demcrop(dem, xlim, ylim)
%
% Crop a DEM to the specified x and y boundaries. Includes cells
% intersected by the boundary (rather than resampling the DEM to conform
% perfectly to the crop boundaries).
%
% 
% Input:    dem         DEM structure
%           xlim        x crop boundaries
%           ylim        y crop boundaries
%           
% Output:   cdem        cropped DEM         
%
% See also dembuild, demresize, demplot.

function cdem = demcrop(dem, varargin)

% Convert x,y clipping bounds to grid indices
mincol = floor((min(xlim) - dem.min(1)) / dem.dx);
maxcol = ceil((max(xlim) - dem.min(1)) / dem.dx);
maxrow = ceil((dem.max(2) - min(ylim)) / dem.dy);
minrow = floor((dem.max(2) - max(ylim)) / dem.dy);

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
cdem = dembuild(cgrid, cxlim, cylim);