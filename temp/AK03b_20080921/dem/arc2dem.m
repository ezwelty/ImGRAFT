%--------------------------------------------------------------------------
% [dem] = arc2dem(arcfile, nodata)
%
% Constructs a DEM structure from an ASCII file exported by ESRI ARC's
% Raster to ASCII tool. NoData grid cells are assigned NaN values. 
%
% See build_dem.m for details on the resulting DEM structure.
%
% 
% Input:    arcfile     path of the ascii dem file
%           nodata      nodata value, if different from file metadata
%           xlim        x clipping bounds [xmin xmax] (pending)
%           ylim        y clipping bounds [ymin ymax] (pending)
%           
% Output:   dem         resulting DEM structure         
%
%
% Author:   Ethan Welty
%
%--------------------------------------------------------------------------

function [dem] = arc2dem(arcfile, nodata)

% Check number of arguments
if (nargin < 1)
    error('Specify at least 1 input argument.')
end

% Check nodata is numeric
if (nargin >= 2)
    if (~isnumeric(nodata))
       error('NoData value must be numeric.')
    end
end

% Read metadata (example below)
% ncols     7905
% nrows     7438
% xllcorner 486426
% yllcorner 676636
% cellsize  4
% NODATA_value  -9999
fid = fopen(arcfile, 'r');
hlines = 0;
while fid
    
    line = fgetl(fid);
    tag = regexp(line, '([A-za-z_]*)', 'match');

    if (~isempty(tag))
        hlines = hlines + 1; 
        value = str2double(regexp(line, '([0-9\.\-\.]*)', 'match'));
        if (strcmp(tag,'ncols'))
            nx = value;
        elseif (strcmp(tag,'nrows'))
            ny = value;
        elseif (strcmp(tag,'xllcorner'))
            xmin = value;
        elseif (strcmp(tag,'yllcorner'))
            ymin = value;
        elseif (strcmp(tag,'cellsize'))
            dx = value;
            dy = value;
        elseif (strcmp(tag,'NODATA_value') && nargin < 2)
            nodata = value;
        end
    else
        break
    end
end
fclose(fid);

% Load raster grid
grid = txt2mat(arcfile, hlines, nx);

% Fill nodata cells with NaN
nullcells = grid == nodata;
grid(nullcells(:)) = NaN;

% Preprocess for build_dem
xlim = [xmin, xmin+dx*nx];
ylim = [ymin, ymin+dy*ny];

% Build DEM
dem = build_dem(grid, xlim, ylim);