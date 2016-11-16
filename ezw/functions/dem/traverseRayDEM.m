% TRAVERSERAYDEM  Find the DEM grid cells traversed by a ray.
%
%   cells = traverseRayDEM(ray, dem)
%
% Traversal algorithm through a 2D grid, adapted from the 3D algorithm
% proposed by Amanatides and Woo (1987).
%
%
% Input:    ray         origin and direction of ray [x y z dx dy dz]
%           dem         DEM structure
%
% Output:   cells       DEM cells intersected by the ray
%           start       entrance point
%           stop        exit point
%
% See also intersectRayBox, intersectRayDEM.

% Adapted to 2D from code by Jesus Mena (MATLAB Central #26834).

function [cells, start, stop] = traverseRayDEM(ray, dem)

% Check number of arguments
if (nargin < 2)
    error('Specify two input arguments.')
end

% Test size of arguments
if (length(ray) ~= 6)
    error('Ray must have length = 6.')
end

% Test for ray - box intersection
[tmin, tmax] = intersectRayBox(ray, [dem.min dem.max]);

% Decompose input
origin = ray(1:3);
direction = ray(4:6);

% Calculate cells
if (~isempty(tmin))
    start = origin + tmin * direction;
    stop = origin + tmax * direction;

    boxSize = dem.max - dem.min;

    x = ceil( (start(1)-dem.min(1)) / dem.dx );
    y = ceil( (start(2)-dem.min(2)) / dem.dy );

    % bc: start == minBound
    if (x == 0) x = 1; end
    if (y == 0) y = 1; end

    mx = ceil( (stop(1)-dem.min(1)) / dem.dx );
    my = ceil( (stop(2)-dem.min(2)) / dem.dy );

    if (x==(dem.nx+1));  x=x-1;  end;
    if (y==(dem.ny+1));  y=y-1;  end;

    if (direction(1)>=0)
        tCellX = (x)/dem.nx;
        stepX = 1;
    else
        tCellX = (x-1)/dem.nx;
        stepX = -1;
    end;

    if (direction(2)>=0)
        tCellY = (y)/dem.ny;
        stepY = 1;
    else
        tCellY = (y-1)/dem.ny;
        stepY = -1;
    end;

    cellMaxX  = dem.min(1) + tCellX*boxSize(1);
    cellMaxY  = dem.min(2) + tCellY*boxSize(2);

    tMaxX      = tmin + (cellMaxX-start(1))/direction(1);
    tMaxY      = tmin + (cellMaxY-start(2))/direction(2);

    cellsizeX = boxSize(1)/dem.nx;
    cellsizeY = boxSize(2)/dem.ny;

    tDeltaX    = cellsizeX/abs(direction(1));
    tDeltaY    = cellsizeY/abs(direction(2));

    % Return list of traversed cells
    cells = [];
    while ((x <= dem.nx) && (x >= 1) && (y <= dem.ny) && (y >= 1))

        if (x == mx && y == my)
            cells(end+1,:) = [x y];
            break
        end

        % return cells
        cells(end+1,:) = [x y];

        if (tMaxX < tMaxY)
            x = x + stepX;
            tMaxX = tMaxX + tDeltaX;
        else
            y = y + stepY;
            tMaxY = tMaxY + tDeltaY;
        end
    end

else

    % Return empty
    cells = [];
    start = [];
    stop = [];
end
