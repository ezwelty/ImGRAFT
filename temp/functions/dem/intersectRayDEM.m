% INTERSECTRAYDEM  Find the intersection of a ray with a DEM.
%
%   intersect = intersectRayDEM(ray, dem)
%
% Traversal algorithm to find the intersection between a 3D ray and DEM,
% adjusted from the algorithm presented in "Using a DEM to determine
% geospatial object trajectories" by Collins et al (1998).
%
% NOTE: Currently intersects with 'side walls' rather than intersecting
% with an interpolated surface.
%
%
% Input:    ray         origin and direction of ray [x y z dx dy dz]
%           dem         DEM structure
%
% Output:   intersect   interpolated point of intersection [x y z]
%                       (-1 if no intersection)
%
% See also intersectRayBox, traverseRayDEM.

function intersect = intersectRayDEM(ray, dem)

% Check number of arguments
if (nargin < 2)
    error('Specify 2 input arguments.')
end

% Test size of arguments
if (length(ray) ~= 6)
    error('Ray must have length = 6.')
end

% Decompose input
origin = ray(1:3);
direction = ray(4:6);

% Cells
[cells, ~, stop] = traverseRayDEM(ray, dem);

% (flag values assigned, but currently unused)
if ~isempty(cells)

    % Convert to upper-left matrix indices (flip y)
    cells(:, 2) = dem.ny - (cells(:, 2) - 1);
    xi = cells(:, 1); yi = cells(:, 2);

    % Retrieve dem.Z based on cell xy
    ind = sub2ind(size(dem.Z), yi, xi);
    trueZ = dem.Z(ind);

    % locate exit points of ray
    cx = dem.x(xi(1:end - 1))'; cy = dem.y(yi(1:end - 1))';
    cx2 = dem.x(xi(2:end))'; cy2 = dem.y(yi(2:end))';
    dir = [-(cy2 - cy), cx2 - cx];
    point = [mean([cx cx2], 2), mean([cy cy2], 2)];
    exit = intersectLines([point dir], [origin(1:2) direction(1:2)]); % geom2d
    exit = [exit; stop(1:2)];

    % calculate Z of line at exit points
    lineZ = origin(3) + (exit(:, 1) - origin(1)) .* (direction(3) / direction(1));
    dZ = lineZ - trueZ;

    % locate intersected grid cell
    i = find(dZ <= 0,1);
    if (i == 1)
        pi = find(dZ >= 0,1);
        if (pi && pi ~= length(dZ))
            i = find(dZ(pi+1:end),1);
            i = i + pi;
        else
            i = [];
        end
    end

    % check for intersection at last cell
    % (tiny rounding error can lead to missed intersection)
    if isempty(i) && stop(3) == dem.min(3)
        i = size(cells,1);
        lineZ(i) = stop(3);
    end

    % interpolate intersection point
    if ~isempty(i) && i > 1

        % check for a "side impact"
        if lineZ(i-1) < trueZ(i)

            % intersect with vertical plane
            % FIXME: intersect with interpolated surface.
            vplane = [exit(i-1,1:2) 0 0 0 1 dir(i-1,1:2)/norm(dir(i-1,1:2)) 0];
            intersect = intersectLinePlane(ray, vplane); % geom2d
        else

            % linear interpolation
            trueX = exit(i,1) + (trueZ(i) - lineZ(i)) * (direction(1) / -abs(direction(3)));
            trueY = exit(i,2) + (trueZ(i) - lineZ(i)) * (direction(2) / -abs(direction(3)));
            intersect = [trueX trueY trueZ(i)];
        end
    else

        % return NaN
        intersect = -1;
    end

else

    % return NaN
    intersect = -1;
end
