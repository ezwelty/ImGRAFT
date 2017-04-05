function [tmin, tmax] = intersectRayBox(ray, box)
  % INTERSECTRAYBOX  Find the intersections of a ray with a box.
  %
  %   [tmin, tmax] = intersectRayBox(ray, box)
  %
  % Finds the intersections of a ray with a box using the algorithm presented
  % in "Efficiency issues for ray tracing" by Smit (1998).
  %
  % Inputs:
  %   ray - Origin and direction of ray [x y z dx dy dz]
  %   box - Two opposite corners of box [x1 y1 z1 x2 y2 z2]
  %
  % Outputs:
  %   tmin - Distance of entrance point from ray origin (zero if origin inside box)
  %   tmax - Distance of exit point from ray origin
  %
  % Adapted from code by Jesus Mena (MATLAB Central #26834).

  % TODO: Support multiple rays.
  % TODO: Optimize for rays of shared origin.

  % Check number of arguments
  if (nargin < 2)
    error('Specify 2 input arguments.')
  end

  % Test size of arguments
  if (length(ray) ~= 6 || length(box) ~= 6)
    error('Both arguments must have length = 6.')
  end

  % Break down input
  origin = ray(1:3);
  direction = ray(4:6);
  boxmin = box(1:3);
  boxmax = box(4:6);

  % Find intersections
  if direction(1) >= 0
    tmin = (boxmin(1) - origin(1)) / direction(1);
    tmax = (boxmax(1) - origin(1)) / direction(1);
  else
    tmin = (boxmax(1) - origin(1)) / direction(1);
    tmax = (boxmin(1) - origin(1)) / direction(1);
  end

  if direction(2) >= 0
    tymin = (boxmin(2) - origin(2)) / direction(2);
    tymax = (boxmax(2) - origin(2)) / direction(2);
  else
    tymin = (boxmax(2) - origin(2)) / direction(2);
    tymax = (boxmin(2) - origin(2)) / direction(2);
  end

  if tmin > tymax || tymin > tmax
    tmin = [];
    tmax = [];
    return
  end

  if tymin > tmin
    tmin = tymin;
  end

  if tymax < tmax
    tmax = tymax;
  end

  if direction(3) >= 0
    tzmin = (boxmin(3) - origin(3)) / direction(3);
    tzmax = (boxmax(3) - origin(3)) / direction(3);
  else
    tzmin = (boxmax(3) - origin(3)) / direction(3);
    tzmax = (boxmin(3) - origin(3)) / direction(3);
  end

  if tmin > tzmax || tzmin > tmax
    tmin = [];
    tmax = [];
    return
  end

  if tzmin > tmin
    tmin = tzmin;
  end

  if tzmax < tmax
    tmax = tzmax;
  end

  % If origin is inside box, start ray at origin
  if tmin < 0
    tmin = 0;
  end
