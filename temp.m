dem = DEM(peaks(100));
n = 1e3;
xy = [rand(n)' rand(n)'] * 100;
z = dem.sample_points_tri(xy);

dem.plot(3), hold on
scatter3(xy(:, 1), xy(:, 2), z, 'r.')

dem = DEM(peaks(10));
box = [min([dem.zmin.xlim', dem.zmin.ylim', dem.zmin.zlim']), max([dem.zmin.xlim', dem.zmin.ylim', dem.zmin.zlim'])];
origin = [1 1 1];
direction = [-1 1 0];
ray = [origin, direction];
box2d = box([1, 4, 2, 5]);
ray2d = [origin(1:2), direction(1:2)];

[tmin, tmax] = intersectRayBox(origin, direction, box(1:3), box(4:6));
[flag, tmin, tmax] = rayBoxIntersection(origin, direction, box(1:3), box(4:6));

start = origin + tmin * direction;
stop = origin + tmax * direction;

% dem.plot(3), hold on, drawRay(ray2d)
% scatter3(start(:, 1), start(:, 2), start(:, 3), 'g*')
% scatter3(stop(:, 1), stop(:, 2), stop(:, 3), 'r*')

% clipRay(RAY, BOX)
% edgePosition
% intersectLineEdge

dem.sample_ray_tri(origin, direction)

img = Image('demos/data/images/AK03b_20080921/AK03b_20080621_121611.JPG');
file = img.file;
tic
for i = 1:1e4
  img.file;
end
toc
tic
for i = 1:1e4
  file;
end
toc
