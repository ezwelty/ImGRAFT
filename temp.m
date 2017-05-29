%% Sample DEM at xy points
dem = DEM(peaks(100)).build();
n = 1e3;
xy = [rand(n)' rand(n)'] * 100;
% Grid interpolation
z_grd = dem.sample_points(xy);
% Intersect with triangular mesh
z_tri = dem.sample_points_tri(xy);
% Results
nanstd(z_grd - z_tri)
figure();
dem.plot(3), hold on
scatter3(xy(:, 1), xy(:, 2), z_grd, 'r.')
scatter3(xy(:, 1), xy(:, 2), z_tri, 'y.')

%% Intersect DEM
dem = DEM(peaks(1e3), [0, 100], [100, 0]).build();
figure();
dem.plot(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
hold on
daspect([1, 1, 0.25]);
box = [sort(dem.xlim), sort(dem.ylim), sort(dem.zlim)];
drawBox3d(box);
origin = [30, -10, 4];
direction = [0.1, 1, -0.05];
scatter3(origin(1), origin(2), origin(3), 'b*')
quiver3(origin(1), origin(2), origin(3), direction(1), direction(2), direction(3), 120, 'k');
% Intersect bounding box
[tmin, tmax] = intersectRayBox(origin, direction, box([1, 3, 5]), box([2, 4, 6]));
start = origin + tmin * direction;
stop = origin + tmax * direction;
scatter3(start(:, 1), start(:, 2), start(:, 3), 'k*')
scatter3(stop(:, 1), stop(:, 2), stop(:, 3), 'k*')
% Intersect DEM (full search)
X = dem.sample_ray_tri(origin, direction, false);
if size(X, 1) > 0
  scatter3(X(:, 1), X(:, 2), X(:, 3), 'r*');
end
% Intersect DEM (bounded search)
for i = 1:size(X, 1)
  x = dem.sample_ray_tri(origin, direction, true, X(i, :), 1);
  scatter3(x(1), x(2), x(3), 'g*');
end
