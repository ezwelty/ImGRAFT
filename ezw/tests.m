tic
x = 1:1e7;
length(unique(diff(x))) == 1;
toc

tic
x = 1:1e7;
dx = diff(x);
any(dx(1) ~= dx);
toc


% Camera parameters
cam.xyz = [0 0 100];
cam.viewdir = [0 pi/2 0];

% Surface
plane = [0 0 1 0];
dem = dembuild(zeros(100), [-1000 1000], [-1000 1000]);

% Prepare points
[xg, yg] = meshgrid(1:500:cam.imgsz(2), 1:500:cam.imgsz(1));
uvi = [xg(:) yg(:)];

% Invproject / project
tic
xyzi = cam.invproject(uvi, plane);
toc
tic
uv = cam.project(xyzi);
toc

%sqrt(sum((uv - uvi).^2, 2))
sqrt(mean(sum((uv - uvi).^2, 2)))

% Plot reprojection errors
plot(uvi(:, 1), uvi(:, 2), 'ko'), hold on
plot(uv(:, 1), uv(:, 2), 'ro')



%%%
cam = camera();
cam.xyz = [446722 7396671 770];
cam.imgsz = [2856 4290];
cam.viewdir = [3.82594030371847 0.101305139574579 -0.0129904722239836];
cam.f = [6016.38351662639 6027.08180377813];
cam.c = [2145.5 1428.5];
cam.k = [-0.117189198645842 0 0 0 0 0];
cam.p = [0 0];

datafolder = 'demos';
dem = load(fullfile(datafolder,'dem')); %load DEM
nx = (dem.X(1, 2) - dem.X(1, 1));
ny = (dem.Y(end-1, 1) - dem.Y(end, 1));
dem = dembuild(dem.filled, [dem.X(1, 1) - nx/2, dem.X(1, end) + nx/2], [dem.Y(end, 1) - ny/2, dem.Y(1, 1) + ny/2]);

[x, y] = meshgrid(dem.cx(1, 1):50:dem.cx(1, end), dem.cy(end, 1):50:dem.cy(1, 1));
z = interp2(dem.cx, dem.cy, dem.z, x, y);
visible = voxelviewshed(dem.cx, dem.cy, dem.z, cam.xyz);
keepers = interp2(dem.cx, dem.cy, double(visible), x(:), y(:)) == 1;
xyzi = [x(keepers) y(keepers) z(keepers)];
[uv, ~, inframe] = cam.project(xyzi);
xyzi = xyzi(inframe, :);

uv = cam.project(xyzi);
xyz = cam.invproject(uv, dem);

sqrt(mean(sum((xyz - xyzi).^2, 2)))
plot(xyzi(:, 1), xyzi(:, 2), 'ko'), hold on
plot(xyz(:, 1), xyz(:, 2), 'ro')
