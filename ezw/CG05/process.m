% CG05_20050811_130000
matlab -display :0.0 -nodesktop
cd ~/sites/ImGRAFT/
addpath(genpath('.'))

% Convert angles
R = [0.06991 0.99730 -0.02237; 0.03493 -0.02486 -0.99908; -0.99694 0.06907 -0.03657];
oblang = rot2oblang(R);
viewdir = [-(oblang(1) - 90) -oblang(2) oblang(3)] * pi / 180;
xyz = [496990.581 6776004.758 152.039];
imgsz = [2848 4288];
f = [5065.99156];

% Initialize anchor camera
cam = camera(xyz, imgsz, viewdir, f);
images = loadimages('ezw/CG05/CG05_20050811_130000.JPG', cam);
cg05 = load('ezw/CG05/CG05_20050811_130000.mat');
gcp = cg05.cam.gcp;
gcp(:, 1:2) = gcp(:, 1:2) + 0.5; % convert UL to Matlab image coordinates

% Optimize camera
[newcam, rmse, ~] = images(1).cam.optimizecam(gcp(:, 3:5), gcp(:, 1:2), '00000111111111000000')

% Plot reprojections vs original pixel coordinates

I = imread('ezw/CG05/CG05_20050811_130000.JPG');
% screen_size = get(0, 'ScreenSize');
% fig = figure();
% set(fig, 'Position', [0 0 screen_size(3) screen_size(4)]);

uv = cam.project(gcp(:, 3:5));
rmse = sqrt(sum(sum((gcp(:, 1:2) - uv).^2, 2)) / size(uv, 1));
figure
imshow(I / 1.5); hold on;
plot(gcp(:, 1), gcp(:, 2), 'yo', uv(:, 1), uv(:, 2), 'ro', [gcp(:, 1)'; uv(:, 1)'], [gcp(:, 2)'; uv(:, 2)'], 'y-');
% s = 1; quiver(gcp(:, 1), gcp(:, 2), s * (uv(:, 1) - gcp(:, 1)), s * (uv(:, 2) - gcp(:, 2)), 0, 'yellow')
title(sprintf('Reprojection errors (RMSE = %.1fpx)', rmse))

uv = newcam.project(gcp(:, 3:5));
rmse = sqrt(sum(sum((gcp(:, 1:2) - uv).^2, 2)) / size(uv, 1));
figure
imshow(I / 1.5); hold on;
plot(gcp(:, 1), gcp(:, 2), 'yo', uv(:, 1), uv(:, 2), 'ro', [gcp(:, 1)'; uv(:, 1)'], [gcp(:, 2)'; uv(:, 2)'], 'y-');
% s = 1; quiver(gcp(:, 1), gcp(:, 2), s * (uv(:, 1) - gcp(:, 1)), s * (uv(:, 2) - gcp(:, 2)), 0, 'yellow')
title(sprintf('Reprojection errors (RMSE = %.1fpx)', rmse))

%% Calibrate images from anchor
% TODO: Prepare best-of-class bedrock DEM
% TODO: Prepare glacier outlines
% TODO: Prepare glacier terminus heights
[Z, x, y] = geoimread('/volumes/science-b/data/columbia/dem/2004 Aerometric/20050811_2m.tif');
dem = DEM(Z, x, y);
cdem = dem.crop([4.9e5 5e5], [6.77e6 6.785e6], [0 Inf]);
smdem = cdem.resize(0.1);

% Viewpoint visibility
visible = voxelviewshed(smdem.X, smdem.Y, smdem.Z, cam.xyz);
% smdem.plot(3, visible); hold on;
% plot3(cam.xyz(1), cam.xyz(2), cam.xyz(3), 'y.', 'markersize', 40)

% Camera visibility
v = visible(:);
vxyz = [smdem.X(v) smdem.Y(v) smdem.Z(v)];
[vuv, ~, inframe] = cam.project(vxyz);
vuv = vuv(inframe, :);
vxyz = vxyz(inframe, :);

% Plot DEM-intersection errors
xyz2 = cam.invproject(vuv);
rays = [repmat(cam.xyz, size(xyz2, 1), 1) bsxfun(@minus, xyz2, cam.xyz)];
rvxyz = nan(size(rays, 1), 3);
for i = 1:size(rays, 1)
  rvxyz(i, :) = intersectRayDEM(rays(i, :), smdem);
end
plot(vxyz(:, 1), vxyz(:, 2), 'ko'); hold on
plot(rvxyz(:, 1), rvxyz(:, 2), 'r.');

figure
imshow(I / 1.5); hold on;
plot(vuv(:, 1), vuv(:, 2), 'g.');

%% Image segmentation

% Non-convex boundary polygon
vbi = boundary(vuv, 1);
vb = vuv(vbi, :);
figure
imshow(I / 1.5); hold on;
plot(vuv(:, 1), vuv(:, 2), 'g.');
plot(vb(:, 1), vb(:, 2), 'y-');

% Alpha shape (not better)
% vb = alphaShape(vuv, 'HoleThreshold', 1e9);
% [~, vbp] = boundaryFacets(vb);
% figure
% imshow(I / 1.5); hold on;
% plot(vb);
% plot(vbp(:, 1), vbp(:, 2), 'y-');
% plot(vuv(:, 1), vuv(:, 2), 'g.');

% ifSAR DEM
[Z, ~, bbox] = geotiffread('/volumes/science/data/columbia/_sort/ifsar/merged_projected_clipped.tif');
sardem = DEM(Z, bbox(:, 1), flip(bbox(:, 2)));
%bdem = bdem.crop([], [], [0 Inf]);
visible = voxelviewshed(bdem.X, bdem.Y, bdem.Z, cam.xyz);
bdem.plot(3); hold on
h = bdem.plot(3, visible * 255);
set(h, 'AlphaData', 0.2);
plot3(cam.xyz(1), cam.xyz(2), cam.xyz(3), 'r*', 'markersize', 20)

% Bed DEM
[Z, x, y] = geoimread('/volumes/science/data/columbia/bed/McNabb_100m.tif');
bdem = DEM(Z, x, y);
%bdem = bdem.crop([], [], [0 Inf]);
visible = voxelviewshed(bdem.X, bdem.Y, bdem.Z, cam.xyz);
bdem.plot(3); hold on
h = bdem.plot(3, visible * 255);
set(h, 'AlphaData', 0.2);
plot3(cam.xyz(1), cam.xyz(2), cam.xyz(3), 'r*', 'markersize', 20)

v = visible(:);
vbxyz = [bdem.X(v) bdem.Y(v) bdem.Z(v)];
[vbuv, ~, inframe] = cam.project(vbxyz);
vbuv = vbuv(inframe, :);
vbxyz = vbxyz(inframe, :);
figure
imshow(I / 1.5); hold on;
plot(vbuv(:, 1), vbuv(:, 2), 'y.', wuv(:, 1), wuv(:, 2), 'b-', vruv(:, 1), vruv(:, 2), 'g.')
plot(wuv(:, 1), wuv(:, 2), 'b-', vb(:, 1), vb(:, 2), 'k-', guv(:, 1), guv(:, 2), 'w-', vruv(:, 1), vruv(:, 2), 'g.')

% Bedrock DEM
[Z, x, y] = geoimread('/volumes/science/data/columbia/dem_analysis/2005_2010_mean.tif');
rdem = DEM(Z, x, y);
crdem = rdem.crop([4.9e5 5e5], [6.77e6 6.785e6], [30 Inf]);
smrdem = crdem.resize(0.1);
% Viewpoint visibility
visible = voxelviewshed(smrdem.X, smrdem.Y, smrdem.Z, cam.xyz);
vr = visible(:);
% Camera visibility
vrxyz = [smdem.X(vr) smdem.Y(vr) smdem.Z(vr)];
[vruv, ~, inframe] = cam.project(vrxyz);
vruv = vruv(inframe, :);
vrxyz = vrxyz(inframe, :);

figure
imshow(I / 1.5); hold on;
plot(wuv(:, 1), wuv(:, 2), 'b-', vb(:, 1), vb(:, 2), 'k-', guv(:, 1), guv(:, 2), 'w-', vruv(:, 1), vruv(:, 2), 'g.')


% Glacier outline
gshp = shaperead('/volumes/science/data/columbia/outline/2002_2010/20050811.shp');
z = interp2(smdem.X, smdem.Y, smdem.Z, gshp(1).X', gshp(1).Y');
gxyz = double([gshp(1).X' gshp(1).Y' z]); % HACK: double() needed for polybool??
guv = cam.project(gxyz);

% Forebay outline
bshp = shaperead('/volumes/science/data/columbia/outline/forebay.shp');
bxyz = [bshp(1).X' bshp(1).Y' repmat(0, length(bshp(1).X), 1)];
buv = cam.project(bxyz);

% Subtractions (bay - glacier)
[xd, yd] = polybool('subtraction', bxyz(:, 1), bxyz(:, 2), gxyz(:, 1), gxyz(:, 2));
wxy = [xd yd]; nans = find(isnan(xd)); wxy = wxy((nans(3) + 1):end, :);
plot(bxyz(:, 1), bxyz(:, 2), 'k-', gxyz(:, 1), gxyz(:, 2), 'k-', wxy(:, 1), wxy(:, 2), 'r-')
wuv = cam.project([wxy repmat(0, size(wxy, 1), 1)]);

figure
imshow(I / 1.5); hold on;
plot(wuv(:, 1), wuv(:, 2), 'b-', vb(:, 1), vb(:, 2), 'k-', guv(:, 1), guv(:, 2), 'w-')
plot(wuv(:, 1), wuv(:, 2), 'b-');



% plot(xyz(:, 1), xyz(:, 2), 'k-', newcam.xyz(1), newcam.xyz(2), 'ro');
% plot3(xyz(:, 1), xyz(:, 2), xyz(:, 3));

S = shaperead('/volumes/science/others/azinck/terminus.shp');
S = S(find([S.DATE]' == 20050811));
xyz = [S(1).X' S(1).Y' repmat(0, length(S(1).X), 1)];
uv = cam.project(xyz);
plot(uv(:, 1), uv(:, 2), 'g-');

%%%%

Ia = imread('ezw/CG05/CG05_20050811_130000.JPG');
Ib = imread('ezw/CG05/CG05_20050812_130000.JPG');

% Approximate estimate of image shift using a single large template
[du0, dv0] = templatematch(Ia, Ib, 1600, 800, 'templatewidth', 261, 'searchwidth', 400, 'supersample', 0.5)

% Gather image correspondences using a grid of probe points.
[pu, pv] = meshgrid(1000:50:1900, 500:50:1200);
[du, dv, C] = templatematch(Ia, Ib, pu, pv, 'templatewidth', 61, 'searchwidth', 81, 'supersample', 3, 'initialdu', du0, 'initialdv', dv0);
figure
imshow(Ia / 1.5); hold on;
s = 1; quiver(pu, pv, s * du, s * dv, 0, 'y')

% % Approximate estimate of image shift using a single large template
% [du0, dv0] = templatematch(Ia, Ib, 3800, 1200, 'templatewidth', 261, 'searchwidth', 400, 'supersample', 0.5)
%
% % Gather image correspondences using a grid of probe points.
% [pu, pv] = meshgrid(3400:10:4200, 1200:10:1400);
% [du, dv, C] = templatematch(Ia, Ib, pu, pv, 'templatewidth', 61, 'searchwidth', 81, 'supersample', 3, 'initialdu', du0, 'initialdv', dv0);
% figure
% imshow(Ia / 1.5); hold on;
% s = 1; quiver(pu, pv, s * du, s * dv, 0, 'y')

% Determine camera rotation from point correspondences
xyzi = cam.invproject([pu(:) pv(:)]);
rays = [repmat(cam.xyz, size(v1, 1), 1) bsxfun(@minus, xyzi, cam.xyz)];
% intersect with DEM
xyz = nan(size(rays, 1), 3);
for i = 1:size(rays, 1)
  xyz(i, :) = intersectRayDEM(rays(i, :), dem);
end

%xyz = cam.invproject([pu(:) pv(:)], dem);
[cam2, rmse, ~] = cam.optimizecam(xyz, [pu(:) + du(:), pv(:) + dv(:)], '00000111000000000000')
(cam2.viewdir - cam.viewdir) * 180 / pi
