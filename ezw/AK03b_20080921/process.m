%%% AK03b_20080921

% matlab -display :0.0 -nodesktop
cd ~/sites/ImGRAFT/
addpath(genpath('.'))
root_dir = fullfile('ezw', 'AK03b_20080921');

%% Load images
c = camera('xyz', [497017.2113 6776153.168 163.52], 'viewdir', [-22.5 -5 3]);
img_files = dir(fullfile(root_dir, 'images', '*.JPG'));
img_paths = arrayfun(@(f) fullfile(f.folder, f.name), img_files, 'UniformOutput', false);
images = loadimages(img_paths, c);

%% Load control points
GCP = readtable(fullfile(root_dir, 'gcp.tsv'), 'filetype', 'text', 'delimiter', '\t');
GCP.Properties.RowNames = lower(GCP.name);
svg_files = dir(fullfile(root_dir, 'images', '*.svg'));
svg_paths = arrayfun(@(f) fullfile(f.folder, f.name), svg_files, 'UniformOutput', false);
svg = cellfun(@svg2struct, svg_paths, 'UniformOutput', false);
for i = 1:length(svg_paths)
  [folder, file, ext] = fileparts(svg_paths{i});
  has_svg = strcmp(fullfile(folder, [file '.JPG']), {images.path});
  images(has_svg).svg = svg{i};
  if isfield(svg{i}, 'gcp')
    gcp_names = fieldnames(svg{i}.gcp);
    gcp_xyz = GCP{lower(gcp_names), {'x_wgs84', 'y_wgs84', 'z_hae'}};
    for j = 1:length(gcp_names)
      images(has_svg).svg.gcp.(gcp_names{j}) = [images(has_svg).svg.gcp.(gcp_names{j}) + 0.5, gcp_xyz(j, :)];
    end
  end
end

%% Calibrate anchor images
% (images with svg.gcp field)
is_anchor = arrayfun(@(img) eval('isfield(img.svg, ''gcp'')'), images)';
is_anchor_ind = find(is_anchor);
cams = {images(is_anchor).cam};
gcp = arrayfun(@(img) cell2mat(struct2cell(img.svg.gcp)), images(is_anchor), 'UniformOutput', false);
uv = cellfun(@(x) x(:, 1:2), gcp, 'UniformOutput', false);
xyz = cellfun(@(x) x(:, 3:5), gcp, 'UniformOutput', false);
[newcams, fit] = camera.optimizeCams(cams, xyz, uv, {'viewdir'}, {'f'})
for i = 1:length(newcams)
  images(is_anchor_ind(i)).cam = newcams{i};
end

%% Assign images to nearest anchor
anchor_dates = [images(is_anchor).date];
for i = find(~is_anchor)
  [~, i_min] = min(abs(anchor_dates - images(i).date));
  images(i).anchor = images(is_anchor_ind(i_min)).path;
end

%% Segment anchor images
for i = find(is_anchor)
  % Water: image frame clipped by terminus (coastline trace)
  c = polyclip(images(i).cam.framepoly, images(i).svg.terminus + 0.5);
  waterpoly = c{1}{1}; % FIXME: Choose automatically
  % Glacier: glacier - water
  [gx, gy] = polybool('subtraction', images(i).svg.glacier(:, 1) + 0.5, images(i).svg.glacier(:, 2) + 0.5, waterpoly(:, 1), waterpoly(:, 2));
  glacierpoly = [gx, gy];
  % Sky: image frame clipped by horizon
  c = polyclip(images(i).cam.framepoly, images(i).svg.horizon + 0.5);
  skypoly = c{1}{2}; % FIXME: Choose automatically
  % Land: image frame - (sky + glacier + water)
  [lx, ly] = polybool('subtraction', images(i).cam.framepoly(:, 1), images(i).cam.framepoly(:, 2), {glacierpoly(:, 1), waterpoly(:, 1), skypoly(:, 1)}, {glacierpoly(:, 2), waterpoly(:, 2), skypoly(:, 2)});
  landpolys = cellfun(@(x, y) [x y], lx, ly, 'UniformOutput', false);
  % Save
  images(i).glacierpolys = {glacierpoly};
  images(i).landpolys = landpolys;
  % % Plot
  % figure()
  % imshow(imread(images(i).path));
  % patch(skypoly(:, 1), skypoly(:, 2), 1, 'Facecolor', 'b');
  % patch(glacierpoly(:, 1), glacierpoly(:, 2), 1, 'Facecolor', 'w');
  % patch(waterpoly(:, 1), waterpoly(:, 2), 1, 'Facecolor', 'b');
  % for j = 1:length(landpolys)
  %   patch(landpolys{j}(:, 1), landpolys{j}(:, 2), 1, 'Facecolor', 'y');
  % end
  % alpha(0.5);
end

%% Calibrate images from assigned anchors
image_paths = {images.path};
for i = find(~is_anchor)
  img = imread(images(i).path);
  anchor_i = find(strcmp(images(i).anchor, image_paths));
  anchor_img = imread(images(anchor_i).path);
  % Generate grid of points in land polygons
  gdu = 100; gdv = 100;
  mpu = []; mpv = []; mdu = []; mdv = [];
  for i_poly = 1:length(landpolys)
    u = 0:gdu:range(landpolys{i_poly}(:, 1)) + min(landpolys{i_poly}(:, 1));
    v = 0:gdv:range(landpolys{i_poly}(:, 2)) + min(landpolys{i_poly}(:, 2));
    [pu pv] = meshgrid(u, v);
    [in on] = inpolygon(pu, pv, landpolys{i_poly}(:, 1), landpolys{i_poly}(:, 2));
    pu = pu(in & ~on); pv = pv(in & ~on);
    [du, dv] = templatematch(Ia, Ib, pu, pv, 'templatewidth', gdu / 2, 'searchwidth', gdu);
    mpu = [mpu ; pu]; mpv = [mpv ; pv]; mdu = [mdu ; du]; mdv = [mdv ; dv];
  end
  figure
  imshow(Ia / 1.5); hold on;
  s = 10; quiver(mpu, mpv, s * mdu, s * mdv, 0, 'y')
  % Idealize (undistort) image coordinates
  nonans = ~(isnan(mdu) | isnan(mdv));
  mpu = mpu(nonans); mpv = mpv(nonans);
  mdu = mdu(nonans); mdv = mdv(nonans);
  A = cam.normalize([mpu mpv]);
  B = cam.normalize([mpu + mdu, mpv + mdv]);
  [F, inliersF] = ransacfitfundmatrix(A', B', 0.00000001); % tricky to set threshold (re-express in pixels?)
  [H, inliersH] = ransacfithomography(A', B', 0.000001); % tricky to set threshold (re-express in pixels?)
  figure
  imshow(Ia / 1.5); hold on;
  s = 10; quiver(mpu(inliers), mpv(inliers), s * mdu(inliers), s * mdv(inliers), 0, 'y')

end


%% Load DEM
[Z, ~, bbox] = geotiffread('/volumes/science/data/columbia/_new/ifsar/merged_projected_clipped.tif');
dem = DEM(Z, bbox(:, 1), flip(bbox(:, 2)));
viewbox = images(1).cam.viewbox(30 * 1e3);
cdem = dem.crop(viewbox(:, 1), viewbox(:, 2));
smdem = cdem.resize(0.2);

%% Segment anchor image
% Horizon polyline
hxyz = cam.horizon(smdem);
huv = cam.project(hxyz);
vxyz = [cam.xyz ; hxyz ; cam.xyz];
% Sky polygon
c = polyclip(cam.framepoly, huv);
skypoly = c{1}{2};
% Forebay polygon
bshp = shaperead('/volumes/science/data/columbia/outline/forebay.shp');
[bx, by] = polybool('intersection', bshp(1).X', bshp(1).Y', flip(vxyz(:, 1)), flip(vxyz(:, 2)));
bxyz = [bx by repmat(0, length(bx), 1)];
buv = cam.project(bxyz);
[bu bv] = polybool('intersection', cam.framepoly(:, 1), cam.framepoly(:, 2), flip(buv(:, 1)), flip(buv(:, 2)));
baypoly = [bu bv];
% Terminus polyline
% S = shaperead('/volumes/science/others/azinck/terminus.shp');
% S = S(find([S.DATE]' == 20080811));
% txyz = [S(1).X' S(1).Y' repmat(0, length(S(1).X), 1)];
tuv = svg.terminus.polyline + 0.5;
% Clip forebay with terminus
% c = polyclip(bxyz(:, 1:2), txyz(:, 1:2));
% wxyz = [c{1}{1} repmat(0, size(c{1}{1}, 1), 1)];
% wuv = cam.project(wxyz);
% waterpoly = clipPolygon(wuv, cam.framebox);
c = polyclip(baypoly, tuv);
waterpoly = flipud(c{1}{1});
% Glacier polygon
glacierpoly = svg.glacier.polygon + 0.5;
[gu, gv] = polybool('subtraction', glacierpoly(:, 1), glacierpoly(:, 2), waterpoly(:, 1), waterpoly(:, 2));
glacierpoly = [gu gv; gu(1) gv(1)];
% Land polygons
[Lx Ly] = polybool('subtraction', cam.framepoly(:, 1), cam.framepoly(:, 2), {glacierpoly(:, 1), waterpoly(:, 1), skypoly(:, 1)}, {glacierpoly(:, 2), waterpoly(:, 2), skypoly(:, 2)});
landpolys = {[Lx{1} Ly{1}], [Lx{2} Ly{2}]};
% Plot
imshow(I / 1.5); hold on;
patch(skypoly(:, 1), skypoly(:, 2), 1, 'Facecolor', 'b');
% patch(baypoly(:, 1), baypoly(:, 2), 1, 'Facecolor', 'b');
% plot(tuv(:, 1), tuv(:, 2), 'y')
patch(glacierpoly(:, 1), glacierpoly(:, 2), 1, 'Facecolor', 'w');
patch(waterpoly(:, 1), waterpoly(:, 2), 1, 'Facecolor', 'b');
patch(landpolys{1}(:, 1), landpolys{1}(:, 2), 1, 'Facecolor', 'y');
patch(landpolys{2}(:, 1), landpolys{2}(:, 2), 1, 'Facecolor', 'y');
alpha(0.5);

%% Load second image
images = loadimages('ezw/AK03b_20080921/images/AK03b_20080621_211611.JPG', cam);
Ia = imread(images(1).path);
images(2) = loadimages('ezw/AK03b_20080921/images/AK03b_20080622_071612.JPG', cam);
Ib = imread(images(2).path);
% Generate grid of points in land polygons
gdu = 100; gdv = 100;
mpu = []; mpv = []; mdu = []; mdv = [];
for i_poly = 1:length(landpolys)
  u = 0:gdu:range(landpolys{i_poly}(:, 1)) + min(landpolys{i_poly}(:, 1));
  v = 0:gdv:range(landpolys{i_poly}(:, 2)) + min(landpolys{i_poly}(:, 2));
  [pu pv] = meshgrid(u, v);
  [in on] = inpolygon(pu, pv, landpolys{i_poly}(:, 1), landpolys{i_poly}(:, 2));
  pu = pu(in & ~on); pv = pv(in & ~on);
  [du, dv] = templatematch(Ia, Ib, pu, pv, 'templatewidth', gdu / 2, 'searchwidth', gdu);
  mpu = [mpu ; pu]; mpv = [mpv ; pv]; mdu = [mdu ; du]; mdv = [mdv ; dv];
end
figure
imshow(Ia / 1.5); hold on;
s = 10; quiver(mpu, mpv, s * mdu, s * mdv, 0, 'y')
% Idealize (undistort) image coordinates
nonans = ~(isnan(mdu) | isnan(mdv));
mpu = mpu(nonans); mpv = mpv(nonans);
mdu = mdu(nonans); mdv = mdv(nonans);
A = cam.normalize([mpu mpv]);
B = cam.normalize([mpu + mdu, mpv + mdv]);
[F, inliersF] = ransacfitfundmatrix(A', B', 0.00000001); % tricky to set threshold (re-express in pixels?)
[H, inliersH] = ransacfithomography(A', B', 0.000001); % tricky to set threshold (re-express in pixels?)
figure
imshow(Ia / 1.5); hold on;
s = 10; quiver(mpu(inliers), mpv(inliers), s * mdu(inliers), s * mdv(inliers), 0, 'y')
% Calibrate with DEM
Xi = cam.normalized2image(A(inliersH, :));
Xw = cam.invproject(Xi, smdem);
Xi2 = cam.normalized2image(B(inliersH, :));
[newcam, rmse, ~] = images(2).cam.optimizecam(Xw, Xi2, '00000111000000000000');
images(2).cam = newcam;

% Generate grid of points in glacier polygon
[Z, ~, bbox] = geotiffread('/volumes/science/data/columbia/dem/2004 Aerometric/20080811_2m.tif');
gdem = DEM(Z, bbox(:, 1), flip(bbox(:, 2)));
viewbox = images(1).cam.viewbox(15 * 1e3);
gcdem = gdem.crop(viewbox(:, 1), viewbox(:, 2), [0 Inf]);
gsmdem = gcdem.resize(0.25);
%
gdu = 100; gdv = 100;
u = 0:gdu:range(glacierpoly(:, 1)) + min(glacierpoly(:, 1));
v = 0:gdv:range(glacierpoly(:, 2)) + min(glacierpoly(:, 2));
[pu pv] = meshgrid(u, v);
[in on] = inpolygon(pu, pv, glacierpoly(:, 1), glacierpoly(:, 2));
pu = pu(in & ~on); pv = pv(in & ~on);
Xw1 = images(1).cam.invproject([pu pv], gsmdem);
Xi2i = images(2).cam.project(Xw1);
[du, dv] = templatematch(Ia, Ib, pu, pv, 'templatewidth', gdu / 2, 'searchwidth', 2 * gdu, 'initialdu', Xi2i(:, 1) - pu, 'initialdv', Xi2i(:, 2) - pv, 'supersample', 1);
Xi2f = [pu + du, pv + dv];
Xw2 = images(2).cam.invproject(Xi2f, smdem);
Xi1f = images(1).cam.project(Xw2);
figure
imshow(Ia / 1.5); hold on;
s = 10; quiver(pu, pv, s * (Xi1f(:, 1) - pu), s * (Xi1f(:, 2) - pv), 0, 'y')
figure
gsmdem.plot(2); hold on;
s = 10; quiver(Xw1(:, 1), Xw1(:, 2), s * (Xw2(:, 1) - Xw1(:, 1)), s * (Xw2(:, 2) - Xw1(:, 2)), 0, 'y')
gsmdem.plot(2); hold on; plot(Xw1(:, 1), Xw2(:, 2), 'r*')


% Calibrate with H
[H, inliersH] = ransacfithomography(A', B', 0.0000001); % tricky to set threshold (re-express in pixels?)
Xi2 = cam.normalized2image(B(inliersH, :));
R = H; % relative rotation
R = ang2rot(rot2ang(H));
pXc2 = inhomog(hnormalise(R * homog(A(inliersH, :)')))';
pXi2 = cam.normalized2image(pXc2);
rmse = sqrt(sum(sum((Xi2(:, 1:2) - pXi2).^2, 2)) / size(pXi2, 1))
% newcam = cam;
% newcam.R = R * newcam.R;
% Xw = bsxfun(@plus, cam.normalized2world(A(inliersH, :)) * 1000, cam.xyz);
% pXi2 = newcam.project(Xw);
% rmse = sqrt(sum(sum((Xi2(:, 1:2) - pXi2).^2, 2)) / size(pXi2, 1))

% Calibrate with LM
R = camera.optimizeR(A(inliersH, :), B(inliersH, :));
pXc2 = inhomog(hnormalise(R * homog(A(inliersH, :)')))';
pXi2 = cam.normalized2image(pXc2);
rmse = sqrt(sum(sum((Xi2(:, 1:2) - pXi2).^2, 2)) / size(pXi2, 1))

% Calibrate with F/E
% [F, inliersF] = ransacfitfundmatrix(A', B', 0.00000001); % tricky to set threshold (re-express in pixels?)
% [U, S, V] = svd(F);
% E = U * diag([1 1 0]) * V';
% [U, S, V] = svd(E);
% W = [0 -1 0; 1 0 0; 0 0 1];
% Ra = U * D * V'
% Rb = U * D' * V'
% ta = U(:, 3)
% tb = -U(:, 3)
% Xi2 = cam.normalized2image(B(inliersF, :));
% pXc2 = inhomog(hnormalise(Rb * (homog(A(inliersF, :)') - repmat(ta * -0.0008, 1, length(inliersF)))))';
% pXi2 = cam.normalized2image(pXc2);
% rmse = sqrt(sum(sum((Xi2(:, 1:2) - pXi2).^2, 2)) / size(pXi2, 1))





figure
imshow(Ib / 1.5); hold on;
uv = newcam.project(Xw);
plot(Xi2(:, 1), Xi2(:, 2), 'yo', uv(:, 1), uv(:, 2), 'ro', [Xi2(:, 1)'; uv(:, 1)'], [Xi2(:, 2)'; uv(:, 2)'], 'y-');


rot2oblang(newcam.R)
rot2oblang(ang2rot(rot2ang(H)) * cam.R)
rot2oblang(cam.R)
% Solve for F to second camera
% B'FA = 0
% B = HA
% R = A-1HA (A is camera matrix)
temp = H * homog(A(inliers, :)'); homog(B(inliers, :)') - hnormalise(temp)
R = H
R2 = R * images(1).cam.K;
temp = homog(B(inliers, :)')' * F * homog(A(inliers, :)');

% rot2oblang(H * cam.R)
% ans =
%      -30.3897907518943         -7.45187045566383        -0.316512082183499

t = -cam.R * cam.xyz'
hnormalise(cam.R * ([1;1;1] - cam.xyz'))
hnormalise(cam.R * [1;1;1] + t)
hnormalise([cam.R t] * homog([1;1;1]))
% equal

[U, S, V] = svd(F);
E = U * [1 0 0; 0 1 0; 0 0 0] * V';
[U, S, V] = svd(E);
D = [0 1 0; -1 0 0; 0 0 1];
Ra = U * D * V'
Rb = U * D' * V'
t = U(:, 3)
R = -vgg_contreps(t) * E;




x1 = homog(A(inliers, :)');
x2 = homog(B(inliers, :)');
[H, inliers] = ransacfithomography(x1, x2, 0.000001); % tricky to set threshold (re-express in pixels?)
x2 - hnormalise(H * x1)
%
images(2).cam.R = H * images(1).cam.R;
figure
imshow(Ib / 1.5), hold on
huv1 = images(1).cam.project(hxyz);
huv2 = images(2).cam.project(hxyz);
plot(huv1(:, 1), huv1(:, 2), 'r-', huv2(:, 1), huv2(:, 2), 'g-')
% % Plot reprojections vs original pixel coordinates
% I = imread(images(1).path);
% screen_size = get(0, 'ScreenSize');
% fig = figure();
% set(fig, 'Position', [0 0 screen_size(3) screen_size(4)]);
% % Original
% subplot(1, 2, 1)
% uv = images(1).cam.project(gcp(:, 3:5));
% rmse = sqrt(sum(sum((gcp(:, 1:2) - uv).^2, 2)) / size(uv, 1));
% imshow(I / 1.5); hold on;
% plot(gcp(:, 1), gcp(:, 2), 'yo', uv(:, 1), uv(:, 2), 'ro', [gcp(:, 1)'; uv(:, 1)'], [gcp(:, 2)'; uv(:, 2)'], 'y-');
% title(sprintf('Reprojection errors (RMSE = %.1fpx)', rmse))
% % Calibrated
% subplot(1, 2, 2)
% uv = newcam.project(gcp(:, 3:5));
% rmse = sqrt(sum(sum((gcp(:, 1:2) - uv).^2, 2)) / size(uv, 1));
% imshow(I / 1.5); hold on;
% plot(gcp(:, 1), gcp(:, 2), 'yo', uv(:, 1), uv(:, 2), 'ro', [gcp(:, 1)'; uv(:, 1)'], [gcp(:, 2)'; uv(:, 2)'], 'y-');
% title(sprintf('Reprojection errors (RMSE = %.1fpx)', rmse))

%% Segment anchor image
[Z, ~, bbox] = geotiffread('/volumes/science/data/columbia/_new/ifsar/merged_projected_clipped.tif');
dem = DEM(Z, bbox(:, 1), flip(bbox(:, 2)));
viewbox = images(1).cam.viewbox(30 * 1e3);
cdem = dem.crop(viewbox(:, 1), viewbox(:, 2));
smdem = cdem.resize(0.2);

visible = voxelviewshed(cdem.X, cdem.Y, cdem.Z, images(1).cam.xyz);
v = visible(:);
vbxyz = [cdem.X(~v) cdem.Y(~v) cdem.Z(~v)];
[vbuv, ~, inframe] = images(1).cam.project(vbxyz);
vbuv = vbuv(inframe, :);
vbxyz = vbxyz(inframe, :);
figure
imshow(I / 1.5); hold on;
plot(vbuv(:, 1), vbuv(:, 2), 'y.')

%% Sky: Horizon detection (works great!)
viewedges = images(1).cam.viewpyramid();
dxy = bsxfun(@minus, viewedges(1:(end - 1), 1:2), viewedges(end, 1:2));
angles = atan2d(dxy(:, 2), dxy(:, 1));
ray_angles = [min(angles):0.1:max(angles)]';
dx = cosd(ray_angles); dy = sind(ray_angles);
rays = [repmat(images(1).cam.xyz, length(ray_angles), 1) dx dy repmat(0, length(ray_angles), 1)];
X = [];
for i = 1:length(ray_angles)
  ray = rays(i, :);
  cells = traverseRayDEM(ray, smdem);
  % Convert to upper-left matrix indices (flip y)
  cells(:, 2) = smdem.ny - (cells(:, 2) - 1);
  xi = cells(:, 1); yi = cells(:, 2);
  % Retrieve true x,y,z based on cell xy
  ind = sub2ind(size(smdem.Z), yi, xi);
  x = smdem.x(xi);
  y = smdem.y(yi);
  z = smdem.Z(ind);
  elevation = atand((z - cam.xyz(3))' ./ sqrt((x - cam.xyz(1)).^2 + (y - cam.xyz(2)).^2));
  [~, i_max] = max(elevation);
  X(i, :) = [x(i_max) y(i_max) z(i_max)];
end
[huv, ~, inframe] = images(1).cam.project(X);
%huv = huv(inframe, :);
figure
imshow(I / 1.5); hold on;
plot(huv(:, 1), huv(:, 2), 'g-');

% Water: Forebay cut by terminus
bshp = shaperead('/volumes/science/data/columbia/outline/forebay.shp');
bxyz = [bshp(1).X' bshp(1).Y' repmat(0, length(bshp(1).X), 1)];
% Terminus
S = shaperead('/volumes/science/others/azinck/terminus.shp');
S = S(find([S.DATE]' == 20080811));
txyz = [S(1).X' S(1).Y' repmat(0, length(S(1).X), 1)];
% Cut forebay with terminus
c = polyclip(bxyz(:, 1:2), txyz(:, 1:2));
wxyz = [c{1}{1} repmat(0, size(c{1}{1}, 1), 1)];
wuv = images(1).cam.project(wxyz);
wuv = clipPolygon(wuv, [0.5 images(1).cam.imgsz(2) 0.5 images(1).cam.imgsz(1)]);
figure
imshow(I / 1.5); hold on
plot(buv(:, 1), buv(:, 2), 'b-');
patch(wuv(:, 1), wuv(:, 2), 1, 'Facecolor', 'b');

% Mountain: Glacier outline
oshp = shaperead('/volumes/science/others/gdurand/outline/20070922-2d.shp');
ind = dem.xy2ind([oshp(1).X', oshp(1).Y']);
z = ind; z(~isnan(ind)) = dem.Z(ind(~isnan(ind))); z(isnan(ind)) = NaN;
oxyz = [oshp(1).X' oshp(1).Y' z];
figure
dem.plot(3); hold on; plot3(oxyz(:, 1), oxyz(:, 2), oxyz(:, 3));
ouv = images(1).cam.project(oxyz);
figure
imshow(I / 1.5); hold on
plot(ouv(:, 1), ouv(:, 2), 'b-');

p = bxyz(:, 1:2);
l = txyz(:, 1:2);
%p = [0 1; 1 1; 1 0; 0 0; 0 1]; % intersect spans polygon closure
p = [0 0; 0 1; 0.25 1; 0.25 1.25; 0.5 1.25; 0.5 1; 1 1; 1 0; 0 0];
%l = [0 1; 0 2; 1 2; 1 1;]; % out only
%l = [-0.5 0.5; 0.5 0.5; 0.5 1.5; 0.75 1.5; 0.75 0.5; 1.5 0.5];
%l = [-0.5 0.5; 0.5 0.5; 0.5 1.5; 0.25 1.5; 0.25 0.75; -1.5 0.75; -1.5 0.25; 0.5 0.25; 1.5 0.25];
%l = [-1 0.25; 0.5 0.25; 1.5 0.25];
%l = [-0.5 0.5; 0.5 0.5; 0.5 1.5; 0.25 1.5; 0.25 0.75; -1 0.75; -1 0.25; 0 0.25];
l = [-0.5 1; 0.5 1; 0.5 0];
c = polyclip(p, l);
plot(p(:, 1), p(:, 2), 'k-', l(:, 1), l(:, 2), 'k-'); hold on
% Xi = intersectPolylines(p, l);
% plot(Xi(:, 1), Xi(:, 2), 'r*')
for i_poly = 1:length(c{1})
  patch(c{1}{i_poly}(:, 1), c{1}{i_poly}(:, 2), 1, 'facecolor', 'g')
  % patch(c{1}{i_poly}(:, 1), c{1}{i_poly}(:, 2), i_poly)
end
for i_poly = 1:length(c{2})
  patch(c{2}{i_poly}(:, 1), c{2}{i_poly}(:, 2), 1, 'facecolor', 'r')
  % patch(c{2}{i_poly}(:, 1), c{2}{i_poly}(:, 2), i_poly + length(c{1}))
end
plot(p(:, 1), p(:, 2), 'k-', l(:, 1), l(:, 2), 'k:');
%% works!
%% TODO:
% Prune polygons (choose largest area?)
% Adapt for polyline-polyline intersection


Xi = intersectPolylines(bxyz(:, 1:2), txyz(:, 1:2));
plot(bxyz(:, 1), bxyz(:, 2), 'b-'); hold on
plot(txyz(:, 1), txyz(:, 2), 'g-');
plot(Xi(:, 1), Xi(:, 2), 'r*')


tpxyz = [txyz; flip(txyz(1:(end - 1), :)) - 100; txyz(1, :)];
[xd, yd] = polybool('xor', bxyz(:, 1), bxyz(:, 2), tpxyz(:, 1), tpxyz(:, 2));
wxy = [xd yd]; nans = find(isnan(xd)); wxy = wxy(1:nans(1)-1, :);
fill(wxy(:, 1), wxy(:, 2), 'red')
patch(xd, yd, 1, 'Facecolor', 'r')
plot(bxyz(:, 1), bxyz(:, 2), 'g-', tpxyz(:, 1), tpxyz(:, 2), 'r-')
wuv = cam.project([wxy repmat(0, size(wxy, 1), 1)]);






% Glacier outline
gshp = shaperead('/volumes/science/data/columbia/outline/2002_2010/20080811.shp');
z = interp2(smdem.X, smdem.Y, smdem.Z, gshp(1).X', gshp(1).Y');
gxyz = double([gshp(1).X' gshp(1).Y' z]); % HACK: double() needed for polybool??
guv = images(1).cam.project(gxyz);
figure
imshow(I / 1.5); hold on;
plot(guv(:, 1), guv(:, 2), 'w-');
plot(gxyz(:, 1), gxyz(:, 2))

% Viewpoint visibility
visible = voxelviewshed(smdem.X, smdem.Y, smdem.Z, images(1).cam.xyz);
% Camera visibility
v = visible(:);
vxyz = [smdem.X(v) smdem.Y(v) smdem.Z(v)];
[vuv, ~, inframe] = images(1).cam.project(vxyz);
vuv = vuv(inframe, :);
vxyz = vxyz(inframe, :);
figure
imshow(I / 1.5); hold on;
plot(vuv(:, 1), vuv(:, 2), 'g.');



% ifSAR DEM
[Z, ~, bbox] = geotiffread('/volumes/science/data/columbia/_sort/ifsar/merged_projected_clipped.tif');
sardem = DEM(Z, bbox(:, 1), flip(bbox(:, 2)));
%bdem = bdem.crop([], [], [0 Inf]);
visible = voxelviewshed(bdem.X, bdem.Y, bdem.Z, cam.xyz);
bdem.plot(3); hold on
h = bdem.plot(3, visible * 255);
set(h, 'AlphaData', 0.2);
plot3(cam.xyz(1), cam.xyz(2), cam.xyz(3), 'r*', 'markersize', 20)

% Forebay outline
bshp = shaperead('/volumes/science/data/columbia/outline/forebay.shp');
bxyz = [bshp(1).X' bshp(1).Y' repmat(0, length(bshp(1).X), 1)];
buv = cam.project(bxyz);

% Glacier outline
gshp = shaperead('/volumes/science/data/columbia/outline/2002_2010/20050811.shp');
z = interp2(smdem.X, smdem.Y, smdem.Z, gshp(1).X', gshp(1).Y');
gxyz = double([gshp(1).X' gshp(1).Y' z]); % HACK: double() needed for polybool??
guv = cam.project(gxyz);



% Subtractions (bay - glacier)
[xd, yd] = polybool('subtraction', bxyz(:, 1), bxyz(:, 2), gxyz(:, 1), gxyz(:, 2));
wxy = [xd yd]; nans = find(isnan(xd)); wxy = wxy((nans(3) + 1):end, :);
plot(bxyz(:, 1), bxyz(:, 2), 'k-', gxyz(:, 1), gxyz(:, 2), 'k-', wxy(:, 1), wxy(:, 2), 'r-')
wuv = cam.project([wxy repmat(0, size(wxy, 1), 1)]);

figure
imshow(I / 1.5); hold on;
plot(wuv(:, 1), wuv(:, 2), 'b-', vb(:, 1), vb(:, 2), 'k-', guv(:, 1), guv(:, 2), 'w-')
plot(wuv(:, 1), wuv(:, 2), 'b-');



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
