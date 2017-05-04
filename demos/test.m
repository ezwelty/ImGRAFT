addpath(genpath('.'))

%% --- Constants ---

% Paths
DATA_DIR = fullfile('demos', 'data');
IMG_PATH = fullfile(DATA_DIR, 'images', 'test', '*.JPG');
DEM_PATH = fullfile(DATA_DIR, 'dems', '20090827_2m.tif');
OUT_DIR = fullfile('demos', 'test');

% Parameters
IMG_SCALE = 0.25;
DEM_SCALE = 0.25;
DEM_DISTANCE = 20e3; % m
DEM_SMOOTH_WIDTH = 300; % m
COAST_HAE = 17; % m
CAM_ID = lower('AK10');
SVG_GCP = 'gcp';
SVG_COAST = 'coast';
SVG_HORIZON = 'horizon';
SVG_LAND = 'land';
SVG_GLACIER = 'glacier';
SVG_LINE_DENSITY = 0.10; % points per pixel

% Datasets
CAM = readtable(fullfile(DATA_DIR, 'cameras.csv'));
CAM.Properties.RowNames = lower(CAM.id);
GCP = readtable(fullfile(DATA_DIR, 'points.csv'));
GCP.Properties.RowNames = lower(GCP.name);
COAST = shaperead(fullfile(DATA_DIR, 'lines', 'coastline.shp'));
COAST = arrayfun(@(s) [s.X' s.Y' repmat(COAST_HAE, length(s.X), 1)], COAST, 'uniform', false);

%% --- Initialize ----

%% Load camera (initial)
cam_xyz = CAM{CAM_ID, {'x_wgs84', 'y_wgs84', 'z_hae'}};
cam_viewdir = CAM{CAM_ID, {'yaw', 'pitch', 'roll'}};
cam = Camera('xyz', cam_xyz, 'viewdir', cam_viewdir);

%% Load images (and SVG markup)
% gcp . <name> : [x y]
% land . polygon_<i> : [x1 y1; ... ; xn yn]
% glacier . polygon_<i> : [x1 y1; ... ; xn yn]
% coast . polyline_<i> : [x1 y1; ... ; xn yn]
% horizon . polyline_<i> : [x1 y1; ... ; xn yn]
% FIXME: xmlread (for svg) very slow
images = Image(IMG_PATH, cam);

%% Format GCP (ground control points)
% gcp.uv : [u1 v1; ... ]
% gcp.xyz : [x1 y1 z1; ... ]
has_gcp = find(arrayfun(@(img) isfield(img.svg, SVG_GCP), images));
for i_img = has_gcp
  gcp_names = fieldnames(images(i_img).svg.gcp);
  images(i_img).gcp.uv = cell2mat(struct2cell(images(i_img).svg.gcp));
  images(i_img).gcp.xyz = GCP{lower(gcp_names), {'x_wgs84', 'y_wgs84', 'z_hae'}};
end

%% Load DEM
% Crop to camera view
% FIXME: Won't work if camera missing initial view direction and position
box = images(1).cam.viewbox(DEM_DISTANCE);
dem = DEM(DEM_PATH).crop(box(:, 1), box(:, 2), [0 Inf]).resize(DEM_SCALE).build();
smdem = dem.smooth(DEM_SMOOTH_WIDTH).build();

%% Load Horizon
% FIXME: Slow
% TODO: Cache results by CAM_ID
% TODO: Traverse grid in concentric circles or triangular wedges
% FIXME: Won't work if camera missing initial position
if isempty(cam.viewdir)
  HORIZON = dem.horizon(cam.xyz, 0:0.1:359);
else
  HORIZON = images(1).cam.horizon(dem, 0.1);
end

%% Format "GCL" (ground control lines)
% gcl.uv : [u1 v1; ... ] (as point matrix)
% gcl.xyz : {[x1 y1 z1; ... ], ... } (as polyline array)
has_lines = find(arrayfun(@(img) any(isfield(img.svg, {SVG_COAST, SVG_HORIZON})), images));
for i = has_lines
  [coast, horizon] = deal({});
  if isfield(images(i).svg, SVG_COAST)
    coast = struct2cell(images(i).svg.(SVG_COAST));
  end
  if isfield(images(i).svg, SVG_HORIZON)
    horizon = struct2cell(images(i).svg.(SVG_HORIZON));
  end
  lines = [coast; horizon];
  % Resample traced lines
  for j = 1:length(lines)
    line_length = polylineLength(lines{j});
    n_points = round(line_length * SVG_LINE_DENSITY);
    lines{j} = resamplePolyline(lines{j}, n_points);
  end
  % Reduce to single matrix of points
  images(i).gcl.uv = cell2mat(lines);
  % Attach control lines
  % FIXME: Redundant and not memory efficient
  images(i).gcl.xyz = [ifelse(isempty(coast), {}, COAST); ifelse(isempty(horizon), {}, HORIZON)];
end

%% Format fixed (land) and free (glacier) polygons
% freepolys = {[x1 y1 z1; ...], ...}
% fixedpolys = {[x1 y1 z1; ...], ...}
for i = 1:length(images)
  if isfield(images(i).svg, SVG_LAND)
    for j = fieldnames(images(i).svg.(SVG_LAND))'
      images(i).fixedpolys{end + 1} = images(i).svg.(SVG_LAND).(j{1});
    end
  end
  if isfield(images(i).svg, SVG_GLACIER)
    for j = fieldnames(images(i).svg.(SVG_GLACIER))'
      images(i).freepolys{end + 1} = images(i).svg.(SVG_GLACIER).(j{1});
    end
  end
end

%% --- Calibrate anchor images ----

is_anchor = arrayfun(@(img) any(~isempty(img.gcp.xyz) & ~isempty(img.gcl.xyz)), images);
anchor_ind = find(is_anchor);
flexparams = {'viewdir'};
fixparams = {'f', 'k', [1 2]};
[anchors, fit] = Camera.optimize_images(images(is_anchor), flexparams, fixparams, 10);

% Review results
for i = 1:length(anchors)
  figure();
  anchors(i).plot(true);
  title(['\fontsize{14} ', num2str(anchor_ind(i)), ': RMSE ', num2str(fit.rmse(i), '%.2f'), 'px']);
end

% Save results
images(is_anchor) = anchors;

%% --- Assign images to anchors ---
% Match each image to the nearest (temporal) anchor.

anchor_dates = [images(is_anchor).date_num];
for i = find(~is_anchor)
  [~, i_min] = min(abs(anchor_dates - images(i).date_num));
  images(i).anchor = anchor_ind(i_min);
end

%% --- Orient images from anchors ---
% Match features between image and anchor, filter with RANSAC, and optimize orientation.
% TODO: Support variable image sizes

fits = cell(size(images));
for i = find(~is_anchor)
  I = images(i).read();
  i0 = images(i).anchor;
  I0 = images(i0).read();
  % Generate grid of points in land polygons
  % TODO: Move outside loop
  % TODO: Switch to keypoints matching
  [gdu, gdv] = deal(50);
  matches = [];
  for j = 1:length(images(i0).fixedpolys)
    pts = polygon2grid(images(i0).fixedpolys{j}, gdu, gdv);
    [du, dv] = templatematch(I0, I, pts(:, 1), pts(:, 2), 'templatewidth', gdu, 'searchwidth', 5 * gdu);
    matches = vertcat(matches, horzcat(pts, du, dv));
  end
  matches(any(isnan(matches), 2), :) = [];
  % Filter matches with RANSAC
  % FIXME: Assumes images are the same size and camera
  xy0 = images(i0).cam.image2camera(matches(:, 1:2));
  xy = images(i0).cam.image2camera(matches(:, 1:2) + matches(:, 3:4));
  % TODO: Express threshold in pixels
  [F, in] = ransacfitfundmatrix(xy0', xy', 0.0000005);
  mean_motion = mean(sqrt(sum(matches(in, 3:4).^2, 2)));
  % Plot filtered matches
  figure
  imshow(I0 / 1.5), hold on
  s = 1;
  quiver(matches(:, 1), matches(:, 2), matches(:, 3), matches(:, 4), 0, 'r');
  quiver(matches(in, 1), matches(in, 2), matches(in, 3), matches(in, 4), 0, 'y');
  title([num2str(i0), ' -> ', num2str(i)]);
  % Orient camera
  % FIXME: Assumes images are the same size and camera
  % TODO: Retrieve rotation from F?
  % TODO: Write dedicated orientation function using normalized camera coordinates
  if mean_motion > 1
    [newcams, fits{i}] = Camera.optimize_bundle(images(i0).cam, matches(in, 1:2) + matches(in, 3:4), images(i0).cam.xyz + images(i0).cam.camera2world(xy0(in, :)), 'viewdir');
    images(i).cam = newcams{1};
  else
    images(i).cam = images(i0).cam;
  end
  % Transform free polygons
  for j = 1:length(images(i0).freepolys)
    images(i).freepolys{j} = images(i).cam.project(images(i0).cam.invproject(images(i0).freepolys{j}), true);
  end
end

%% --- Save aligned images ---
% TODO: Undistort images at the same time?
% TODO: Support variable image sizes
% TODO: Wrap into function

mkdir(OUT_DIR, 'aligned');
for i0 = find(is_anchor)
  % Copy file
  [~, filename, ext] = fileparts(images(i0).file);
  outfile = fullfile(OUT_DIR, 'aligned', [num2str(i0), '-', num2str(i0), '-', filename, ext]);
  copyfile(images(i0).file, outfile);
  % Prepare reference grid
  [u, v] = meshgrid(0.5:(images(i0).cam.imgsz(1) - 0.5), 0.5:(images(i0).cam.imgsz(2) - 0.5));
  uv = [u(:), v(:)];
  dxyz = images(i0).cam.invproject(uv);
  for i = find([images.anchor])
    % Project reference grid to new image
    puv = images(i).cam.project(dxyz, true);
    % Interpolate image at projected grid points
    I = double(images(i).read());
    I0 = uint8(nan(size(I)));
    for channel = 1:size(I, 3)
      temp = interp2(u, v, I(:, :, channel), puv(:, 1), puv(:, 2), '*linear');
      I0(:, :, channel) = reshape(temp, flip(images(i0).cam.imgsz));
    end
    % Save to file
    [~, filename, ext] = fileparts(images(i).file);
    outfile = fullfile(OUT_DIR, 'aligned', [num2str(i), '-', num2str(i0), '-', filename, ext]);
    imwrite(I0, outfile, 'Quality', 100);
  end
end

%% --- Match moving features ----
% Uses the Farneback optical flow method
% see: https://www.mathworks.com/help/vision/ref/opticalflowfarneback-class.html

% Match glacier features between consecutive image pairs
obj = opticalFlowFarneback('NumPyramidLevels', 5, 'NeighborhoodSize', 3, 'FilterSize', 25);
flow = cell(length(images));
for i = 1:length(images)
  I = rgb2gray(images(i).read());
  flow{i} = obj.estimateFlow(I);
end

% Generate glacier points
smdem = dem.resize(0.125).smooth(300).build();
[gdx, gdy] = deal(10);
gxyz = [];
i0 = 1;
for poly = images(i0).freepolys
  xyz = images(i0).cam.invproject(poly{:});
  for i_pt = 1:size(xyz, 1)
    xyz(i_pt, :) = smdem.sample_ray_tri(images(i0).cam.xyz, xyz(i_pt, :), true);
  end
  gxy = polygon2grid(xyz, gdx, gdy);
  z = smdem.sample_points_tri(gxy);
  gxyz0 = [gxy, z];
  uv = images(i0).cam.project(gxyz0);
  gxyz1 = images(i0).cam.invproject(uv);
  for i_pt = 1:size(gxyz1, 1)
    gxyz1(i_pt, :) = smdem.sample_ray_tri(images(i0).cam.xyz, gxyz1(i_pt, :), true);
  end
  visible = sqrt(sum((gxyz1 - gxyz0).^2, 2)) < 1;
  gxyz = [gxyz; [gxyz0(visible, :)]];
end

figure()
smdem.plot(2)
hold on
plot(gxyz(:, 1), gxyz(:, 2), 'r*');
guv = images(i0).cam.project(gxyz);
figure()
images(i0).plot();
hold on
plot(guv(:, 1), guv(:, 2), 'r*')

% Sample motion at glacier points
for i0 = 1:(length(images) - 1)
  i = i0 + 1;
  % Starting positions (i0)
  uv0 = images(i0).cam.project(gxyz);
  [u, v] = meshgrid(0.5:(images(i0).cam.imgsz(1) - 0.5), 0.5:(images(i0).cam.imgsz(2) - 0.5));
  du = interp2(u, v, flow{i}.Vx, uv0(:, 1), uv0(:, 2), '*cubic');
  dv = interp2(u, v, flow{i}.Vy, uv0(:, 1), uv0(:, 2), '*cubic');
  % Ending positions (i)
  uv = uv0 + [du, dv];
  xyz = images(i).cam.invproject(uv);
  for i_pt = 1:size(xyz, 1)
    xyz(i_pt, :) = smdem.sample_ray_tri(images(i).cam.xyz, xyz(i_pt, :), true, gxyz(i_pt, 1:2), 100);
  end
  figure();
  showimg(dem.x, dem.y, hillshade(dem.Z, dem.x, dem.y));
  hold on
  ddays = images(i).date_num - images(i0).date_num;
  v = sqrt(sum((xyz - gxyz).^2, 2)) / ddays;
  [X, Y, V] = pts2grid(gxyz(:, 1), gxyz(:, 2), v, gdx * 5, gdy * 5);
  alphawarp(X, Y, V, 1);
  caxis([0 10]);
  colormap jet;
  colorbar
  [~, ~, DX] = pts2grid(gxyz(:, 1), gxyz(:, 2), xyz(:, 1) - gxyz(:, 1), gdx * 5, gdy * 5);
  [~, ~, DY] = pts2grid(gxyz(:, 1), gxyz(:, 2), xyz(:, 2) - gxyz(:, 2), gdx * 5, gdy * 5);
  s = 1;
  quiver(X, Y, s * DX / ddays, s * DY / ddays, 0, 'k');
  figure();
  ind = 1:size(uv, 1);
  ind(isnan(v(ind))) = [];
  ind(v(ind) > 10) = [];
  images(i0).plot();
  hold on
  quiver(uv0(ind, 1), uv0(ind, 2), du(ind), dv(ind), 0, 'w');
  colors = jet(101);
  scatter(uv0(ind, 1), uv0(ind, 2), 20, colors(ceil(v(ind) * 10), :), 'markerFaceColor', 'flat');
  caxis([0 10]);
  colormap jet;
  colorbar
end

% Alternate: plot all results

for i0 = 1:(length(images) - 1)
  % Starting positions (image i0)
  mask = zeros(flip(images(i0).cam.imgsz));
  for j = 1:length(images(i0).freepolys)
    mask = mask | poly2mask(images(i0).freepolys{j}(:, 1), images(i0).freepolys{j}(:, 2), images(i0).cam.imgsz(2), images(i0).cam.imgsz(1));
  end
  [yi, xi] = find(mask);
  uv0 = [xi, yi];
  % Ending positions (image i)
  i = i0 + 1;
  du = flow{i}.Vx;
  dv = flow{i}.Vy;
  duv = [du(mask), dv(mask)];
  uv = uv0 + duv;

  % figure();
  % imshow(images(i).read());
  % hold on
  % quiver(uv0(ind, 1), uv0(ind, 2), duv(ind, 1), duv(ind, 2), 0, 'y');

  % Reduce and convert to world coordinates
  ind = round(rand(1e6, 1) * size(uv0, 1));
  plane = [0 0 -1 100];
  dxyz0 = images(i0).cam.invproject(uv0(ind, :));
  xyz0 = intersectRayPlane(images(i0).cam.xyz, dxyz0, plane);
  dxyz = images(i).cam.invproject(uv(ind, :));
  xyz = intersectRayPlane(images(i).cam.xyz, dxyz, plane);

  % smdem.plot(2);
  figure();
  showimg(dem.x, dem.y, hillshade(dem.Z, dem.x, dem.y));
  hold on
  ddays = images(i).date_num - images(i0).date_num;
  v = sqrt(sum((xyz - xyz0).^2, 2)) / ddays;
  [X, Y, V] = pts2grid(xyz0(:, 1), xyz0(:, 2), v, 100, 100);
  alphawarp(X, Y, V, 1);
  caxis([0 10]);
  colormap jet;
  colorbar
  [~, ~, DX] = pts2grid(xyz0(:, 1), xyz0(:, 2), xyz(:, 1) - xyz0(:, 1), 100, 100);
  [~, ~, DY] = pts2grid(xyz0(:, 1), xyz0(:, 2), xyz(:, 2) - xyz0(:, 2), 100, 100);
  s = 2;
  quiver(X, Y, s * DX / ddays, s * DY / ddays, 0, 'k');
end

%% --- Match moving features (deprecated) ----

% Set image used as coordinate reference
i_ref = 1;

% Generate grid of points in glacier polygons
guv = [];
[du, dv] = deal(50);
for j = 1:length(images(i_ref).freepolys)
  guv = [guv; polygon2grid(images(i_ref).freepolys{j}, du, dv)];
end
gdxyz = images(i_ref).cam.invproject(guv);

% Match glacier features between consecutive image pairs
matches = cell(length(images) - 1, 1);
for i0 = 1:(length(images) - 1)
  % Prepare first image
  I0 = images(i0).read();
  p0 = images(i0).cam.project(gdxyz, true);
  % Prepare second image
  i = i0 + 1;
  I = images(i).read();
  p = images(i).cam.project(gdxyz, true);
  % Match points
  % TODO: Use reference velocity field to initialize and constrict search
  % TODO: Try matching methods that account for rotation
  [du, dv] = templatematch(I0, I, p0(:, 1), p0(:, 2), 'templatewidth', 50, 'searchwidth', 250,'initialdu', p(:, 1) - p0(:, 1), 'initialdv', p(:, 2) - p0(:, 2), 'supersample', 1, 'method', 'OC');
  dxyz = images(i).cam.invproject(p0 + [du, dv]);
  matches{i0} = images(i_ref).cam.project(dxyz, true);
end

i = 1;
figure();
imshow(images(i_ref).read());
hold on
duv = matches{i} - guv;
quiver(guv(:, 1), guv(:, 2), duv(:, 1), duv(:, 2), 0, 'y');
