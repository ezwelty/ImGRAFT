%% Velocity from Time-lapse
% Images from a time-lapse camera are calibrated, aligned, and then used to
% measure glacier motion.

%% Assign constants

%%
% Paths
IMGRAFT_PATH = '.';
DATA_DIR = fullfile('demos', 'data');
IMG_PATH = fullfile(DATA_DIR, 'images', 'AK10b_20150324-radar', '*.JPG');
GLACIER_DEM_PATH = fullfile(DATA_DIR, 'dems', 'SETSM_WV01_20150423_102001003C038400_102001003C64D800_seg1_2m_v1.0_dem_proj.tif');
HORIZON_DEM_PATH = fullfile(DATA_DIR, 'dems', '20090827_2m.tif');
OUT_DIR = fullfile('demos', 'AK10b_20150324-radar');
addpath(genpath(IMGRAFT_PATH));

%%
% Parameters
IMG_SCALE = 1; % FIXME: currently unused
HORIZON_DEM_DISTANCE = 20e3; % m
HORIZON_DEM_SCALE = 0.25;
GLACIER_DEM_DISTANCE = 5e3; % m
GLACIER_DEM_SCALE = 0.25;
GLACIER_DEM_SMOOTH_WIDTH = 300; % m
COAST_HAE = 17; % m
CAM_ID = lower('AK10b');
SVG_GCP = 'gcp';
SVG_COAST = 'coast';
SVG_HORIZON = 'horizon';
SVG_LAND = 'land';
SVG_GLACIER = 'glacier';
SVG_LINE_DENSITY = 0.05; % points per pixel
LDMAX = 10; % pixels (max distance to assign when calibrating with lines)

%%
% Datasets
CAM = readtable(fullfile(DATA_DIR, 'cameras.csv'));
CAM.Properties.RowNames = lower(CAM.id);
GCP = readtable(fullfile(DATA_DIR, 'points.csv'));
GCP.Properties.RowNames = lower(GCP.name);
COAST = shaperead(fullfile(DATA_DIR, 'lines', 'coastline.shp'));
COAST = arrayfun(@(s) [s.X' s.Y' repmat(COAST_HAE, length(s.X), 1)], COAST, 'uniform', false);

%% Prepare variables

%%
% Initialize the camera with a surveyed position and estimated direction.
xyz = CAM{CAM_ID, {'x_wgs84', 'y_wgs84', 'z_hae'}};
viewdir = CAM{CAM_ID, {'yaw', 'pitch', 'roll'}};
cam = Camera('xyz', xyz, 'viewdir', viewdir);

%%
% Load images and SVG markup from file.
%
%  gcp . <name> : [x y]
%  land . polygon_<i> : [x1 y1; ... ; xn yn]
%  glacier . polygon_<i> : [x1 y1; ... ; xn yn]
%  coast . polyline_<i> : [x1 y1; ... ; xn yn]
%  horizon . polyline_<i> : [x1 y1; ... ; xn yn]
%
%  FIXME: xmlread (for svg) very slow.
%
images = Image(IMG_PATH, cam);

%%
% Format ground control points (GCP).
%
%  gcp.uv : [u1 v1; ... ]
%  gcp.xyz : [x1 y1 z1; ... ]
%
has_gcp = find(arrayfun(@(img) isfield(img.svg, SVG_GCP), images));
for i_img = has_gcp
  gcp_names = fieldnames(images(i_img).svg.gcp);
  images(i_img).gcp.uv = cell2mat(struct2cell(images(i_img).svg.gcp));
  images(i_img).gcp.xyz = GCP{lower(gcp_names), {'x_wgs84', 'y_wgs84', 'z_hae'}};
end

%%
% Load digital elevation model (DEM) from file, crop to camera, and smooth
% to fill crevasses.
%
%  FIXME: Won't work if camera missing initial view direction and position
%
box = images(1).cam.viewbox(HORIZON_DEM_DISTANCE);
HORIZON_DEM = DEM(HORIZON_DEM_PATH).crop(box(:, 1), box(:, 2), [0 Inf]).resize(HORIZON_DEM_SCALE).build();
box = images(1).cam.viewbox(GLACIER_DEM_DISTANCE);
GLACIER_DEM = DEM(GLACIER_DEM_PATH).crop(box(:, 1), box(:, 2), [0 Inf]).resize(GLACIER_DEM_SCALE).smooth(GLACIER_DEM_SMOOTH_WIDTH).build();

%%
% Compute the world coordinates of the horizon as seen from this camera.
%
%  FIXME: Slow
%  TODO: Cache results by CAM_ID
%  TODO: Traverse grid in concentric circles or triangular wedges
%
if isempty(cam.viewdir)
  HORIZON = HORIZON_DEM.horizon(cam.xyz, 0:0.1:359);
else
  HORIZON = images(1).cam.horizon(HORIZON_DEM, 0.1);
end

%%
% Compile control lines (same for all images).
LXYZ = [ifelse(isempty(coast), {}, COAST); ifelse(isempty(horizon), {}, HORIZON)];

%%
% Format "ground control lines" (GCL).
%
%  gcl.uv : [u1 v1; ... ] (as point matrix)
%  gcl.xyz : {[x1 y1 z1; ... ], ... } (as polyline array)
%
has_lines = find(arrayfun(@(img) any(isfield(img.svg, {SVG_COAST, SVG_HORIZON})), images));
for i = has_lines

  %%
  % Collect line traces from SVG data.
  [coast, horizon] = deal({});
  if isfield(images(i).svg, SVG_COAST)
    coast = struct2cell(images(i).svg.(SVG_COAST));
  end
  if isfield(images(i).svg, SVG_HORIZON)
    horizon = struct2cell(images(i).svg.(SVG_HORIZON));
  end
  lines = [coast; horizon];

  %%
  % Sample points along lines.
  for j = 1:length(lines)
    line_length = polylineLength(lines{j});
    n_points = round(line_length * SVG_LINE_DENSITY);
    lines{j} = resamplePolyline(lines{j}, n_points);
  end
  images(i).gcl.uv = cell2mat(lines);

  %%
  % Attach world coordinates of lines (now only used by image.plot).
  images(i).gcl.xyz = LXYZ;
end

%%
% Format "fixed" (land) and "free" (glacier) polygons.
%
%  freepolys = {[x1 y1 z1; ...], ...}
%  fixedpolys = {[x1 y1 z1; ...], ...}
%
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

%% Calibrate anchor images
% Use available ground control points and lines to optimize all anchor
% images at once.
is_anchor = arrayfun(@(img) any(~isempty(img.gcp.xyz) & ~isempty(img.gcl.xyz)), images);
anchor_ind = find(is_anchor);
flexparams = {'viewdir'};
fixparams = {'f', 'k', 1};
[anchors, fit] = Camera.optimize_images(images(is_anchor), flexparams, fixparams, LDMAX, LXYZ);

%%
% Plot the results.
for i = 1:length(anchors)
  figure();
  anchors(i).plot(true);
  title(['\fontsize{14} ', num2str(anchor_ind(i)), ': RMSE ', num2str(fit.rmse(i), '%.2f'), 'px']);
end

%%
% Save the results.
images(is_anchor) = anchors;

%% Orient images relative to anchors
% Assign each image to the nearest (temporal) anchor.
anchor_dates = [images(is_anchor).date_num];
for i = find(~is_anchor)
  [~, i_min] = min(abs(anchor_dates - images(i).date_num));
  images(i).anchor = anchor_ind(i_min);
end

%%
% Match features between image and anchor, filter with RANSAC, and optimize orientation.
%
%  TODO: Support variable image size
%
for i = find(~is_anchor)
  I = images(i).read();
  i0 = images(i).anchor;
  I0 = images(i0).read();

  %%
  % Generate grid of points in land polygons and match between images.
  %
  %  TODO: Move point generation to anchor-level loop
  %
  [gdu, gdv] = deal(10);
  pts = [];
  for j = 1:length(images(i0).fixedpolys)
    pts = [pts; polygon2grid(images(i0).fixedpolys{j}, gdu, gdv)];
  end
  [du, dv, correlation, signal, pu, pv] = templatematch(I0, I, pts(:, 1), pts(:, 2), 'templatewidth', gdu, 'searchwidth', 2 * gdu, 'method', 'NCC', 'super', 2);
  is_strong = correlation > quantile(correlation, 0.8) & signal > quantile(signal, 0.8);
  matches = horzcat(pu, pv, du, dv);
  matches = matches(is_strong, :);

  %%
  % Filter matches with RANSAC.
  %
  %  FIXME: Assumes images are the same size and camera
  %  TODO: Express threshold in pixels
  %
  xy0 = images(i0).cam.image2camera(matches(:, 1:2));
  xy = images(i0).cam.image2camera(matches(:, 1:2) + matches(:, 3:4));
  [F, in] = ransacfitfundmatrix(xy0', xy', 1e-6);
  mean_motion = mean(sqrt(sum(matches(in, 3:4).^2, 2)));

  %%
  % Plot filtered matches.
  figure
  imshow(I0 / 1.5), hold on
  s = 5;
  quiver(matches(:, 1), matches(:, 2), s * matches(:, 3), s * matches(:, 4), 0, 'r');
  quiver(matches(in, 1), matches(in, 2), s * matches(in, 3), s * matches(in, 4), 0, 'y');
  title([num2str(i0), ' -> ', num2str(i)]);

  %%
  % Orient camera
  %
  %  FIXME: Assumes images are the same size and camera
  %  TODO: Retrieve rotation from F
  %  TODO: Write dedicated orientation function using normalized camera coordinates
  %
  uv = matches(in, 1:2) + matches(in, 3:4);
  xyz = images(i0).cam.xyz + images(i0).cam.camera2world(xy0(in, :));
  images(i).cam = Camera.optimize_bundle(images(i0).cam, uv, xyz, 'viewdir'){1};

  %%
  % Transform free polygons.
  for j = 1:length(images(i0).freepolys)
    images(i).freepolys{j} = images(i).cam.project(images(i0).cam.invproject(images(i0).freepolys{j}), true);
  end
end

%%
% Save aligned images to file.
%
%  TODO: Also undistort images
%  TODO: Support variable image sizes
%
mkdir(OUT_DIR, 'aligned');
for i0 = find(is_anchor)

  %%
  % Copy anchor image.
  [~, filename, ext] = fileparts(images(i0).file);
  outfile = fullfile(OUT_DIR, 'aligned', [num2str(i0), '-', num2str(i0), '-', filename, ext]);
  copyfile(images(i0).file, outfile);

  %%
  % Precompute reference grid.
  [u, v] = meshgrid(0.5:(images(i0).cam.imgsz(1) - 0.5), 0.5:(images(i0).cam.imgsz(2) - 0.5));
  uv = [u(:), v(:)];
  dxyz = images(i0).cam.invproject(uv);

  %%
  % For each image, project grid, interpolate at points, and save to file.
  for i = find([images.anchor] == i0)
    % puv = images(i).cam.project(dxyz, true);
    % I = double(images(i).read());
    % I0 = uint8(nan(size(I)));
    % for channel = 1:size(I, 3)
    %   temp = interp2(u, v, I(:, :, channel), puv(:, 1), puv(:, 2), '*linear');
    %   I0(:, :, channel) = reshape(temp, flip(images(i0).cam.imgsz));
    % end
    I0 = images(i).project(images(i0).cam, dxyz);
    [~, filename, ext] = fileparts(images(i).file);
    outfile = fullfile(OUT_DIR, 'aligned', [num2str(i), '-', num2str(i0), '-', filename, ext]);
    imwrite(I0, outfile, 'Quality', 100);
  end
end

%% Track moving features
% Match features between consecutive image pairs using the Farneback
% optical flow method:
% <https://www.mathworks.com/help/vision/ref/opticalflowfarneback-class.html>.
%
%  TODO: Set pyramid levels based on expected motion.
%
obj = opticalFlowFarneback('NumPyramidLevels', 3, 'NeighborhoodSize', 3, 'FilterSize', 25);
flow = cell(length(images));
for i = 1:length(images)
  I = rgb2gray(images(i).read());
  flow{i} = obj.estimateFlow(I);
end
% for i = 2:length(images)
%   images(i).plot();
%   hold on;
%   plot(flow{i}, 'Decimation', [20, 20], 'Scale', 20);
%   pause;
% end

%%
% Generate a regular grid of glacier points.
[gdx, gdy] = deal(10);
gxyz = [];
i0 = 1;
for poly = images(i0).freepolys
  xyz = images(i0).cam.invproject(poly{:});
  for i_pt = 1:size(xyz, 1)
    xyz(i_pt, :) = GLACIER_DEM.sample_ray_tri(images(i0).cam.xyz, xyz(i_pt, :), true);
  end
  gxy = polygon2grid(xyz, gdx, gdy);
  z = GLACIER_DEM.sample_points_tri(gxy);
  gxyz0 = [gxy, z];
  uv = images(i0).cam.project(gxyz0);
  gxyz1 = images(i0).cam.invproject(uv);
  for i_pt = 1:size(gxyz1, 1)
    gxyz1(i_pt, :) = GLACIER_DEM.sample_ray_tri(images(i0).cam.xyz, gxyz1(i_pt, :), true);
  end
  visible = sqrt(sum((gxyz1 - gxyz0).^2, 2)) < 1;
  gxyz = [gxyz; [gxyz0(visible, :)]];
end

%%
% Plot glacier points on map.
figure()
HORIZON_DEM.plot(2);
hold on
plot(gxyz(:, 1), gxyz(:, 2), 'r.');

%%
% Plot glacier points on image.
figure()
images(i0).plot();
hold on
guv = images(i0).cam.project(gxyz);
plot(guv(:, 1), guv(:, 2), 'y.');

%%
% Sample motion at glacier points.
motion = struct();
for i0 = 1:(length(images) - 1)
  i = i0 + 1;

  %%
  % Calculate starting positions (i0).
  uv0 = images(i0).cam.project(gxyz);
  [u, v] = meshgrid(0.5:(images(i0).cam.imgsz(1) - 0.5), 0.5:(images(i0).cam.imgsz(2) - 0.5));
  du = interp2(u, v, flow{i}.Vx, uv0(:, 1), uv0(:, 2), '*cubic');
  dv = interp2(u, v, flow{i}.Vy, uv0(:, 1), uv0(:, 2), '*cubic');

  %%
  % Calculate ending positions (i).
  uv = uv0 + [du, dv];
  xyz = images(i).cam.invproject(uv);
  for i_pt = 1:size(xyz, 1)
    xyz(i_pt, :) = GLACIER_DEM.sample_ray_tri(images(i).cam.xyz, xyz(i_pt, :), true);
  end

  xyz = images(i0).cam.invproject(uv);
  for i_pt = 1:size(xyz, 1)
    xyz(i_pt, :) = GLACIER_DEM.sample_ray_tri(images(i).cam.xyz, xyz(i_pt, :), true);
  end

  %%
  % Store in data structure.
  motion(i0).t0 = images(i0).date_num;
  motion(i0).t = images(i).date_num;
  motion(i0).uv0 = uv0;
  motion(i0).uv = uv;
  motion(i0).xyz0 = gxyz;
  motion(i0).xyz = xyz;
end

%%
% Visualize results.
for i0 = 7%1:(length(images) - 1)
  [dx, dy] = deal(100);

  %%
  % Plot glacier motion on map.
  %
  %  TODO: Incorporate into DEM.plot function
  %
  figure();
  showimg(HORIZON_DEM.x, HORIZON_DEM.y, hillshade(HORIZON_DEM.Z, HORIZON_DEM.x, HORIZON_DEM.y));
  hold on
  ddays = motion(i0).t - motion(i0).t0;
  v = sqrt(sum((motion(i0).xyz - motion(i0).xyz0).^2, 2)) / ddays;
  [X, Y, V] = pts2grid(motion(i0).xyz0(:, 1), motion(i0).xyz0(:, 2), v, dx, dy);
  alphawarp(X, Y, V, 1);
  % caxis([0 100]);
  colormap jet;
  colorbar
  [~, ~, DX] = pts2grid(motion(i0).xyz0(:, 1), motion(i0).xyz0(:, 2), motion(i0).xyz(:, 1) - motion(i0).xyz0(:, 1), dx, dy);
  [~, ~, DY] = pts2grid(motion(i0).xyz0(:, 1), motion(i0).xyz0(:, 2), motion(i0).xyz(:, 2) - motion(i0).xyz0(:, 2), dx, dy);
  % s = 1;
  % quiver(X, Y, s * DX / ddays, s * DY / ddays, 0, 'k');

  %%
  % Plot glacier motion on image.
  figure();
  ind = 1:size(motion(i0).uv, 1);
  % ind(isnan(v(ind))) = [];
  % ind(v(ind) > 10) = [];
  images(i0).plot();
  hold on
  s = 5;
  quiver(motion(i0).uv0(ind, 1), motion(i0).uv0(ind, 2), s * (motion(i0).uv(ind, 1) - motion(i0).uv0(ind, 1)), s * (motion(i0).uv(ind, 2) - motion(i0).uv0(ind, 2)), 0, 'r');
  % colors = jet(101);
  % scatter(motion(i0).uv0(ind, 1), motion(i0).uv0(ind, 2), 20, colors(ceil(v(ind) * 10), :), 'markerFaceColor', 'flat');
  % caxis([0 10]);
  % colormap jet;
  % colorbar
end
