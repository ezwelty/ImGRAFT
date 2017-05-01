addpath(genpath('.'))

%% --- Constants ---

% Paths
DATA_DIR = fullfile('demos', 'data');
IMG_PATH = fullfile(DATA_DIR, 'images', 'test', '*.JPG');
DEM_PATH = fullfile(DATA_DIR, 'dems', '20090827_2m.tif');
OUT_DIR = fullfile('demos', 'test');

% Scalars
IMG_SCALE = 0.25;
DEM_SCALE = 0.25;
DEM_DISTANCE = 20e3; % m
DEM_SMOOTH_WIDTH = 300; % m
COAST_HAE = 17; % m
CAM_ID = lower('AK10');
SVG_GCP = 'gcp';
SVG_COAST = 'coast';
SVG_HORIZON = 'horizon';
SVG_FIXED_AREAS = {'land'};
SVG_MOVING_AREAS = {'glacier'};
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
% FIXME: Won't work if camera missing initial view direction
box = images(1).cam.viewbox(DEM_DISTANCE);
dem = DEM(DEM_PATH).crop(box(:, 1), box(:, 2), [0 Inf]).resize(DEM_SCALE).build();
smdem = dem.smooth(DEM_SMOOTH_WIDTH).build();

%% Load Horizon
% FIXME: Slow
% TODO: Cache results by CAM_ID
% TODO: Traverse grid in concentric circles or triangular wedges
% HORIZON = dem.horizon(cam.xyz, 0:0.1:359);
HORIZON = images(1).cam.horizon(dem, 0.1);

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

%% --- Calibrate anchor images ----
is_anchor = find(arrayfun(@(img) any(isfield(img.svg, {SVG_GCP, SVG_COAST, SVG_HORIZON})), images));

%% Gather calibration data
% cams: {cam, ...}, gcp_uv: {[], ...}, gcp_xyz: {[], ...}, gcl_uv: {[], ...}, gcl_xyz: {{[], ...}, ...}
cams = {images(is_anchor).cam};
gcp_uv = arrayfun(@(img) img.gcp.uv, images(is_anchor), 'UniformOutput', false);
gcp_xyz = arrayfun(@(img) img.gcp.xyz, images(is_anchor), 'UniformOutput', false);
gcl_uv = arrayfun(@(img) img.gcl.uv, images(is_anchor), 'UniformOutput', false);
gcl_xyz = arrayfun(@(img) img.gcl.xyz, images(is_anchor), 'UniformOutput', false);

%% Optimize
% FIXME: Solutions using lines are unstable!
% Adjust xtol (see camera.optimizeCams)
% HACK: Converge on full model:
[newcams, fit] = Camera.optimizeCams(cams, gcp_uv, gcp_xyz, {'viewdir'})
% [newcams, fit] = Camera.optimizeCams(cams, {[]}, {[]}, {{'viewdir'}}, {}, gcl_uv, gcl_xyz)
% [newcams, fit] = Camera.optimizeCams(newcams, gcp_uv, gcp_xyz, {{'viewdir', 'f'}}, {}, gcl_uv, gcl_xyz)

% Then, iterate as points (slow, but effective):
% FIXME: SLOW (but effective)
% (for best results, set dmax threshold)
dmax = 10;
freeparams = {{'viewdir', 'f', 'k', [1 2]}};
freeparams = {{'viewdir'}};
for i = 1:50
  previous_cams = newcams;
  previous_rmse = fit.rmse;
  duv = newcams{1}.projerror_nearest(gcl_xyz{1}, gcl_uv{1});
  lxyz2 = newcams{1}.xyz + newcams{1}.invproject(gcl_uv{1} + duv) * 1000;
  [newcams, fit] = Camera.optimizeCams(newcams{1}, [gcp_uv{1}; gcl_uv{1}], [gcp_xyz{1}; lxyz2], freeparams, [], [], [], dmax)
  % [newcams, fit] = cams{1}.optimizeR(luv{1} + duv, luv{1}), newcams = {newcams};
  if fit.rmse > previous_rmse
    newcams = previous_cams;
    break
  end
end

%% Plot results
figure()
for i = 1:length(newcams)
  imshow(imread(images(is_anchor_ind(i)).path)); hold on
  % Point errors
  duv = newcams{i}.projerror(xyz{i}, uv{i});
  plot(uv{i}(:, 1), uv{i}(:, 2), 'g*');
  s = 1; quiver(uv{i}(:, 1), uv{i}(:, 2), s * duv(:, 1), s * duv(:, 2), 0, 'y');
  % Line errors
  plot(luv{i}(:, 1), luv{i}(:, 2), 'g*');
  duv = newcams{i}.projerror_nearest(lxyz{i}, luv{i});
  s = 1; quiver(luv{i}(:, 1), luv{i}(:, 2), s * duv(:, 1), s * duv(:, 2), 0, 'y');
  for j = 1:length(lxyz{i})
    pluv = newcams{i}.project(lxyz{i}{j});
    plot(pluv(:, 1), pluv(:, 2), 'r-');
  end
  if length(newcams) > 1
    pause();
  end
end
%% Save results
for i = 1:length(newcams)
  images(is_anchor_ind(i)).cam = newcams{i};
end
