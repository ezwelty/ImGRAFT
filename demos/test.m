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

%% Format fixed (land) and free (glacier) polygons
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
fit
for i = 1:length(anchors)
  figure();
  imshow(anchors(i).read());
  hold on
  % Point errors
  duv = anchors(i).cam.projerror(uv{i}, xyz{i});
  plot(uv{i}(:, 1), uv{i}(:, 2), 'g*');
  s = 1; quiver(uv{i}(:, 1), uv{i}(:, 2), s * duv(:, 1), s * duv(:, 2), 0, 'r');
  % Line errors
  plot(luv{i}(:, 1), luv{i}(:, 2), 'g.');
  duv = anchors(i).cam.projerror_lines(luv{i}, lxyz{i});
  s = 1; quiver(luv{i}(:, 1), luv{i}(:, 2), s * duv(:, 1), s * duv(:, 2), 0, 'r');
  for j = 1:length(lxyz{i})
    pluv = anchors(i).cam.project(lxyz{i}{j});
    plot(pluv(:, 1), pluv(:, 2), 'y-');
  end
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

%% --- Calibrate images from anchors ---

rmse = [];
for i = find(~is_anchor)
  I = images(i).read();
  i0 = images(i).anchor;
  I0 = images(i0).read;
  % Generate grid of points in land polygons
  % TODO: Move to outside loop
  gdu = 50; gdv = 50;
  mpu = []; mpv = []; mdu = []; mdv = [];
  
  for j = 1:length(images(i0).land)

    [du, dv] = templatematch(I0, I, pu, pv, 'templatewidth', gdu, 'searchwidth', 4 * gdu);
    mpu = [mpu ; pu]; mpv = [mpv ; pv]; mdu = [mdu ; du]; mdv = [mdv ; dv];
  end
  % Plot matches
  figure
  imshow(I0 / 1.5); hold on;
  s = 2; quiver(mpu, mpv, s * mdu, s * mdv, 0, 'y')
  % Filter matches with RANSAC Fundamental Matrix
  nonans = ~(isnan(mdu) | isnan(mdv));
  mpu = mpu(nonans); mpv = mpv(nonans);
  mdu = mdu(nonans); mdv = mdv(nonans);
  A = images(i0).cam.image2camera([mpu, mpv]);
  B = images(i0).cam.image2camera([mpu + mdu, mpv + mdv]);
  % TODO: Tricky to set threshold (re-express in pixels?)
  [F, inliersF] = ransacfitfundmatrix(A', B', 0.0000005);
  % [H, inliersH] = ransacfithomography(A', B', 0.00000005);
  % Plot filtered matches
  figure
  imshow(I0 / 1.5); hold on;
  s = 2; quiver(mpu, mpv, s * mdu, s * mdv, 0, 'r')
  s = 2; quiver(mpu(inliersF), mpv(inliersF), s * mdu(inliersF), s * mdv(inliersF), 0, 'y')
  % s = 100; quiver(mpu(inliersH), mpv(inliersH), s * mdu(inliersH), s * mdv(inliersH), 0, 'y')
  % Orient image
  [newcam, fit] = images(i0).cam.optimizeR([mpu(inliersF), mpv(inliersF)], [mpu(inliersF) + mdu(inliersF), mpv(inliersF) + mdv(inliersF)]);
  % [newcam, fit] = images(i0).cam.optimizeR([mpu(inliersH), mpv(inliersH)], [mpu(inliersH) + mdu(inliersH), mpv(inliersH) + mdv(inliersH)])
  rmse(i) = fit.rmse;
  images(i).cam = newcam;
  % Plot transformed traces
  figure()
  imshow(I / 1.5), hold on
  for j = fieldnames(images(i0).svg.horizon)'
    puv = images(i).cam.project(images(i0).cam.invproject(images(i0).svg.horizon.(j{1})), true);
    plot(puv(:, 1), puv(:, 2), 'r-')
  end
  for j = fieldnames(images(i0).svg.coast)'
    puv = images(i).cam.project(images(i0).cam.invproject(images(i0).svg.coast.(j{1})), true);
    plot(puv(:, 1), puv(:, 2), 'r-')
  end
  % Transform glacier polygons
  images(i).glacier = {};
  for j = 1:length(images(i0).glacier)
    images(i).glacier{j} = images(i).cam.project(images(i0).cam.invproject(images(i0).glacier{j}), true);
    % plot(images(i).glacier{j}(:, 1), images(i).glacier{j}(:, 2), 'b-')
  end
end

%% Visualize motion correction
i0 = find(is_anchor, 1);
scale = 0.2;
cam0 = images(i0).cam.resize(scale);
[x0, y0] = meshgrid(1:cam0.imgsz(1), 1:cam0.imgsz(2));
Xi0 = [x0(:), y0(:)];
imwrite(imresize(rgb2gray(imread(images(i0).path)), flip(cam0.imgsz)), [num2str(i0) '-anchor.jpg']);
for i = find(~is_anchor)
  cam = images(i).cam.resize(scale);
  % Project reference grid to new image
  Xi = cam.project(cam0.invproject(Xi0), true);
  % Interpolate image at points
  I = imresize(rgb2gray(imread(images(i).path)), flip(cam.imgsz));
  Zi = interp2(double(I), Xi(:, 1), Xi(:, 2));
  I0 = uint8(reshape(Zi, flip(cam0.imgsz)));
  imwrite(I0, [num2str(i) '.jpg']);
end
