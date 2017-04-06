%%% AK10_20090825

% /applications/matlab_r2016b.app/bin/matlab -display :0.0 -nodesktop -display :0.0 -nodesktop
cd ~/sites/ImGRAFT/
addpath(genpath('.'))
project_root = fullfile('ezw', 'AK10_20090825-tide');
% data_root = fullfile('/', 'volumes', 'science', 'data', 'columbia');
data_root = fullfile('ezw', 'data');

%% Load images
c = camera('xyz', [499211.336 6783755.954 478.96], 'viewdir', [-165 -10 -4]);
img_files = dir(fullfile(project_root, 'images', '*.JPG'));
img_paths = arrayfun(@(f) fullfile(f.folder, f.name), img_files, 'UniformOutput', false);
images = loadimages(img_paths, c);

%% Load control data
GCP = readtable(fullfile('ezw', 'gcp.tsv'), 'filetype', 'text', 'delimiter', '\t');
GCP.Properties.RowNames = lower(GCP.name);
svg_files = dir(fullfile(project_root, 'images', '*.svg'));
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

%% Load DEM
[Z, ~, bbox] = geotiffread(fullfile(data_root, 'dem', '20090827_2m.tif'));
gdem = DEM(double(Z), bbox(:, 1), flip(bbox(:, 2)));
viewbox = images(1).cam.viewbox(20 * 1e3);
gcdem = gdem.crop(viewbox(:, 1), viewbox(:, 2), [0 Inf]);
gsmdem = gcdem.resize(0.25);
%% Crevasse-filling
% per Messerli & Grinsted 2015
wm = 300; % m
w = 2 .* round((wm / gsmdem.dx + 1) / 2) - 1; % nearest odd integer
Zg = imgaussfilt(gsmdem.Z, 6, 'FilterSize', w, 'Padding', 'replicate', 'FilterDomain', 'spatial');
disk = fspecial('disk', w);
a = nanstd(gsmdem.Z(:) - Zg(:));
Zd = imfilter(exp((gsmdem.Z - Zg) ./ a), disk, 'replicate');
Zf = a * log(Zd) + Zg;
gsmdem.Z = Zf;
% gsmdemf = DEM(gsmdem.Z, gsmdem.xlim, gsmdem.ylim);
% gsmdemf.Z = Zf;
% plot(gsmdem.Z(:, 600), 'k--'); hold on;
% plot(Zg(:, 600), 'r-');
% plot(Zf(:, 600), 'b-');
%% Draw circle around camera
[x0, y0] = gsmdem.xy2ind(images(1).cam.xyz);
r = round(100 / gsmdem.dx);
[xc, yc] = getmidpointcircle(x0, y0, r);
% Fill circle
ind = [];
y = unique(yc); yin = ~(y < 1 | y > gsmdem.ny);
for yi = reshape(y(yin), 1, [])
  xb = xc(yc == yi);
  xi = max(min(xb), 1):min(max(xb), gsmdem.nx);
  ind((end + 1):(end + length(xi))) = sub2ind([gsmdem.ny, gsmdem.nx], repmat(yi, 1, length(xi)), xi);
end
% Apply to DEM
gsmdem.Z(ind) = NaN;
% gsmdem.plot(2);

%% Calibrate anchor images
% (images with svg.gcp field)
is_anchor = arrayfun(@(img) eval('isfield(img.svg, ''gcp'')'), images)';
is_anchor_ind = find(is_anchor);
cams = {images(is_anchor).cam.idealize};
%% Points
gcp = arrayfun(@(img) cell2mat(struct2cell(img.svg.gcp)), images(is_anchor), 'UniformOutput', false);
uv = cellfun(@(x) x(:, 1:2), gcp, 'UniformOutput', false);
xyz = cellfun(@(x) x(:, 3:5), gcp, 'UniformOutput', false);
%% Lines
% Horizon
hxyz = {images(is_anchor_ind(1)).cam.horizon(gcdem, 0.1)};
% Coastline
shp = shaperead(fullfile(data_root, 'outline', 'coastline.shp'));
cxyz = arrayfun(@(s) [s.X' s.Y' repmat(17, length(s.X), 1)], shp, 'UniformOutput', false);
lxyz = {[hxyz; cxyz]};
% Traces
luv = {};
for i = 1:length(is_anchor_ind)
  traces = {};
  if isfield(images(is_anchor_ind(i)).svg, 'horizon')
    traces = [traces; struct2cell(images(is_anchor_ind(i)).svg.horizon)];
  end
  if isfield(images(is_anchor_ind(i)).svg, 'coast')
    traces = [traces; struct2cell(images(is_anchor_ind(i)).svg.coast)];
  end
  for j = 1:length(traces)
    temp = traces{j} + 0.5;
    l = polylineLength(temp);
    traces{j} = resamplePolyline(temp, round(l / 10));
  end
  luv{i} = cell2mat(traces);
end
%% Optimize
% FIXME: Solutions using lines are unstable!
% Adjust xtol (see camera.optimizeCams)
% HACK: Converge on full model:
[newcams, fit] = camera.optimizeCams(cams, xyz, uv, {'viewdir'})
[newcams, fit] = camera.optimizeCams(newcams, {[]}, {[]}, {{'viewdir'}}, {}, lxyz, luv)
[newcams, fit] = camera.optimizeCams(newcams, xyz, uv, {{'viewdir', 'f'}}, {}, lxyz, luv)
% Then, iterate as points (slow, but effective):
% (for best results, set dmax threshold)
dmax = 10;
freeparams = {{'viewdir', 'f', 'k', [1 2]}};
for i = 1:50
  previous_cams = newcams;
  previous_rmse = fit.rmse;
  duv = newcams{1}.projerror_nearest(lxyz{1}, luv{1});
  lxyz2 = newcams{1}.xyz + newcams{1}.invproject(luv{1} + duv) * 1000;
  [newcams, fit] = camera.optimizeCams(newcams{1}, [xyz{1}; lxyz2], [uv{1}; luv{1}], freeparams, [], [], [], dmax)
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

%% Assign images to nearest anchor
anchor_dates = [images(is_anchor).date];
for i = find(~is_anchor)
  [~, i_min] = min(abs(anchor_dates - images(i).date));
  images(i).anchor = images(is_anchor_ind(i_min)).path;
end

%% Segment anchor images
for i = find(is_anchor)
  images(i).glacier = {};
  for j = fieldnames(images(i).svg.glacier)'
    images(i).glacier{end + 1} = images(i).svg.glacier.(j{1});
  end
  images(1).land = {};
  for j = fieldnames(images(i).svg.land)'
    images(i).land{end + 1} = images(i).svg.land.(j{1});
  end
  % % Plot
  % figure
  % imshow(imread(images(i).path));
  % for j = fieldnames(images(i).svg.glacier)'
  %   patch(images(i).svg.glacier.(j{1})(:, 1), images(i).svg.glacier.(j{1})(:, 2), 1, 'Facecolor', 'b');
  % end
  % for j = fieldnames(images(i).svg.land)'
  %   patch(images(i).svg.land.(j{1})(:, 1), images(i).svg.land.(j{1})(:, 2), 1, 'Facecolor', 'y');
  % end
  % alpha(0.5);
end

%% Calibrate images from assigned anchors
image_paths = {images.path};
rmse = [];
for i = find(~is_anchor)
  i0 = find(strcmp(images(i).anchor, image_paths));
  I0 = imread(images(i0).path);
  I = imread(images(i).path);

  % Generate grid of points in land polygons
  % FIXME: Pre-process anchor images first?
  gdu = 50; gdv = 50;
  mpu = []; mpv = []; mdu = []; mdv = [];
  for j = 1:length(images(i0).land)
    u = min(images(i0).land{j}(:, 1)):gdu:max(images(i0).land{j}(:, 1));
    v = min(images(i0).land{j}(:, 2)):gdv:max(images(i0).land{j}(:, 2));
    [pu pv] = meshgrid(u, v);
    % FIXME: Remove points too close to edge?
    [in on] = inpolygon(pu, pv, images(i0).land{j}(:, 1), images(i0).land{j}(:, 2));
    pu = pu(in & ~on); pv = pv(in & ~on);
    [du, dv] = templatematch(I0, I, pu, pv, 'templatewidth', gdu, 'searchwidth', 2 * gdu);
    mpu = [mpu ; pu]; mpv = [mpv ; pv]; mdu = [mdu ; du]; mdv = [mdv ; dv];
  end

  % % Plot matches
  % figure
  % imshow(I0 / 1.5); hold on;
  % s = 2; quiver(mpu, mpv, s * mdu, s * mdv, 0, 'y')

  % Filter matches with RANSAC Fundamental Matrix
  nonans = ~(isnan(mdu) | isnan(mdv));
  mpu = mpu(nonans); mpv = mpv(nonans);
  mdu = mdu(nonans); mdv = mdv(nonans);
  A = images(i0).cam.image2camera([mpu, mpv]);
  B = images(i0).cam.image2camera([mpu + mdu, mpv + mdv]);
  [F, inliersF] = ransacfitfundmatrix(A', B', 0.0000005); % TODO: tricky to set threshold (re-express in pixels?)
  % [H, inliersH] = ransacfithomography(A', B', 0.00000005); % TODO: tricky to set threshold (re-express in pixels?)

  % % Plot filtered matches
  % figure
  % imshow(I0 / 1.5); hold on;
  % s = 2; quiver(mpu, mpv, s * mdu, s * mdv, 0, 'r')
  % s = 2; quiver(mpu(inliersF), mpv(inliersF), s * mdu(inliersF), s * mdv(inliersF), 0, 'y')
  % % s = 100; quiver(mpu(inliersH), mpv(inliersH), s * mdu(inliersH), s * mdv(inliersH), 0, 'y')

  % Orient image
  [newcam, fit] = images(i0).cam.optimizeR([mpu(inliersF), mpv(inliersF)], [mpu(inliersF) + mdu(inliersF), mpv(inliersF) + mdv(inliersF)]);
  % [newcam, fit] = images(i0).cam.optimizeR([mpu(inliersH), mpv(inliersH)], [mpu(inliersH) + mdu(inliersH), mpv(inliersH) + mdv(inliersH)])
  rmse(i) = fit.rmse;
  images(i).cam = newcam;

  % % Plot transformed traces
  % figure()
  % imshow(I / 1.5), hold on
  % for j = fieldnames(images(i0).svg.horizon)'
  %   puv = images(i).cam.project(images(i0).cam.invproject(images(i0).svg.horizon.(j{1})), true);
  %   plot(puv(:, 1), puv(:, 2), 'r-')
  % end
  % for j = fieldnames(images(i0).svg.coast)'
  %   puv = images(i).cam.project(images(i0).cam.invproject(images(i0).svg.coast.(j{1})), true);
  %   plot(puv(:, 1), puv(:, 2), 'r-')
  % end

  % Transform glacier polygons
  images(i).glacier = {};
  for j = 1:length(images(i0).glacier)
    images(i).glacier{j} = images(i).cam.project(images(i0).cam.invproject(images(i0).glacier{j}), true);
    % plot(images(i).glacier{j}(:, 1), images(i).glacier{j}(:, 2), 'b-')
  end
end

%% Visualize motion correction
i0 = find(is_anchor, 1);
scale = 0.5;
cam0 = images(i0).cam.resize(scale);
[x0, y0] = meshgrid(1:cam0.imgsz(1), 1:cam0.imgsz(2));
Xi0 = [x0(:), y0(:)];
imwrite(imresize(rgb2gray(imread(images(i0).path)), flip(cam0.imgsz)), fullfile(project_root, [num2str(i0) '-anchor.jpg']));
for i = find(~is_anchor)
  cam = images(i).cam.resize(scale);
  % Project reference grid to new image
  Xi = cam.project(cam0.invproject(Xi0), true);
  % Interpolate image at points
  I = imresize(rgb2gray(imread(images(i).path)), flip(cam.imgsz));
  Zi = interp2(double(I), Xi(:, 1), Xi(:, 2));
  I0 = uint8(reshape(Zi, flip(cam0.imgsz)));
  imwrite(I0, fullfile(project_root, [num2str(i) '.jpg']));
end

%% Glacier points (v1)
i = find(is_anchor, 1);
gxyz = images(i).cam.invproject(images(i).glacier{1}, gsmdem);
gxyz(any(isnan(gxyz), 2), :) = [];
gdx = 10; gdy = 10;
x = min(gxyz(:, 1)):gdx:max(gxyz(:, 1));
y = min(gxyz(:, 2)):gdy:max(gxyz(:, 2));
[gx gy] = meshgrid(x, y);
gv = false(size(gx));
% Discard points outside glacier polygon
% FIXME: Remove points too close to edge?
[in on] = inpolygon(gx, gy, gxyz(:, 1), gxyz(:, 2));
gv(in & ~on) = true;
% Discard points not visible to camera
g = cell2mat(cellfun(@gsmdem.sample, num2cell([gx(gv(:)) gy(gv(:))], 2), 'UniformOutput', false));
% gz = nan(size(gx));
% gz(gv(:)) = g(:, 3);
gpts = images(i).cam.invproject(images(i).cam.project(g), gsmdem);
gv(gv(:)) = sqrt(sum((gpts - g).^2, 2)) < 1;
% Discard points in mostly invisible regions
% gv = bwmorph(gv, 'clean');
gv = bwmorph(gv, 'open');
gv = bwmorph(gv, 'shrink', 1);
% Final point selection
gpts = cell2mat(cellfun(@gsmdem.sample, num2cell([gx(gv(:)) gy(gv(:))], 2), 'UniformOutput', false));
gsmdem.plot(2); hold on; plot(gxyz(:, 1), gxyz(:, 2), 'y-'); plot(g(:, 1), g(:, 2), 'y*'); plot(gpts(:, 1), gpts(:, 2), 'r*');

% (v2: dense pixel grid)
i0 = find(is_anchor, 1);
% Generate grid of points in glacier polygons
gdu = 50; gdv = 50;
p0 = [];
for j = 1:length(images(i0).glacier)
  u = min(images(i0).glacier{j}(:, 1)):gdu:max(images(i0).glacier{j}(:, 1));
  v = min(images(i0).glacier{j}(:, 2)):gdv:max(images(i0).glacier{j}(:, 2));
  [pu pv] = meshgrid(u, v);
  [in on] = inpolygon(pu, pv, images(i0).glacier{j}(:, 1), images(i0).glacier{j}(:, 2));
  pu = pu(in & ~on); pv = pv(in & ~on);
  p0 = [p0; pu pv];
end
gpts = images(i0).cam.invproject(p0, gsmdem);

% (v3: benchmark locations) +/-
V = importdata(fullfile(data_root, 'velocity', '20090803_20090827_100m.tab'));
% x0, y0, x1, x2 | v (m/day), h (degrees cw from north)
% Vb = [V.data(:, 1) - V.data(:, 3) / 2, V.data(:, 2) - V.data(:, 4) / 2, V.data(:, 1) + V.data(:, 3) / 2, V.data(:, 2) + V.data(:, 4) / 2, V.data(:, 8), V.data(:, 9)];
% cx, cy | vx, vy (m/day) | h (degrees cw from north)
Vb = [V.data(:, 1:2), V.data(:, 5:6), V.data(:, 9)];
has_nan = any(isnan(Vb), 2);
out_of_bounds = Vb(:, 1) < gsmdem.min(1) | Vb(:, 1) > gsmdem.max(1) | Vb(:, 2) < gsmdem.min(2) | Vb(:, 2) > gsmdem.max(2);
Vb(has_nan | out_of_bounds, :) = [];
% gsmdem.plot(2); hold on; quiver(Vb(:, 1), Vb(:, 2), Vb(:, 3), Vb(:, 4));
gpts0 = nan(size(Vb, 1), 3);
for i = 1:size(gpts0, 1)
  gpts0(i, :) = gsmdem.sample(Vb(i, 1:2));
end
[Xi, infront] = images(1).cam.project(gpts0);
inframe = images(1).cam.inframe(Xi);
gpts1 = images(1).cam.invproject(Xi, gsmdem);
% plot(gpts0(:, 1), gpts0(:, 2), 'ko', gpts1(:, 1), gpts1(:, 2), 'r*')
d = sqrt(sum((gpts1 - gpts0).^2, 2));
valid = d < 1 & infront & inframe;
Vb = Vb(valid, :);
gpts = gpts0(valid, :);
% gsmdem.plot(2); hold on; plot(gpts(:, 1), gpts(:, 2), 'y.')
% Xi = images(1).cam.project(gpts);
% imshow(imread(images(1).path) / 1.25); hold on; plot(Xi(:, 1), Xi(:, 2), 'y*')
% (add +/-)
d = 50;
gpts = [gpts(:, 1:2); gpts(:, 1:2) + d; gpts(:, 1:2) - d; gpts(:, 1) + d, gpts(:, 2); gpts(:, 1) - d, gpts(:, 2); gpts(:, 1), gpts(:, 2) + d; gpts(:, 1), gpts(:, 2) - d; gpts(:, 1) - d, gpts(:, 2) + d; gpts(:, 1) + d, gpts(:, 2) - d];
gpts(:, 3) = NaN;
for i = 1:size(gpts, 1)
  gpts(i, :) = gsmdem.sample(gpts(i, 1:2));
end
% gsmdem.plot(2); hold on; plot(gpts(:, 1), gpts(:, 2), 'y.')
% Xi = images(1).cam.project(gpts);
% imshow(imread(images(1).path) / 1.25); hold on; plot(Xi(:, 1), Xi(:, 2), 'y*')

%% Velocities
i_start = 1;
Mi = {};
Mw = {};
for i0 = i_start:(length(images) - 1)
  i = i0 + 1;
  I0 = imread(images(i0).path);
  I = imread(images(i).path);
  p0 = images(i0).cam.project(gpts);
  p = images(i).cam.project(gpts);
  % imshow(I0 / 1.5); hold on; plot(p0(:, 1), p0(:, 2), 'y*');
  [du, dv] = templatematch(I0, I, p0(:, 1), p0(:, 2), 'templatewidth', 40, 'searchwidth', 60, 'initialdu', p(:, 1) - p0(:, 1), 'initialdv', p(:, 2) - p0(:, 2), 'supersample', 2);
  % pXw = images(i).cam.invproject(p0 + [du, dv], gsmdem);
  d = images(i).cam.invproject(p0 + [du, dv]);
  pXw = nan(size(d));
  for j = 1:size(d, 1)
    % X = gsmdem.sample(images(i).cam.xyz, d(j, :), false);
    % [~, imin] = min(sum(bsxfun(@minus, X, gpts(j, :)).^2, 2));
    % pXw(j, :) = X(imin, :);
    if ~any(isnan(d(j, :)))
      pXw(j, :) = gsmdem.sample(images(i).cam.xyz, d(j, :), true);
    end
  end
  Mw{end + 1} = pXw;
  Mi{end + 1} = images(i_start).cam.project(images(i).cam.invproject(p0 + [du, dv]), true);
  % gpts = pXw;
  % figure
  % imshow(I0 / 1.5); hold on;
  % pXi0 = images(i0).cam.project(images(i).cam.invproject(p0 + [du, dv]), true);
  % s = 5; quiver(p0(:, 1), p0(:, 2), s * (pXi0(:, 1) - p0(:, 1)), s * (pXi0(:, 2) - p0(:, 2)), 0, 'y');
  % figure
  % ddays = abs(images(i0).date - images(i).date);
  % gsmdem.plot(2); hold on;
  % s = 5; quiver(Vb(:, 1), Vb(:, 2), s * Vb(:, 3), s * Vb(:, 4), 0, 'k');
  % s = 5; quiver(gpts(:, 1), gpts(:, 2), s * (pXw(:, 1) - gpts(:, 1)) / ddays, s * (pXw(:, 2) - gpts(:, 2)) / ddays, 0, 'r');
end

%% Compute average daily velocity for each point
ddays = diff([images.date]);

% (world)
Xw = reshape(cell2mat(Mw), [size(gpts, 1), 3, length(Mw)]);
dXw = Xw - gpts;
Vw = nan(size(dXw(:, :, 1)));
for i = 1:length(ddays)
  Vw(:, :, i) = dXw(:, :, i) ./ ddays(i);
end
Vwm = nanmean(Vw, 3);
Vwmf = permute(reshape(permute(Vwm, [2 1 3]), [3, size(Vb, 1), 9]), [2 1 3]);
Vwmf = nanmean(Vwmf, 3);
%
figure
gsmdem.plot(2); hold on;
s = 8;
quiver(Vb(:, 1), Vb(:, 2), s * Vb(:, 3), s * Vb(:, 4), 0, 'g', 'linewidth', 2);
dv = Vwmf(:, 1:2) - Vb(:, 3:4); ind = sqrt(sum(dv.^2, 2)) < 100;
quiver(Vb(ind, 1), Vb(ind, 2), s * dv(ind, 1), s * dv(ind, 2), 0, 'r', 'linewidth', 2)
% quiver(Vb(:, 1), Vb(:, 2), s * Vwmf(:, 1), s * Vwmf(:, 2), 0, 'r', 'linewidth', 2);
% quiver(gpts(:, 1), gpts(:, 2), s * Vwm(:, 1), s * Vwm(:, 2), 0, 'r');
% di = 2; quiver(gpts(:, 1), gpts(:, 2), s * dXw(:, 1, di) / ddays(di), s * dXw(:, 2, di) / ddays(di), 0, 'r');
%
figure
sp1 = subaxis(2, 2, 1)
showimg(gsmdem.x, gsmdem.y, hillshade(gsmdem.Z, gsmdem.x, gsmdem.y));
hold on;
[x, y, Z] = pts2grid(Vb(:, 1), Vb(:, 2), sqrt(sum(Vb(:, 3:4).^2, 2)), 100, 100);
[X, Y] = meshgrid(x, y);
alphawarp(X, Y, Z, 1);
caxis([0 20]);
colormap jet;
colorbar
sp2 = subaxis(2, 2, 2)
showimg(gsmdem.x, gsmdem.y, hillshade(gsmdem.Z, gsmdem.x, gsmdem.y));
hold on;
[x, y, Z] = pts2grid(gpts(:, 1), gpts(:, 2), sqrt(sum(Vwm(:, 1:2).^2, 2)), 100, 100);
[X, Y] = meshgrid(x, y);
alphawarp(X, Y, Z, 1);
caxis([0 20]);
colormap jet;
colorbar
sp3 = subaxis(2, 2, 3)
showimg(gsmdem.x, gsmdem.y, hillshade(gsmdem.Z, gsmdem.x, gsmdem.y));
hold on;
[x, y, Z] = pts2grid(Vb(:, 1), Vb(:, 2), Vwm(1:size(Vb, 1), 1) - Vb(:, 3), 100, 100);
[X, Y] = meshgrid(x, y);
alphawarp(X, Y, Z, 1);
caxis([0 20]);
colormap jet;
colorbar
sp4 = subaxis(2, 2, 4)
showimg(gsmdem.x, gsmdem.y, hillshade(gsmdem.Z, gsmdem.x, gsmdem.y));
hold on;
[x, y, Z] = pts2grid(Vb(:, 1), Vb(:, 2), Vwm(1:size(Vb, 1), 2) - Vb(:, 4), 100, 100);
[X, Y] = meshgrid(x, y);
alphawarp(X, Y, Z, 1);
caxis([0 20]);
colormap jet;
colorbar
linkaxes([sp1, sp2, sp3, sp4], 'xy')
%
figure
% (error vs distance)
subaxis(2, 2, 1)
camera_distances = sqrt(sum((Vb(:, 1:2) - images(1).cam.xyz(1:2)).^2, 2));
plot(camera_distances, sqrt(sum(dv.^2, 2)), 'k.')
% (error vs distance - percent)
subaxis(2, 2, 2)
camera_distances = sqrt(sum((Vb(:, 1:2) - images(1).cam.xyz(1:2)).^2, 2));
plot(camera_distances, 100 * sqrt(sum(dv.^2, 2)) ./ sqrt(sum(Vb(:, 3:4).^2, 2)), 'k.')
ylim([0 100]);
% (error vs angle)
subaxis(2, 2, 3)
dz = images(1).cam.xyz(3) - gpts(1:size(Vb, 1), 3);
ray_angles = atand(dz ./ camera_distances);
plot(ray_angles, sqrt(sum(dv.^2, 2)), 'k.')
set(gca,'XDir','reverse')
% (error vs angle - percent)
subaxis(2, 2, 4)
camera_distances = sqrt(sum((Vb(:, 1:2) - images(1).cam.xyz(1:2)).^2, 2));
plot(ray_angles, 100 * sqrt(sum(dv.^2, 2)) ./ sqrt(sum(Vb(:, 3:4).^2, 2)), 'k.')
set(gca,'XDir','reverse')
ylim([0 100]);
%
% (velocity vs time vs tides)
figure
subaxis(2, 1, 1)
times = [images.date];
Vw_temp = permute(reshape(permute(Vw, [2 1 3]), [3, size(Vb, 1), 9, size(Vw, 3)]), [2 1 3 4]);
Vwf = squeeze(nanmean(Vw_temp, 3));
hold on
for i = 1:size(Vwf, 1)
  if (ray_angles(i) > 10 & camera_distances(i) < 3e3)
    plot(times(1:end-1), sqrt(sum(squeeze(Vwf(i, 1:2, :))'.^2, 2)), '-', 'color', [0.5, 0.5, 0.5]);
  else
    % plot(times(1:end-1), sqrt(sum(squeeze(Vwf(i, 1:2, :))'.^2, 2)), 'k-');
  end
end
ylim([0 10])
xlim(range(times))
plot(repmat(times, 2, 1), repmat([0; 1000], 1, length(times)), 'k:')
subaxis(2, 1, 2)
temp = importdata(fullfile(data_root, 'tides.tab'));
tides = struct();
tides.date = temp(:, 5);
tides.height = temp(:, 6);
ind = tides.date > min(times) & tides.date < max(times);
plot(tides.date(ind), tides.height(ind), 'b-'); hold on;
plot(repmat(times, 2, 1), repmat([0; 1000], 1, length(times)), 'k:')
xlim([min(times), max(times)])
ylim(1.1 * [min(tides.height(ind)), max(tides.height(ind))])

% (image)
gpts1_xy = [Vb(:, 1) + Vb(:, 3), Vb(:, 2) + Vb(:, 4)];
gpts1 = nan(size(gpts1_xy, 1), 3);
for i = 1:size(gpts1, 1)
  gpts1(i, :) = gsmdem.sample(gpts1_xy(i, :));
end
p0 = images(i_start).cam.project(gpts(1:size(Vb, 1), :));
p1 = images(i_start).cam.project(gpts1);
%
p1i = images(i_start).cam.project(gpts(1:size(Vb, 1), :) + Vwmf);
%
figure
imshow(imread(images(i_start).path) / 1.25); hold on
s = 5;
quiver(p0(:, 1), p0(:, 2), s * (p1(:, 1) - p0(:, 1)), s * (p1(:, 2) - p0(:, 2)), 0, 'g', 'linewidth', 2);
quiver(p0(:, 1), p0(:, 2), s * (p1i(:, 1) - p1(:, 1)), s * (p1i(:, 2) - p1(:, 2)), 0, 'r', 'linewidth', 2);
% quiver(p0(:, 1), p0(:, 2), s * (p1i(:, 1) - p0(:, 1)), s * (p1i(:, 2) - p0(:, 2)), 0, 'r');
% quiver(p0(:, 1), p0(:, 2), s * Vim(:, 1), s * Vim(:, 2), 0, 'r')
% di = 2; quiver(p0(:, 1), p0(:, 2), s * dXi(:, 1, di) / ddays(di), s * dXi(:, 2, di) / ddays(di), 0, 'r');




%% Test Camera-World conversion
[ortho, ref] = geotiffread(fullfile(data_root, 'ortho', '20090827_2m.tif'));

Xi = nan(1, 3);
Xw = images(i_start).cam.invproject(Xi, gsmdem);

image_fig = figure
imshow(imread(images(i_start).path)); hold on
plot(Xi(:, 1), Xi(:, 2), 'y-', 'XDataSource', 'Xi(:, 1)', 'YDataSource', 'Xi(:, 2)'); hold off
linkdata on

world_fig = figure
% gsmdem.plot(2); hold on;
mapshow(ortho, ref); hold on;
plot(Xw(:, 1), Xw(:, 2), 'y-', 'XDataSource', 'Xw(:, 1)', 'YDataSource', 'Xw(:, 2)'); hold off
linkdata on

figure(image_fig)
Xi = ginput();
Xw = images(i_start).cam.invproject(Xi, gsmdem);

figure(world_fig)
temp = ginput();
temp2 = nan(size(temp, 1), 3);
for i = 1:size(temp2, 1)
  temp2(i, :) = gsmdem.sample(temp(i, :));
end
Xw = temp2;
Xi = images(i_start).cam.project(Xw);




% % previous results
% w = reshape(cell2mat(oldMw), [size(gpts, 1), 3, length(oldMw)]);
% dw = bsxfun(@minus, w, w(:, :, 1));
% dwm = nanmedian(dw(:, :, 2:end), 3);
% s = 5; quiver(gpts(:, 1), gpts(:, 2), s * dwm(:, 1), s * dwm(:, 2), 0, 'g');

% %% Plot each frame in turn
w = reshape(cell2mat(Mw), [size(gpts, 1), 3, length(Mw)]);
dw = bsxfun(@minus, w, w(:, :, 1));
figure
gsmdem.plot(2); hold on;
s = 5; quiver(repmat(V.Xm, size(V.Ym, 1), 1), repmat(V.Ym, 1, size(V.Xm, 2)), s * V.Gdxx, s * V.Gdyy, 0, 'k');
s = 5; h = quiver(gpts(:, 1), gpts(:, 2), s * dw(:, 1, 2), s * dw(:, 2, 2), 0, 'r');
for i = 3:size(dw, 3)
  pause();
  set(h, 'udata', s * dw(:, 1, i), 'vdata', s * dw(:, 2, i));
end

figure
imshow(imread(images(i_start).path)); hold on
s = 3;
% plot(Mi{1}(:, 1), Mi{1}(:, 2), 'g*');
for i = 2:length(Mi)
  quiver(Mi{i - 1}(:, 1), Mi{i - 1}(:, 2), s * (Mi{i}(:, 1) - Mi{i - 1}(:, 1)), s * (Mi{i}(:, 2) - Mi{i - 1}(:, 2)), 0, 'r');
  quiver(Mi{i - 1}(j, 1), Mi{i - 1}(j, 2), s * (Mi{i}(j, 1) - Mi{i - 1}(j, 1)), s * (Mi{i}(j, 2) - Mi{i - 1}(j, 2)), 0, 'y');
end
quiver(p0(:, 1), p0(:, 2), s * (p1(:, 1) - p0(:, 1)), s * (p1(:, 2) - p0(:, 2)), 0, 'g');
%
vx0 = repmat(V.Xm, size(V.Ym, 1), 1); vx0 = vx0(:);
vy0 = repmat(V.Ym, 1, size(V.Xm, 2)); vy0 = vy0(:);
dvx = V.Gdxx(:);
dvy = V.Gdyy(:);
vx1 = vx0 + dvx;
vy1 = vy0 + dvy;
has_nan = isnan(vx0) | isnan(vy0) | isnan(vx1) | isnan(vy1);
out_of_bounds = vx0 < gsmdem.min(1) | vx0 > gsmdem.max(1) | vy0 < gsmdem.min(2) | vy0 > gsmdem.max(2) | vx1 < gsmdem.min(1) | vx1 > gsmdem.max(1) | vy1 < gsmdem.min(2) | vy1 > gsmdem.max(2);
vx0(has_nan | out_of_bounds) = [];
vy0(has_nan | out_of_bounds) = [];
vx1(has_nan | out_of_bounds) = [];
vy1(has_nan | out_of_bounds) = [];
X0 = nan(length(vx0), 3);
X1 = X0;
for i = 1:size(X0, 1)
  X0(i, :) = gsmdem.sample([vx0(i) vy0(i)]);
  X1(i, :) = gsmdem.sample([vx1(i) vy1(i)]);
end
p0 = images(i_start).cam.project(X0);
p1 = images(i_start).cam.project(X1);

figure
gsmdem.plot(2); hold on
s = 3;
quiver(repmat(V.Xm, size(V.Ym, 1), 1), repmat(V.Ym, 1, size(V.Xm, 2)), s * V.Gdxx, s * V.Gdyy, 0, 'k');
% plot(Mw{1}(:, 1), Mw{1}(:, 2), 'g*');
for i = 2:length(Mw)
  directions = images(i_start).cam.invproject(Mi{i});
  temp = nan(size(directions));
  for j = 1:size(directions, 1)
    if all(~isnan(directions(j, :)))
      temp(j, :) = intersectRayPlane(images(i_start).cam.xyz, directions(j, :), [0 0 -1 Mw{1}(j, 3)]);
    end
  end
  quiver(Mw{i - 1}(:, 1), Mw{i - 1}(:, 2), s * (temp(:, 1) - Mw{i - 1}(:, 1)), s * (temp(:, 2) - Mw{i - 1}(:, 2)), 0, 'y');
  % quiver(Mw{i - 1}(:, 1), Mw{i - 1}(:, 2), s * (Mw{i}(:, 1) - Mw{i - 1}(:, 1)), s * (Mw{i}(:, 2) - Mw{i - 1}(:, 2)), 0, 'y');
end
