%%% AK10_20090825

% /applications/matlab_r2016b.app/bin/matlab -display :0.0 -nodesktop -display :0.0 -nodesktop
cd ~/sites/ImGRAFT/
addpath(genpath('.'))
project_root = fullfile('ezw', 'AK10_20090825');
data_root = fullfile('/', 'volumes', 'science-b', 'data', 'columbia');

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
[Z, ~, bbox] = geotiffread(fullfile(data_root, 'dem', '2004 Aerometric', '20090827_2m.tif'));
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
  % figure()
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
  gdu = 100; gdv = 100;
  mpu = []; mpv = []; mdu = []; mdv = [];
  for j = 1:length(images(i0).land)
    u = min(images(i0).land{j}(:, 1)):gdu:max(images(i0).land{j}(:, 1));
    v = min(images(i0).land{j}(:, 2)):gdv:max(images(i0).land{j}(:, 2));
    [pu pv] = meshgrid(u, v);
    % FIXME: Remove points too close to edge?
    [in on] = inpolygon(pu, pv, images(i0).land{j}(:, 1), images(i0).land{j}(:, 2));
    pu = pu(in & ~on); pv = pv(in & ~on);
    [du, dv] = templatematch(I0, I, pu, pv, 'templatewidth', gdu / 2, 'searchwidth', gdu);
    mpu = [mpu ; pu]; mpv = [mpv ; pv]; mdu = [mdu ; du]; mdv = [mdv ; dv];
  end
  % % Plot matches
  % figure
  % imshow(I0 / 1.5); hold on;
  % s = 10; quiver(mpu, mpv, s * mdu, s * mdv, 0, 'y')
  % Filter matches with RANSAC Fundamental Matrix
  nonans = ~(isnan(mdu) | isnan(mdv));
  mpu = mpu(nonans); mpv = mpv(nonans);
  mdu = mdu(nonans); mdv = mdv(nonans);
  A = images(i0).cam.image2camera([mpu, mpv]);
  B = images(i0).cam.image2camera([mpu + mdu, mpv + mdv]);
  [F, inliersF] = ransacfitfundmatrix(A', B', 0.00000005); % TODO: tricky to set threshold (re-express in pixels?)
  % [H, inliersH] = ransacfithomography(A', B', 0.00000005); % TODO: tricky to set threshold (re-express in pixels?)
  % % Plot filtered matches
  % figure
  % imshow(I0 / 1.5); hold on;
  % s = 100; quiver(mpu(inliersF), mpv(inliersF), s * mdu(inliersF), s * mdv(inliersF), 0, 'y')
  % % s = 100; quiver(mpu(inliersH), mpv(inliersH), s * mdu(inliersH), s * mdv(inliersH), 0, 'y')
  % Orient image
  [newcam, fit] = images(i0).cam.optimizeR([mpu(inliersF), mpv(inliersF)], [mpu(inliersF) + mdu(inliersF), mpv(inliersF) + mdv(inliersF)]);
  % [newcam, fit] = images(i0).cam.optimizeR([mpu(inliersH), mpv(inliersH)], [mpu(inliersH) + mdu(inliersH), mpv(inliersH) + mdv(inliersH)])
  rmse(i) = fit.rmse;
  images(i).cam = newcam;
  % Transform glacier polygons
  images(i).glacier = {};
  for j = 1:length(images(i0).glacier)
    images(i).glacier{j} = images(i).cam.project(images(i0).cam.invproject(images(i0).glacier{j}), true);
  end
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
end

%% Glacier points
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

%% Velocities
i_start = i0;
Mi = {images(i_start).cam.project(gpts)};
Mw = {gpts};
for i0 = i_start%:length(images)
  i = i0 + 1;
  I0 = imread(images(i0).path);
  I = imread(images(i).path);
  p0 = images(i0).cam.project(gpts);
  p = images(i).cam.project(gpts);
  % imshow(I0 / 1.5); hold on; plot(p0(:, 1), p0(:, 2), 'y*');
  [du, dv] = templatematch(I0, I, p0(:, 1), p0(:, 2), 'templatewidth', 10, 'searchwidth', 80,'initialdu', p(:, 1) - p0(:, 1), 'initialdv', p(:, 2) - p0(:, 2), 'supersample', 2);
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
  % gsmdem.plot(2); hold on;
  % s = 5; quiver(repmat(V.Xm, size(V.Ym, 1), 1), repmat(V.Ym, 1, size(V.Xm, 2)), s * V.Gdxx, s * V.Gdyy, 0, 'k');
  % s = 5; quiver(gpts(:, 1), gpts(:, 2), s * (pXw(:, 1) - gpts(:, 1)), s * (pXw(:, 2) - gpts(:, 2)), 0, 'r');
end

[ortho, ref] = geotiffread(fullfile(data_root, 'ortho', '2002 Aerometric', '20090827_2m.tif'));

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


%% Compute average vector for each point
V = load(fullfile(data_root, 'velocity', '_yushin', 'airborne', '20090803_20090827_GRD.mat'));
V = V.GRD;
%
w = reshape(cell2mat(Mw), [size(gpts, 1), 3, length(Mw)]);
dw = bsxfun(@minus, w, w(:, :, 1));
dwm = nanmedian(dw(:, :, 2:end), 3);
%
figure
gsmdem.plot(2); hold on;
s = 5; quiver(repmat(V.Xm, size(V.Ym, 1), 1), repmat(V.Ym, 1, size(V.Xm, 2)), s * V.Gdxx, s * V.Gdyy, 0, 'k');
s = 5; quiver(gpts(:, 1), gpts(:, 2), s * dwm(:, 1), s * dwm(:, 2), 0, 'r');
j = find(abs(s * dwm(:, 1) - -24.95) < 1e-1 & abs(s * dwm(:, 2) - 54.4) < 1e-1);

% % previous results
% w = reshape(cell2mat(oldMw), [size(gpts, 1), 3, length(oldMw)]);
% dw = bsxfun(@minus, w, w(:, :, 1));
% dwm = nanmedian(dw(:, :, 2:end), 3);
% s = 5; quiver(gpts(:, 1), gpts(:, 2), s * dwm(:, 1), s * dwm(:, 2), 0, 'g');

% %% Plot each frame in turn
% w = reshape(cell2mat(Mw), [size(gpts, 1), 3, length(Mw)]);
% dw = bsxfun(@minus, w, w(:, :, 1));
% figure
% gsmdem.plot(2); hold on;
% s = 5; quiver(repmat(V.Xm, size(V.Ym, 1), 1), repmat(V.Ym, 1, size(V.Xm, 2)), s * V.Gdxx, s * V.Gdyy, 0, 'k');
% s = 5; h = quiver(gpts(:, 1), gpts(:, 2), s * dw(:, 1, 2), s * dw(:, 2, 2), 0, 'r');
% for i = 3:size(dw, 3)
%   pause();
%   set(h, 'udata', s * dw(:, 1, i), 'vdata', s * dw(:, 2, i));
% end

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




%%%%%%%%%%%

tic
directions = images(i0).cam.invproject(p0);
gpts2 = nan(size(directions));
for i = 1:size(directions, 1)
  gpts2(i, :) = gsmdem.sample(images(i0).cam.xyz, directions(i, :), true, gpts(i, :), 100);
end
toc

tic
directions = images(i0).cam.invproject(p0);
gpts3 = nan(size(directions));
for i = 1:size(directions, 1)
  gpts3(i, :) = intersectRayPlane(images(i0).cam.xyz, directions(i, :), [0 0 -1 gpts(i, 3)]);
end
toc



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


%%%%%%%%%%%%

%% Earth Curvature and refraction!?
% e.g. http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?topicname=how_viewshed_works
% Zactual = Zsurface - 0.87 * Dist^2 / DiamEarth
% Max dz = Dist^2 / DiamEarth (m) = f * Dist / DiamEarth (px)
De = 12740000; % m
dist = [0:1:20] * 1e3; % m
dz_m = dist .^ 2 ./ De % m
dz_px = mean(images(1).cam.f) .* (dist ./ De) % px
plot(dist, dz_m, 'k-', dist, dz_px, 'r-')
