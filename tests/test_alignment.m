%% Tests of image alignment
% Synthesize a new image from a reference assuming a known rotation, estimate the relative rotation from feature matches, and project the image into the reference to test the alignment.

%% Prepare test environment

%%
% Paths
IMGRAFT_PATH = '.';
PROJECT_DIR = fullfile('tests');
IMG_PATH = fullfile(PROJECT_DIR, 'AK10b_20141012_213336.JPG');
SVG_LAND = 'land';

%%
% MATLAB environment
format shortG;
addpath(genpath(IMGRAFT_PATH));

%%
% Images and cameras
cam = Camera();
images = Image(IMG_PATH, cam);
for i = 1:length(images)
  if isfield(images(i).svg, SVG_LAND)
    for j = fieldnames(images(i).svg.(SVG_LAND))'
      images(i).fixedpolys{end + 1} = images(i).svg.(SVG_LAND).(j{1});
    end
  end
end
img0 = images(1);
img0.I = img0.read();
distimg0 = img0;
distimg0.cam.c = distimg0.cam.c + [4, 4];
distimg0.cam.k = [0.2, 0.1];
distimg0.cam.p = [0.001];
distimg0.cam.plot_distortion();

%% Test 1: No rotation, Ideal camera, Project only [PASSED]

%%
% Project image into itself
newimg = img0;
pI = newimg.project(img0.cam);

%%
% Test that images are equal
isequal(img0.I, pI)

%% Test 2: No rotation, Distorted camera, Project only [PASSED]

%%
% Project image into itself
newimg = distimg0;
pI = newimg.project(distimg0.cam);

%%
% Test that images are equal.
% (Leaving out edge pixels)
isequal(img0.I(2:end-1, 2:end-1, :), pI(2:end-1, 2:end-1, :))

%% Test 3: No rotation, Ideal camera [PASSED]

%%
% Solve for camera motion
[newimg, fit, F] = solve_rotation(img0, img0, 20);

%%
% Check rotation errors (degrees)
newimg.cam.viewdir - img0.cam.viewdir

%%
% Check reprojection errors (pixels)
fit.rmse

%%
% Project image to original
pI = newimg.project(img0.cam);

%%
% Test that images are equal (leaving out edge pixels)
isequal(img0.I(2:end-1, 2:end-1, :), pI(2:end-1, 2:end-1, :))

%% Test 4: No rotation, Distorted camera [PASSED]

%%
% Solve for camera motion
[newimg, fit, F] = solve_rotation(distimg0, distimg0, 20);

%%
% Check rotation errors (degrees)
newimg.cam.viewdir - distimg0.cam.viewdir

%%
% Check reprojection errors (pixels)
fit.rmse

%%
% Project image to original
pI = newimg.project(distimg0.cam);

%%
% Test that images are equal (leaving out edge pixels)
isequal(img0.I(2:end-1, 2:end-1, :), pI(2:end-1, 2:end-1, :))

%% Test 5: Known rotation, Ideal camera [PASSED]

%%
% Generate synthetic image from known rotation
viewdir = [0.1, 0.1, 0.1];
img = img0;
img.cam.viewdir = viewdir;
img.I = img0.project(img.cam);
figure();
imshow(img.I);

%%
% Solve for camera motion
[newimg, fit, F] = solve_rotation(img0, img, 20);

%%
% Check rotation errors (degrees)
newimg.cam.viewdir - img.cam.viewdir

%%
% Check reprojection errors (pixels)
fit.rmse

%%
% Project synthetic image to original
%
% Alignment is ~perfect (although softening in the projected image is apparent).
newimg.I = img.I;
pI = newimg.project(img0.cam);
% flip_images(img0.I, pI);

%% Test 6: Known rotation, Distorted camera [PASSED]

%%
% Generate synthetic image from known rotation
viewdir = [0.1, 0.1, 0.1];
img = distimg0;
img.cam.viewdir = viewdir;
img.I = distimg0.project(img.cam);
figure();
imshow(img.I);

%%
% Solve for camera motion
[newimg, fit, F] = solve_rotation(distimg0, img, 20);

%%
% Check rotation errors (degrees)
newimg.cam.viewdir - img.cam.viewdir

%%
% Check reprojection errors (pixels)
fit.rmse

%%
% Project synthetic image to original
%
% Alignment is ~perfect (although softening in the projected image is apparent).
newimg.I = img.I;
pI = newimg.project(distimg0.cam);
% flip_images(distimg0.I, pI);

%%
% Visualize effect of wrong distortion
errimg0 = distimg0;
errimg0.cam.k = errimg0.cam.k * 0.5;
errimg = newimg;
errimg.cam.k = errimg0.cam.k;
pI = errimg.project(errimg0.cam);
% Use optical flow
obj = opticalFlowFarneback('NumPyramidLevels', 3, 'NeighborhoodSize', 3, 'FilterSize', 25);
obj.estimateFlow(rgb2gray(errimg0.I));
flow = obj.estimateFlow(rgb2gray(pI));
% Plot results
figure();
imshow(errimg0.I / 2);
hold on;
plot(flow, 'Decimation', [20, 20], 'Scale', 20);
