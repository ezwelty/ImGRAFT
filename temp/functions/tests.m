% Tests

%% Focal length (OK)
fmm = 90;
wmm = 35;
imgsize = [400 200];
fpx = fmm2fpx(90, wmm, imgsize);
fmm - fpx2fmm(fpx, wmm, imgsize)

%% Rotations (ROUNDING ERROR)
angles = [30 45 0];
R = ang2rot(angles);
angles - rot2ang(R)

oblangles = [30 45 0];
R = oblang2rot(oblangles);
oblangles - rot2oblang(R)

%% Reprojection (ISSUES)

clear camera
camera.fpx = 3230.7955146593;
camera.c = [1890.13337953846 1279.27135285767];
camera.P = [514632.16 8548770.44 169.53];
camera.R = oblang2rot([95 -10 0]);
sealevel = 0;
seasurface = [0 0 sealevel 1 0 0 0 1 0];
Xi = [0 1500 3888 0 1500 3888 0 1500 3888; 0 0 0 1500 1500 1500 2592 2592 2592];

% Image <=> World <=> Image (ISSUES)

camera.k = [-0.2 0.2 -0.1 -0.1 -0.01]; % OK
% camera.k = [-0.2 0.2 0 0 0]; % OK
% camera.k = [-0.275 0 0 0 0]; % OK
% camera.k = [20 0 0 0 0]; % OK
% camera.k = [0 0 0 0 0]; % OK
% camera.k = [-20 0 0 0 0]; % FAILS
% camera.k = [-1 0 0 0 0]; % FAILS
% camera.k = [1 1 1 1 1]; % FAILS
Xw = image2world(Xi, camera, seasurface);
sum((Xi - world2image(Xw, camera)).^2, 1) < 1

% Image <=> Camera <=> Image
Xc = image2camera(Xi, camera);
sum((Xi - camera2image(Xc, camera)).^2, 1)

% Camera <=> World <=> Camera
Xc = image2camera(Xi, camera);
X = camera2world(Xc, camera, seasurface);
sum((Xc - world2camera(X, camera)).^2, 1)

