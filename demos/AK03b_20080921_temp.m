visible = voxelviewshed(gsmdem.X, gsmdem.Y, gsmdem.Z, images(1).cam.xyz);
ind = gsmdem.xy2ind(gpts);
v = visible(ind);

guv = images(1).cam.project(gpts);
% imshow(imread(images(1).path) / 1.5); hold on; plot(images(1).glacierpolys{1}(:, 1), images(1).glacierpolys{1}(:, 2), 'y-'); plot(guv(:, 1), guv(:, 2), 'r*');
d = images(1).cam.invproject(guv);
% gpts2 = images(1).cam.invproject(guv(1, :), gsmdem)
gpts2 = [];
for j = 1:size(d, 1)
  X = gsmdem.sample(images(1).cam.xyz, d(j, :), false);
  [~, imin] = min(sum(bsxfun(@minus, X, gpts(j, :)).^2, 2));
  gpts2(j, :) = X(imin, :);
end
gpts2 - gpts
plot(gpts(:, 1), gpts(:, 2), 'ko'); hold on;
plot(gpts2(:, 1), gpts2(:, 2), 'r*');


gsmdem.plot(3); hold on
plot3(images(1).cam.xyz(:, 1), images(1).cam.xyz(:, 2), images(1).cam.xyz(:, 3), 'k*');
quiver3(images(1).cam.xyz(:, 1), images(1).cam.xyz(:, 2), images(1).cam.xyz(:, 3), d(i, 1), d(i, 2), d(i, 3), 7e3, 'g');
plot3(gpts(i, 1), gpts(i, 2), gpts(i, 3), 'g*');
plot3(gpts2(i, 1), gpts2(i, 2), gpts2(i, 3), 'r*');


dem = DEM(peaks(25));
dem.Z(1:15, 1:10) = NaN;
o = [7 10 3];
d = [1 0 -0.25];
X = dem.sample(o, d)
figure
dem.plot(3); hold on
plot3(o(:, 1), o(:, 2), o(:, 3), 'k*');
quiver3(o(:, 1), o(:, 2), o(:, 3), d(:, 1), d(:, 2), d(:, 3), 10, 'g');
plot3(X(:, 1), X(:, 2), X(:, 3), 'r*');
% plane
P0 = [16 16 5];
N = [-0.5 -0.5 1];
plane = createPlane(P0, N);
Xp = intersectRayPlane(o, d, [N, -sum(P0 .* N)]);
drawPlane3d(plane)
plot3(Xp(:, 1), Xp(:, 2), Xp(:, 3), 'y*');
