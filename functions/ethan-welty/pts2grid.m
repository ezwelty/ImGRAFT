function [xvec, yvec, zgrid] = pts2grid(x, y, z, dx, dy)

% Process point set
[xmin, xi] = min(x);
[ymax, yi] = max(y);

xlim = [xmin - dx / 2, max(x) + dx / 2];
ylim = [min(y) - dy / 2, ymax + dy / 2];

nx = ceil(range(xlim) / dx);
ny = ceil(range(ylim) / dy);
xlim = [xlim(1), xlim(1) + nx * dx];
ylim = [ylim(1), ylim(1) + ny * dy];

subx = ceil((x - xlim(1)) / dx); subx(xi) = 1;
suby = ceil((ylim(2) - y) / dy); suby(yi) = 1;

% Calculate cell statistic
ind = sub2ind([ny nx], suby, subx);
[sorted, si] = sort(ind);
sz = z(si);
[u, I, J] = unique(sorted, 'rows', 'first');
I = [I ; size(subx, 1) + 1];
N = length(u);
uz = zeros(N, 1);
for i = 1:N
  uz(i) = mean(sz(I(i):(I(i + 1) - 1)));
end

% Grid points
zgrid = nan(ny, nx);
zgrid(u) = uz;
xvec = min(xlim):dx:(max(xlim));
yvec = max(ylim):-dy:min(ylim)';
