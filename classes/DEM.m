classdef DEM < handle
  % DEM Digital elevation model.
  %
  % DEM Properties:
  % Z     - Grid of values on a regular xy grid
  % xlim  - Outer bounds of the grid in x [left, right]
  % ylim  - Outer bounds of the grid in y [top, bottom]
  %
  % DEM Properties (read-only):
  % nx,ny - Dimensions of grid [size(Z, 2)], [size(Z, 1)]
  % dx,dy - Grid cell size [range(x) / nx], [range(y) / ny]
  %
  % DEM Methods:
  % DEM    - Construct a new DEM object
  % crop   - Return a cropped DEM
  % resize - Return a resized DEM
  % plot   - Plot a DEM as a hillshade in 2D or 3D
  %
  % DEM Methods (cached):
  % zlim   - Outer bounds of Z [min(Z), max(Z)]
  % x,y    - Cell center coordinates as row vectors [left to right], [top to bottom]
  % X,Y    - Cell center coordinates as grids, equivalent to meshgrid(x, y)
  % clear  - Clear object cache
  
  properties (SetAccess = private)
    Z
    xlim
    ylim
  end
  
  properties (Dependent, SetAccess = private)
    nx, ny
    dx, dy
  end

  properties (SetAccess = private, Hidden = true)
    cached_zlim
    cached_x, cached_y
    cached_X, cached_Y
    cached_tin
  end

  methods

    function dem = DEM(varargin)
      % DEM  Construct a new DEM object.
      %
      %   dem = DEM(Z, xlim = [0 size(Z, 2)], ylim = [0 size(Z, 1)])
      %   dem = DEM(Z, x, y)
      %   dem = DEM(Z, X, Y)
      %
      % Inputs: 
      %   Z         - Grid of elevations
      %   xlim,ylim - Outer x,y limis of the grid [left, right], [top, bottom]
      %   x,y       - Cell center coordinate vectors [left to right], [top to bottom]
      %               Monotonic and evenly spaced
      %   X,Y       - Cell center coordinate grids
      %               Equivalent to meshgrid(x, y)
      %
      % Assumes that Z and x,y coordinates are supplied such that:
      % left   @  xlim(1)  , x(1)  , or X(1, 1)
      % right  @  xlim(end), x(end), or X(1, end)
      % top    @  ylim(1)  , y(1)  , or Y(1, 1)
      % bottom @  ylim(end), y(end), or Y(end, 1)
      %
      % See also meshgrid

      if nargin < 1
        return
      end
      dem.Z = varargin{1};
      if nargin < 2
        dem.xlim = [0 size(dem.Z, 2)];
      else
        dem.xlim = varargin{2};
      end
      if nargin < 3
        dem.ylim = [size(dem.Z, 1) 0];
      else
        dem.ylim = varargin{3};
      end
    end

    function set.Z(dem, value)
      if isempty(value) || ~isnumeric(value)
        error('Value must be non-empty and numeric.')
      end
      if size(value, 3) > 1
        warning('Dimensions greater than 2 currently not supported. Using first 2D layer.')
        value = value(:, :, 1);
      end
      dem.Z = value;
    end
    
    function set.xlim(dem, value)
      dem.xlim = dem.parse_xlim(value);
    end
    
    function xlim = parse_xlim(dem, value)
      if isempty(value) || ~isnumeric(value)
        error('Value must be non-empty and numeric.');
      end
      x = [];
      X = [];
      isgrid = all(size(value) > 1);
      isvector = ~isgrid && any(size(value) > 2);
      if isgrid
        if any(value(1, 1) ~= value(:, 1))
          error('X does not have all equal rows.');
        end
        if nargout > 2
          X = value;
        end
        value = value(1, :);
      end
      if size(value, 1) > 1
        value = reshape(value, 1, length(value));
      end
      if isgrid || isvector
        d = diff(value);
        if any(d(1) ~= d)
          error('x is not equally spaced monotonic.');
        end
        if nargout > 1
          x = value;
        end
        xlim = [value(1) - diff(value(1:2)) / 2, value(end) + diff(value(end-1:end)) / 2];
      else
        xlim = value;
      end
      dem.cached_x = x;
      dem.cached_X = X;
    end

    function set.ylim(dem, value)
      dem.ylim = dem.parse_ylim(value);
    end
    
    function ylim = parse_ylim(dem, value)
      if isempty(value) || ~isnumeric(value)
        error('Value must be non-empty and numeric.');
      end
      y = []; Y = [];
      isgrid = all(size(value) > 1);
      isvector = ~isgrid && any(size(value) > 2);
      if isgrid
        if any(value(1, 1) ~= value(1, :))
          error('Y does not have all equal rows.');
        end
        if nargout > 2
          Y = value;
        end
        value = value(:, 1);
      end
      if size(value, 1) > 1
        value = reshape(value, 1, length(value));
      end
      if isgrid || isvector
        d = diff(value);
        if any(d(1) ~= d)
          error('y is not equally spaced monotonic.');
        end
        if nargout > 1
          y = value;
        end
        ylim = [value(1) - diff(value(1:2)) / 2, value(end) + diff(value(end-1:end)) / 2];
      else
        ylim = value;
      end
      dem.cached_y = y;
      dem.cached_Y = Y;
    end
    
    function value = get.nx(dem)
     value = size(dem.Z, 2);
    end
    
    function value = get.ny(dem)
     value = size(dem.Z, 1);
    end
    
    function value = get.dx(dem)
      value = abs(diff(dem.xlim)) / dem.nx;
    end
    
    function value = get.dy(dem)
      value = abs(diff(dem.ylim)) / dem.ny;
    end
    
    function value = zlim(dem)
      if isempty(dem.cached_zlim)
        value = [min(min(dem.Z)), max(max(dem.Z))];
        dem.cached_zlim = value;
      else
        value = dem.cached_zlim;
      end
    end
    
    function value = x(dem)
      if isempty(dem.cached_x)
        value = (min(dem.xlim) + dem.dx / 2):dem.dx:(max(dem.xlim) - dem.dx / 2);
        if dem.xlim(1) > dem.xlim(2)
          value = fliplr(value);
        end
        dem.cached_x = value;
      else
        value = dem.cached_x;
      end
    end
    
    function value = y(dem)
      if isempty(dem.cached_y)
        value = (min(dem.ylim) + dem.dy / 2):dem.dy:(max(dem.ylim) - dem.dy / 2);
        if dem.ylim(1) > dem.ylim(2)
          value = fliplr(value);
        end
        dem.cached_y = value;
      else
        value = dem.cached_y;
      end
    end
    
    function value = X(dem)
      if isempty(dem.cached_X)
        value = repmat(dem.x, size(dem.Z, 1), 1);
        dem.cached_X = value;
      else
        value = dem.cached_X;
      end
    end
    
    function value = Y(dem)
      if isempty(dem.cached_Y)
        value = repmat(dem.y', 1, size(dem.Z, 2));
        dem.cached_Y = value;
      else
        value = dem.cached_Y;
      end
    end
    
    function value = tin(dem)
      if isempty(dem.cached_tin)
        % --- Compute cell z min & max z ---
        zming = ordfilt2(dem.Z, 1, ones(2, 2));
        zmaxg = ordfilt2(dem.Z, 4, ones(2, 2));
        % Strip padded rows and columns
        zming = zming(1:(end - 1), 1:(end - 1));
        zmaxg = zmaxg(1:(end - 1), 1:(end - 1));
        % Trim boundaries to outer center coordinates
        tin_xlim = dem.x([1, end]);
        tin_ylim = dem.y([1, end]);
        % Build offset DEMs
        zmin = DEM(zming, tin_xlim, tin_ylim);
        zmax = DEM(zmaxg, tin_xlim, tin_ylim);
        % --- Build offset DEMs ---
        tin = struct();
        tin.zmin = zmin;
        tin.zmax = zmax;
        value = tin;
        dem.cached_tin = value;
      else
        value = dem.cached_tin;
      end
    end

    function dem = crop(dem, xlim, ylim, zlim)
      % CROP Return a cropped DEM.
      %
      %   dem = dem.crop(xlim = dem.xlim, ylim = dem.ylim, zlim = dem.zlim)
      %
      % Crops a DEM to the specified x, y, and z limits. Includes cells
      % intersected by the boundary (rather than resampling the DEM to conform
      % to the crop boundaries). Values of Z outside zlim are set to NaN.
      %
      % Inputs:
      %   xlim - x crop boundaries
      %   ylim - y crop boundaries
      %   zlim - z crop boundaries

      % --- Crop xlim ---
      if nargin < 2 || isempty(xlim)
        xlim = dem.xlim;
      else
        xlim = range_intersection(sort(dem.xlim), sort(xlim));
        if range(xlim) == 0 || any(size(xlim) == 0)
          error('Crop bounds do not intersect DEM.')
        end
        if sign(diff(xlim)) ~= sign(diff(dem.xlim))
          xlim = flip(xlim);
        end
      end 
      % --- Crop ylim ---
      if nargin < 3 || isempty(ylim)
        ylim = dem.ylim;
      else
        ylim = range_intersection(sort(dem.ylim), sort(ylim));
        if range(ylim) == 0 || any(size(ylim) == 0)
          error('Crop bounds do not intersect DEM.')
        end
        if sign(diff(ylim)) ~= sign(diff(dem.ylim))
          ylim = flip(ylim);
        end
      end
      % --- Convert limits to grid indices ---
      mincol = floor(abs(xlim(1) - dem.xlim(1)) / dem.dx) + 1;
      maxcol = ceil(abs(xlim(2) - dem.xlim(1)) / dem.dx);
      minrow = floor(abs(ylim(1) - dem.ylim(1)) / dem.dy) + 1;
      maxrow = ceil(abs(ylim(2) - dem.ylim(1)) / dem.dy);
      % --- Crop xy ---
      z = dem.Z(minrow:maxrow, mincol:maxcol);
      xlim = interp1([0 dem.nx], dem.xlim, [mincol - 1, maxcol]);
      ylim = interp1([0 dem.ny], dem.ylim, [minrow - 1, maxrow]);
      % --- Crop z ---
      if nargin > 3 && ~isempty(zlim)
        z(z > max(zlim) | z < min(zlim)) = NaN;
      end
      % --- Return cropped DEM ---
      dem = DEM(z, xlim, ylim);
    end

    function dem = resize(dem, scale, method)
      % RESIZE Return a resized DEM.
      %
      %   dem = dem.resize(scale, method = 'bicubic')
      %   dem = dem.resize(size, method = 'bicubic')
      %
      % Resamples a DEM grid by the specified scale factor and interpolation
      % method. Cell size is changed, but xy boundaries remain the same.
      %
      % Inputs:
      %   scale  - Scale factor
      %   size   - Target grid size [nx, ny]
      %   method - Interpolation method, passed to imresize
      %
      % See also imresize
      
      % --- Check inputs ---
      if nargin < 3
        method = 'bicubic';
      end
      % --- Resize grid ---
      z = imresize(dem.Z, flip(scale), method);
      % --- Return resized DEM ---
      dem = DEM(z, dem.xlim, dem.ylim);
    end

    function h = plot(dem, dim, shade)
      % PLOT Plot a DEM.
      %
      %   h = dem.plot(dim = 2, shade - hillshade(dem.Z, dem.x, dem.y))
      %
      % Plots a DEM in either 2 or 3 dimensions with a custom overlay.
      %
      % Inputs:
      %   dim   - Dimension of plot (2 or 3)
      %   shade - Custom shading image, same size as dem.Z
      %
      % Outputs:
      %   h - Figure handle
      % 
      % See also imagesc, surf, hillshade

      % --- Choose dimension ---
      if nargin < 2
        dim = 2;
      end
      % --- Resize large DEM ---
      if (dem.nx * dem.ny > 2e6)
        scale = sqrt(2e6 / (dem.nx * dem.ny));
        dem = dem.resize(scale);
        if nargin > 2
          shade = imresize(shade, scale);
        end
        warning('DEM automatically downsized for faster plotting');
      end
      % --- Plot DEM ---
      if nargin < 3
        shade = hillshade(dem.Z, dem.x, dem.y);
      end
      if dim == 2
        h = imagesc(dem.x, dem.y, shade);
        axis image
        set(gca, 'ydir', 'normal')
      elseif dim == 3
        h = surf(dem.x, dem.y', dem.Z, double(shade), 'EdgeColor', 'none');
        shading interp
        axis equal
      end
      colormap gray
    end

    function [xi, yi] = xy2ind(dem, xy)
      xi = ceil((xy(:, 1) - dem.xlim(1)) ./ diff(dem.xlim) * dem.nx);
      xi(xy(:, 1) == dem.xlim(1)) = 1;
      yi = ceil((xy(:, 2) - dem.ylim(1)) ./ diff(dem.ylim) * dem.ny);
      yi(xy(:, 2) == dem.ylim(1)) = 1;
      nans = xi < 1 | xi > dem.nx | yi < 1 | yi > dem.ny;
      xi(nans) = NaN;
      yi(nans) = NaN;
      if nargout == 1
        xi = sub2ind([dem.ny dem.nx], yi, xi);
      end
    end

    function xyz = ind2xyz(dem, xi, yi)
      if nargin == 2
        ind = xi;
        [yi, xi] = ind2sub([dem.ny dem.nx], ind);
      else
        ind = sub2ind([dem.ny dem.nx], yi, xi);
      end
      xyz = [dem.x(xi)' dem.y(yi)' dem.Z(ind)];
    end

    function [X, edge] = horizon(dem, xyz, angles)
      if nargin < 3
        angles = [0:1:359]';
      end
      n_angles = length(angles);
      dx = cosd(angles); dy = sind(angles);
      rays = [repmat(xyz, n_angles, 1) dx dy repmat(0, n_angles, 1)];
      X = nan(n_angles, 3);
      in = nan(n_angles, 1);
      for i = 1:n_angles
        ray = rays(i, :);
        cells = traverseRayDEM(ray, dem);
        % Convert to upper-left matrix indices (flip y)
        xi = cells(:, 1); yi = dem.ny - (cells(:, 2) - 1);
        % Retrieve true x,y,z based on cell xy
        x = dem.ind2xyz(sub2ind([dem.ny dem.nx], yi, xi));
        elevation = atand((x(:, 3) - xyz(3)) ./ sqrt((x(:, 1) - xyz(1)).^2 + (x(:, 2) - xyz(2)).^2));
        [~, i_max] = nanmax(elevation);
        if nargout > 1
          if i_max == length(elevation) || all(isnan(elevation((i_max + 1):end)))
            edge(i) = true;
          else
            edge(i) = false;
          end
        end
        X(i, :) = x(i_max, :);
      end
    end

    function vis = visible(dem, xyz)
      X = dem.X(:) - xyz(1);
      Y = dem.Y(:) - xyz(2);
      Z = dem.Z(:) - xyz(3);
      d = sqrt(X.^2 + Y.^2 + Z.^2);
      x = (atan2(Y, X) + pi) / (pi * 2);
      y = Z ./ d;

      [~, ix] = sortrows([round(sqrt((X / dem.dx).^2 + (Y / dem.dy).^2)) x]); %round

      loopix = find(diff(x(ix)) < 0);
      vis = true(size(X, 1), 1);

      maxd = max(d); % TODO: optimize
      N = ceil(2 * pi / (dem.dx / maxd)); % number of points in voxel horizon

      voxx = (0:N)' / N;
      voxy = zeros(size(voxx)) - Inf;

      for k = 1:length(loopix) - 1
          lp = ix((loopix(k) + 1):loopix(k + 1));
          lp = lp([end 1:end 1]);
          yy = y(lp); xx = x(lp);
          xx(1) = xx(1) - 1; xx(end) = xx(end) + 1; % TODO: why?
          vis(lp(2:end - 1)) = interp1q(voxx, voxy, xx(2:end - 1)) < yy(2:end - 1);
          voxy = max(voxy, interp1q(xx, yy, voxx));
      end
      vis = reshape(vis, [dem.ny dem.nx]);
    end

    function X = sample(dem, origin, direction, first, xy0, r0)

      tin = dem.tin;

      if nargin < 4 || isempty(first)
        first = true;
      end

      % TODO: filter out points outside DEM
      if nargin < 3
        % sample at origin
        [xi, yi] = tin.Dzmin.xy2ind(origin);
        dx = origin(1) - (tin.Dzmin.xlim(1) + (xi - 1) * sign(diff(tin.Dzmin.xlim)) * dem.dx);
        dy = origin(2) - (tin.Dzmin.ylim(1) + yi * sign(diff(tin.Dzmin.ylim)) * dem.dy);
        if dy >= dx * (tin.Dzmin.dy / tin.Dzmin.dx)
          % upper triangle
          ind = sub2ind(size(dem.Z), [yi + 1; yi; yi], [xi; xi; xi + 1]);
        else
          % lower triangle
          ind = sub2ind(size(dem.Z), [yi; yi + 1; yi + 1], [xi + 1; xi + 1; xi]);
        end
        tri = [dem.X(ind), dem.Y(ind), dem.Z(ind)];
        u = tri(1, :) - tri(3, :);
        v = tri(2, :) - tri(3, :);
        n = cross(u, v);
        d = n(1) * tri(1, 1) + n(2) * tri(1, 2) + n(3) * tri(1, 3);
        z = (d - (n(1) * origin(1) + n(2) * origin(2))) / n(3);
        X = [origin(1), origin(2), z];
        return
      end

      % TODO: compare to intersectLineMesh3d in mesh3d | intersectLineTriangle3d, intersectRayPolygon3d in gemo3d
      % TODO: vectorize for shared origin (as when emanating from camera)
      % TODO: Pre-compute TIN and cache in DEM object.
      % TODO: Pre-compute all triangles?
      % TODO: Split DEM into multiple boxes and test each box in turn? (quadtree?)

      % Test for ray - box intersection
      if nargin < 6
        [tmin, tmax] = intersectRayBox([origin direction], [tin.min tin.max]);
        imin = [1 1];
        imax = [tin.nx tin.ny];
      else
        % Apply search radius (if specified)
        % FIXME: Broken
        % Snap search radius to grid cell boundaries
        imin = floor((xy0(1:2) - r0 - tin.min(1:2)) ./ [tin.dx, tin.dy]);
        imax = ceil((xy0(1:2) + r0 - tin.min(1:2)) ./ [tin.dx, tin.dy]);
        smin = tin.min(1:2) + imin .* [tin.dx, tin.dy];
        smax = tin.min(1:2) + imax .* [tin.dx, tin.dy];
        % Intersect with outer grid boundaries
        smin = max([tin.min(1:2) ; smin]);
        smax = min([tin.max(1:2) ; smax]);
        [tmin, tmax] = intersectRayBox([origin direction], [smin tin.min(3) smax tin.max(3)]);
        % imin = round((smin - tin.min(1:2)) ./ [tin.dx, tin.dy]);
        % imax = round((smax - tin.min(1:2)) ./ [tin.dx, tin.dy]);
      end

      % Traverse grid (2D)
      X = [];
      if ~isempty(tmin)

        % Compute endpoints of ray within grid
        start = origin + tmin * direction;
        stop = origin + tmax * direction;

        % Find starting voxel coordinates
        x = ceil((start(1) - tin.min(1)) / tin.dx);
        y = ceil((start(2) - tin.min(2)) / tin.dy);
        % Snap to 1 (from 0)
        if x < 1
          x = 1;
        end
        if y < 1
          y = 1;
        end
        % Snap to nx, ny (from above)
        % (Necessary for rounding errors)
        if x > tin.nx
          x = tin.nx;
        end
        if y > tin.ny
          y = tin.ny;
        end

        % Set x,y increments based on ray slope
        if direction(1) >= 0
          % Increasing x
          tCellX = x / tin.nx;
          stepX = 1;
        else
          % Decreasing x
          tCellX = (x - 1) / tin.nx;
          stepX = -1;
        end
        if direction(2) >= 0
          % Increasing y
          tCellY = y / tin.ny;
          stepY = 1;
        else
          % Decreasing y
          tCellY = (y-1) / tin.ny;
          stepY = -1;
        end

        % TODO: ?
        boxSize = tin.max - tin.min;
        cellMaxX = tin.min(1) + tCellX * boxSize(1);
        cellMaxY = tin.min(2) + tCellY * boxSize(2);

        % Compute values of t at which ray crosses vertical, horizontal voxel boundaries
        tMaxX = tmin + (cellMaxX - start(1)) / direction(1);
        tMaxY = tmin + (cellMaxY - start(2)) / direction(2);

        % Width and height of voxel in t
        tDeltaX = tin.dx / abs(direction(1));
        tDeltaY = tin.dy / abs(direction(2));

        % Find ending voxel coordinates
        mx = ceil((stop(1) - tin.min(1)) / tin.dx);
        my = ceil((stop(2) - tin.min(2)) / tin.dy);
        % TODO: Snap to grid boundaries?
        % (e.g. 0 -> 1)

        % Return list of traversed voxels
        %voxels = [];
        z_in = start(3);
        i_cell = 1;
        while (x <= tin.nx) && (x >= 1) && (y <= tin.ny) && (y >= 1)

          % TODO: Check if origin is below surface

          % Check if within current voxel (in Z)
          % either entered voxel at tMaxX or tMaxY
          % if tMaxX: x != voxels(end, 1)
          % if tMaxY: y != voxels(end, 2)
          z_out = origin(3) + min(tMaxY, tMaxX) * direction(3);

          % Test for intersection of both possible triangles
          % Convert to upper-left matrix indices (flip y)
          % TODO: avoidable?
          yi = tin.ny - y + 1;
          if ~(isnan(tin.zmax(yi, x)) || isnan(tin.zmax(yi, x))) && ~(z_in > tin.zmax(yi, x) && z_out > tin.zmax(yi, x) || z_in < tin.zmin(yi, x) && z_out < tin.zmin(yi, x))
            sqi = sub2ind(size(dem.Z), [yi + 1; yi; yi; yi + 1], [x; x; x + 1; x + 1]);
            square = [dem.X(sqi), dem.Y(sqi), dem.Z(sqi)];
            tri1 = square([1 2 3], :);
            tri2 = square([3 4 1], :);
            [intersection, ~, ~, t] = rayTriangleIntersection(origin, direction, tri1(1, :), tri1(2, :), tri1(3, :));
            if ~intersection
              [intersection, ~, ~, t] = rayTriangleIntersection(origin, direction, tri2(1, :), tri2(2, :), tri2(3, :));
            end
            if intersection
              X(end + 1, :) = origin + t * direction;
              if first
                break
              end
            end
          end

          if x == mx && y == my
            % X = NaN;
            % voxels(end + 1, :) = [x y];
            break
          end

          % return voxels
          % voxels(end + 1, :) = [x y];

          if tMaxX < tMaxY
            tMaxX = tMaxX + tDeltaX;
            x = x + stepX;
          else
            tMaxY = tMaxY + tDeltaY;
            y = y + stepY;
          end
          z_in = z_out;
        end
      end
      if isempty(X)
        X = NaN;
      end
    end

  end % methods

  methods (Static)

  end % static methods

  methods (Access = private)

  end % private methods

end % classdef
