% TODO: Compute x, y, X, Y only as needed, and expire instead of automatically updating.

classdef DEM < handle
  % DEM Digital elevation model
  %
  % This model implements a regularly gridded digital elevation model. Read-only
  % properties are cached for efficiency and updated only when needed. Hidden
  % properties are cached and only computed when requested.
  %
  % DEM Properties:
  % Z    - Grid of values on a regular xy grid
  % xlim - Outer bounds of the grid in x [left, right]
  % ylim - Outer bounds of the grid in y [top, bottom]
  %
  % DEM Properties (read-only):
  % min   - Minimum corner of bounding box [min(x), min(y), min(z)]
  % max   - Maximum corner of bounding box [max(x), max(y), max(z)]
  % nx,ny - Dimensions of grid
  % dx,dy - Cell size in x and y
  %
  % DEM Properties (read-only and hidden):
  % x,y   - Cell center coordinates as row vectors [left to right], [top to bottom]
  % X,Y   - Cell center coordinates as grids, equivalent to meshgrid(x, y)
  %
  % DEM Methods:
  % DEM    - Construct a new DEM object
  % crop   - Crop a DEM to new x, y, and z boundaries
  % resize - Resize a DEM by a given scale factor and interpolaton method
  % plot   - Plot a DEM as a hillshade in 2D or 3D

  properties
    xlim = [0 3];
    ylim = [3 0];
    Z = zeros(3);
  end

  properties (SetAccess = private)
    nx, ny
    dx, dy
    min, max
  end

  properties (SetAccess = private, Hidden = true, GetObservable)
    x, y
    X, Y
  end

  methods

    function dem = DEM(varargin)
      % DEM  Construct a new DEM object.
      %
      %   dem = DEM(Z, xlim = [0 size(Z, 2)], ylim = [0 size(Z, 1)])
      %   dem = DEM(Z, x, y)
      %   dem = DEM(Z, X, Y)
      %
      % Constructs a DEM object from a grid of elevations and either the outer
      % x,y limits of the grid (xlim, ylim) or vectors (x, y) or grids (X, Y) of
      % cell center coordinates. x,y must be monotonic and evenly spaced, and
      % X,Y must be equivalent to meshgrid(x, y).
      %
      % Assumes that Z and x,y outer boundaries are supplied such that:
      % left   @  xlim(1)  , x(1)  , or X(1, 1)
      % right  @  xlim(end), x(end), or X(1, end)
      % top    @  ylim(1)  , y(1)  , or Y(1, 1)
      % bottom @  ylim(end), y(end), or Y(end, 1)
      %
      % See also: meshgrid, DEM.parseZ, DEM.parseXlim, DEM.parseYlim

      % ()
      if nargin == 0
        dem = dem.updateXYZ;
      end
      % (Z)
      if nargin == 1
        dem.Z = DEM.parseZ(varargin{1});
        dem.xlim = [0 size(dem.Z, 2)];
        dem.ylim = [0 size(dem.Z, 1)];
        dem = dem.updateXYZ;
      end
      % (Z, xlim, ylim)
      % (Z, x, y)
      % (Z, X, Y)
      if nargin == 3
        dem.Z = DEM.parseZ(varargin{1});
        [dem.xlim, x, X] = DEM.parseXlim(varargin{2});
        [dem.ylim, y, Y] = DEM.parseYlim(varargin{3});
        dem = dem.updateX(x, X).updateY(y, Y).updateZ;
      end

      % Add listeners
      addlistener(dem, 'x', 'PreGet', @dem.PreGetProperty);
      addlistener(dem, 'y', 'PreGet', @dem.PreGetProperty);
      addlistener(dem, 'X', 'PreGet', @dem.PreGetProperty);
      addlistener(dem, 'Y', 'PreGet', @dem.PreGetProperty);
    end

    function PreGetProperty(dem, property, varargin)
      if isempty(dem.(property.Name))
        dem = dem.(['compute_' property.Name]);
      end
    end

    % Independent properties

    function dem = set.xlim(dem, value)
      value = DEM.parseXlim(value);
      if any(value ~= dem.xlim)
        dem.xlim = value;
        dem = dem.updateX;
      end
    end

    function dem = set.ylim(dem, value)
      value = DEM.parseYlim(value);
      if any(value ~= dem.ylim)
        dem.ylim = value;
        dem = dem.updateY;
      end
    end

    function dem = set.Z(dem, value)
      value = DEM.parseZ(value);
      dem.Z = value;
      if size(dem.Z, 2) ~= dem.nx, dem = dem.updateX; end
      if size(dem.Z, 1) ~= dem.ny, dem = dem.updateY; end
      dem = dem.updateZ;
    end

    % Dependent properties
    % function value = get.nx(dem)
    %   value = size(dem.Z, 2);
    % end
    % function value = get.ny(dem)
    %   value = size(dem.Z, 1);
    % end
    % function value = get.dx(dem)
    %   value = abs(diff(dem.xlim)) / dem.nx;
    % end
    % function value = get.dy(dem)
    %   value = abs(diff(dem.ylim)) / dem.ny;
    % end
    % function value = get.min(dem)
    %   value = [min(dem.xlim), min(dem.ylim), min(min(dem.Z))];
    % end
    % function value = get.max(dem)
    %   value = [max(dem.xlim), max(dem.ylim), max(max(dem.Z))];
    % end
    % function value = get.x(dem)
    %   value = (dem.min(1) + dem.dx / 2):dem.dx:(dem.max(1) - dem.dx / 2);
    %   if (dem.xlim(1) > dem.xlim(2))
    %     value = flip(value);
    %   end
    % end
    % function value = get.y(dem)
    %   value = (dem.min(2) + dem.dy / 2):dem.dy:(dem.max(2) - dem.dy / 2);
    %   if (dem.ylim(1) > dem.ylim(2))
    %     value = flip(value);
    %   end
    % end
    % function value = get.X(dem)
    %   value = repmat(dem.x, dem.ny, 1);
    % end
    % function value = get.Y(dem)
    %   value = repmat(dem.y', 1, dem.nx);
    % end

    function dem = crop(dem, xlim, ylim, zlim)
      % CROP  Crop a DEM.
      %
      %   dem = dem.crop(xlim, ylim = [], zlim = [])
      %
      % Crops a DEM to the specified x, y, and z limits. Includes cells
      % intersected by the boundary (rather than resampling the DEM to conform
      % to the crop boundaries). Values of Z outside zlim are set to NaN.
      %
      % Inputs:
      %   xlim - x crop boundaries
      %   ylim - y crop boundaries
      %   zlim - z crop boundaries

      % Check xlim
      if isempty(xlim)
        % Use existing xlim
        xlim = dem.xlim;
        mincol = 1;
        maxcol = dem.nx;
      else
        % Parse custom xlim
        xlim = DEM.parseXlim(xlim);
        if sign(diff(xlim)) ~= sign(diff(dem.xlim)), xlim = flip(xlim); end
        % Fit to current limits
        dx1 = sign(xlim - dem.xlim(1)); dx2 = sign(xlim - dem.xlim(2));
        % No overlap: if both points are on same side of both limits
        if isequal(dx1(1), dx1(2), dx2(1), dx2(2))
          error('Crop bounds do not intersect DEM.')
        end
        % Beyond current limits: if either point same side of both limits
        xbeyond = dx1 == dx2;
        xlim(xbeyond) = dem.xlim(xbeyond);
        % Convert limits to column indices
        mincol = floor(abs(xlim(1) - dem.xlim(1)) / dem.dx) + 1;
        maxcol = ceil(abs(xlim(2) - dem.xlim(1)) / dem.dx);
      end

      % Check ylim
      if nargin < 3 || isempty(ylim)
        % Use existing ylim
        ylim = dem.ylim;
        minrow = 1;
        maxrow = dem.ny;
      else
        % Parse custom ylim
        ylim = DEM.parseYlim(ylim);
        if sign(diff(ylim)) ~= sign(diff(dem.ylim)), ylim = flip(ylim); end
        % Fit to current limits
        dy1 = sign(ylim - dem.ylim(1)); dy2 = sign(ylim - dem.ylim(2));
        % No overlap: if both points are on same side of both limits
        if isequal(dy1(1), dy1(2), dy2(1), dy2(2))
          error('Crop bounds do not intersect DEM.')
        end
        % Beyond current limits: if either point same side of both limits
        ybeyond = dy1 == dy2;
        ylim(ybeyond) = dem.ylim(ybeyond);
        % Convert limits to row indices
        minrow = floor(abs(ylim(1) - dem.ylim(1)) / dem.dy) + 1;
        maxrow = ceil(abs(ylim(2) - dem.ylim(1)) / dem.dy);
      end

      % Perform clip
      cZ = dem.Z(minrow:maxrow, mincol:maxcol);
      cxlim = interp1([0 dem.nx], dem.xlim, [mincol - 1 maxcol]);
      cylim = interp1([0 dem.ny], dem.ylim, [minrow - 1 maxrow]);

      % Check zlim
      if nargin > 3 && ~isempty(zlim)
        % Set values outside range to NaN
        cZ(cZ > max(zlim) | cZ < min(zlim)) = NaN;
      end

      % Create clipped DEM
      dem = DEM(cZ, cxlim, cylim);
    end

    function dem = resize(dem, scale, method)
      % RESIZE  Resize a DEM.
      %
      %   dem = dem.resize(scale, method = 'bicubic')
      %
      % Resamples a DEM grid by the specified scale factor and interpolation
      % method. The cell size is changed, but grid xy boundaries remain the same.
      %
      % Inputs:
      %   scale  - Scale factor [s] or [srows, scols], or output size [nrows, ncols]
      %   method - Interpolation method passed to imresize
      %
      % See also: imresize

      if nargin < 3
        method = 'bicubic';
      end

      % Scale grid
      Z = imresize(dem.Z, scale, method);
      dem = DEM(Z, dem.xlim, dem.ylim);
    end

    function h = plot(dem, dim, shade)
      % PLOT  Plot a DEM in 2 or 3 dimensions.
      %
      %   h = dem.plot(dim)
      %
      % Plots a DEM at the specified scale in either 2 or 3 dimensions.
      %
      % Inputs:
      %   dim - Dimension of plot (2 or 3)
      %
      % Outputs:
      %   h   - Figure handle

      % Force scaling of large DEM for 3D plotting
      if (dem.nx * dem.ny > 2e6)
        scale = sqrt(2e6 / (dem.nx * dem.ny));
        dem = dem.resize(scale);
        if (nargin > 2)
          shade = imresize(shade, scale);
        end
        warning('DEM automatically downsized for fast plotting')
      end

      % Plot DEM
      if (nargin < 3)
        shade = hillshade(dem.Z, dem.x, dem.y);
      end
      if (dim == 2)
        h = imagesc(dem.x, dem.y, shade);
        axis image
        set(gca, 'ydir', 'normal')
      elseif (dim == 3)
        h = surf(dem.x, dem.y, dem.Z, double(shade), 'EdgeColor','none');
        shading interp
        axis equal
      end
      colormap gray
    end

    function ind = xy2ind(dem, xy)
      xi = ceil((xy(:, 1) - dem.xlim(1)) ./ diff(dem.xlim) * dem.nx);
      xi(xy(:, 1) == dem.xlim(1)) = 1;
      yi = ceil((xy(:, 2) - dem.ylim(1)) ./ diff(dem.ylim) * dem.ny);
      yi(xy(:, 2) == dem.ylim(1)) = 1;
      nans = xi < 1 | xi > dem.nx | yi < 1 | yi > dem.ny;
      xi(nans) = NaN;
      yi(nans) = NaN;
      ind = sub2ind(size(dem.Z), yi, xi);
    end

  end % methods

  methods (Static)

    function [xlim, x, X] = parseXlim(value)
      if (isempty(value) || ~isnumeric(value)), error('value must be non-empty and numeric'); end
      if nargout > 1, x = []; X = []; end
      isgrid = all(size(value) > 1);
      isvector = ~isgrid && any(size(value) > 2);
      if isgrid
        if any(value(1, 1) ~= value(:, 1)), error('X does not have all equal rows'); end
        if nargout > 2, X = value; end
        value = value(1, :);
      end
      if size(value, 1) > 1
        value = reshape(value, 1, length(value));
      end
      if isgrid || isvector
        dx = diff(value);
        if any(dx(1) ~= dx), error('x is not equally spaced monotonic'); end
        if nargout > 1, x = value; end
        xlim = [value(1) - diff(value(1:2)) / 2, value(end) + diff(value(end-1:end)) / 2];
      else
        xlim = value;
      end
    end

    function [ylim, y, Y] = parseYlim(value)
      if (isempty(value) || ~isnumeric(value)), error('value must be non-empty and numeric'); end
      if nargout > 1, y = []; Y = []; end
      isgrid = all(size(value) > 1);
      isvector = ~isgrid && any(size(value) > 2);
      if isgrid
        if any(value(1, 1) ~= value(1, :)), error('Y does not have all equal rows'); end
        if nargout > 2, Y = value; end
        value = value(:, 1);
      end
      if size(value, 1) > 1
        value = reshape(value, 1, length(value));
      end
      if isgrid || isvector
        dy = diff(value);
        if any(dy(1) ~= dy), error('y is not equally spaced monotonic'); end
        if nargout > 1, y = value; end
        ylim = [value(1) - diff(value(1:2)) / 2, value(end) + diff(value(end-1:end)) / 2];
      else
        ylim = value;
      end
    end

    function value = parseZ(value)
      if (isempty(value) || ~isnumeric(value)), error('value must be non-empty and numeric'); end
      if (size(value, 3) > 1)
        warning('3D+ Z currently not supported. Using first 2D layer.')
        value = value(:, :, 1);
      end
    end

  end % static methods

  methods (Access = private)

    function dem = updateXYZ(dem)
      dem = dem.updateX.updateY.updateZ;
    end

    function dem = updateX(dem, x, X)
      dem.min(1) = min(dem.xlim);
      dem.max(1) = max(dem.xlim);
      dem.nx = size(dem.Z, 2);
      dem.dx = abs(diff(dem.xlim)) / dem.nx;
      if nargin > 2 && ~isempty(x)
        dem.x = x;
      else
        dem.x = [];
      end
      if nargin > 3 && ~isempty(X)
        dem.X = X;
      else
        dem.X = [];
      end
    end

    function dem = compute_x(dem)
      dem.x = (dem.min(1) + dem.dx / 2):dem.dx:(dem.max(1) - dem.dx / 2);
      if dem.xlim(1) > dem.xlim(2)
        dem.x = fliplr(dem.x);
      end
    end

    function dem = compute_X(dem)
      dem.X = repmat(dem.x, size(dem.Z, 1), 1);
    end

    function dem = updateY(dem, y, Y)
      dem.min(2) = min(dem.ylim);
      dem.max(2) = max(dem.ylim);
      dem.ny = size(dem.Z, 1);
      dem.dy = abs(diff(dem.ylim)) / dem.ny;
      if nargin > 2 && ~isempty(y)
        dem.y = y;
      else
        dem.y = [];
      end
      if nargin > 3 && ~isempty(Y)
        dem.Y = Y;
      else
        dem.Y = [];
      end
    end

    function dem = compute_y(dem)
      dem.y = (dem.min(2) + dem.dy / 2):dem.dy:(dem.max(2) - dem.dy / 2);
      if dem.ylim(1) > dem.ylim(2)
        dem.y = fliplr(dem.y);
      end
    end

    function dem = compute_Y(dem)
      dem.Y = repmat(dem.y', 1, size(dem.Z, 2));
    end

    function dem = updateZ(dem)
      dem.min(3) = min(min(dem.Z));
      dem.max(3) = max(max(dem.Z));
    end

  end % private methods

end % classdef
