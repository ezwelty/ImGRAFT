function ZEval = trival(TRI, XEval, YEval, interpMethod, extrap)

%TRIVAL Evaluates a triangulated surface created by <a href="matlab:help grid2tri">grid2tri</a> at specified
%points
%
% Z = TRIVAL(TRI,X,Y) evaluates the triangulated surface TRI created by
% <a href="matlab:help grid2tri">grid2tri</a> in the points given by X and Y. X and Y can be scalars, vectors or
% matrices as the ones returned by <a href="matlab:help meshgrid">meshgrid</a> or <a href="matlab:help ndgrid">ndgrid</a>.
%
% Z = TRIVAL(TRI,X,Y,METHOD) optionally specifies the interpolation method.
% Available methods are:
%   'nearest'  - nearest neighbor interpolation
%   'linear'   - linear interpolation (default)
%
% Z = TRIVAL(TRI,X,Y,METHOD,EXTRAP) optionally specifies a method or a value
% use for extrapolation, default is NaN. Available are:
%   'nearest'       - nearest neighbor extrapolation
%   'linear'        - linear extrapolation (only within convex hull)
%   'linearnearest' - linear extrapolation within convex hull and nearest
%                     neighbor outside.
%   scalar          - value used for all extrapolated points
%
%  VERSION 1.0
%  Author: Patrik Forssén 2015
%
%   See also <a href="matlab:help grid2tri">grid2tri</a>, ndgrid, meshgrid.
%


% Default output
ZEval = [];

% Default input
if (nargin < 5)
  extrap = NaN;
end
if (nargin < 4)
  interpMethod = 'linear';
end
if (isempty(interpMethod)), interpMethod = 'linear'; end
if (isempty(extrap)), extrap = NaN; end

% Check input
if (nargin < 3)
  error('Requires at least three input parameters!')
end
if (~isstruct(TRI))
  error('The first input TRI must be a structure as created by ''grid2tri''!')
end
if (~isnumeric(XEval) || ~isreal(XEval))
  error('Second input X must be a real scalar, vector or matrix!')
end
if (~isnumeric(YEval) || ~isreal(YEval))
  error('Third input Y must be a real scalar, vector or matrix!')
end
if (~ischar(interpMethod) || (~strcmpi(interpMethod, 'linear') && ...
    ~strcmpi(interpMethod, 'nearest')))
  error('Fourth input interpMethod must be a string ''linear'' or ''nearest''!')
end
if ((~ischar(extrap) && ~isnumeric(extrap)) || ...
    (ischar(extrap) && ~strcmpi(extrap, 'nearest') && ...
    ~strcmpi(extrap, 'linear') && ~strcmpi(extrap, 'linearnearest')) || ...
    (isnumeric(extrap) && ~isscalar(extrap)))
  error(['Fifth input extrapMethod must be a string ''linear'', ',...
    '''nearest'', ''linearnearest'' or a scalar!'])
end

% Convert to column vectors
if (isvector(XEval))
  XEval = XEval(:);
end
if (isvector(YEval))
  YEval = YEval(:);
end
% Check sizes
if (~isequal(size(XEval), size(YEval)))
  error('Size mismatch between inputs X and Y!')
end

% Rip out fields from TRI structure
try
  x2row     = TRI.interp.x2row;
  y2col     = TRI.interp.y2col;
  ZVec      = TRI.Z;
  DT        = TRI.DT;
  constFlag = 0;
  if (~isempty(DT.Constraints))
    constFlag = 1;
  end
catch
  error('Input TRI is not a valid structure as created by ''grid2tri''!')
end

% Convert to row/column coordinates
rowVec = x2row(XEval(:));
colVec = y2col(YEval(:));

% Only consider interior triangles for interpolation
if (constFlag)
  tf = find(isInterior(TRI.DT));
end

gridPts = [rowVec, colVec];
[triInd, bCenter] = pointLocation(DT, gridPts);
switch interpMethod
  case 'nearest'
    % Nearest neighbor interpolation
    ZEval = NaN*ones(length(rowVec), 1);
    % Have to remove constraints
    DT.Constraints = [];
    ptInd = nearestNeighbor(DT, gridPts);
    if (constFlag && ~strcmpi(extrap, 'nearest'))
      % Remove triangles in undefined areas in the original grid if
      % nearest neighbor extrapolation is not required
      tmpInd        = ismember(triInd, tf);
      ZEval(tmpInd) = ZVec(ptInd(tmpInd));
    else
      % Keep all triangles
      ZEval = ZVec(ptInd);
    end
    
  case 'linear'
    % Linear interpolation
    triVals = NaN*ones(length(rowVec), 3);
    if (constFlag && ~strcmpi(extrap, 'linearnearest') && ...
        ~strcmpi(extrap, 'linear'))
      % Remove triangles in undefined areas in the original grid if
      % linear extrapolation is not required
      triInd(~ismember(triInd, tf)) = NaN;
    end
    triVals(~isnan(triInd), :) = ZVec(DT(triInd(~isnan(triInd)), :));
    ZEval = dot(bCenter', triVals')';
end

% Possibly do extrapolation
if (~strcmpi(extrap, interpMethod))
  % Extrapolation method is different from interpolation method
  switch extrap
    case 'nearest'
      % Linear interpolation used, do nearest neighbor extrapolation
      extrapInd = isnan(ZEval);
      tmpGridPts   = gridPts(extrapInd, :);
      % Have to remove constraints
      DT.Constraints = [];
      ptInd = nearestNeighbor(DT, tmpGridPts);
      ZEval(extrapInd) = ZVec(ptInd);

    case {'linear', 'linearnearest'}
      % Nearest neighbor interpolation used
      if (~constFlag)
        % Do nothing, linear extrapolation cannot be performed
      else
        % Do linear extrapolation
        extrapInd  = isnan(ZEval);
        tmpGridPts = gridPts(extrapInd, :);
        triVals    = NaN*ones(size(tmpGridPts, 1), 3);
        [triInd, bCenter] = pointLocation(DT, tmpGridPts);
        triVals(~isnan(triInd), :) = ZVec(DT(triInd(~isnan(triInd)), :));
        ZEval(extrapInd) = dot(bCenter', triVals')';
      end
      
    otherwise
      % A scalar used in all extrapolated points
      extrapInd = isnan(ZEval);
      ZEval(extrapInd) = extrap;
  end
end
if (strcmpi(extrap, 'linearnearest'))
  % Possibly do nearest neighbor extrapolation outside convex hull
  extrapInd = isnan(ZEval);
  tmpGridPts   = gridPts(extrapInd, :);
  % Have to remove constraints
  DT.Constraints = [];
  ptInd = nearestNeighbor(DT, tmpGridPts);
  ZEval(extrapInd) = ZVec(ptInd);
end
    
% Possibly reshape output
if (~isvector(XEval))
  ZEval = reshape(ZEval, size(XEval, 1), size(XEval, 2));
end