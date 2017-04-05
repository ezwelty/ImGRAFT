function varargout = grid2tri(varargin)

%GRID2TRI Triangulates a gridded surface with given error tolerances
%
%  Works by iteratively adding vertices to a Delaunay triangulation at grid
%  points until the error between the triangulation and the supplied surface is
%  less than given error tolerances. Undefined areas in the supplied surface are
%  supported, the number boundary points between defined and undefined areas are
%  reduced by the recursive Douglas-Peucker Polyline Simplification as
%  implemented by Wolfgang Schwanghart in <a href="matlab:help dpsimplify">dpsimplify</a>.
%
%  GRID2TRI(Z) triangulates the surface in the matrix Z, NaN-values in
%  Z indicates areas that should not be triangulated, e.g. holes.
%
%  GRID2TRI(X,Y,Z) triangulates Z on the grid given by X and Y.
%  X and Y can be matrices as the ones returned by <a href="matlab:help meshgrid">meshgrid</a> or <a href="matlab:help ndgrid">ndgrid</a>, or
%  they can be vectors containing the x- and y-coordinates for the grid.
%
%  TRI = GRID2TRI(...) returns a structure TRI containing the
%  triangulation. TRI has the following fields,
%    DT       - triangulation object as returned by <a href="matlab:help delaunayTriangulation">delaunayTriangulation</a>
%    Z        - the z-coordinates for the vertices given in DT
%    XGrid    - a vector with x-coordinates for the used grid
%    YGrid    - a vector with y-coordinates for the used grid
%    interp   - substructure used by <a href="matlab:help trival">trival</a> to evaluate the triangulation
%    hotStart - substructure for hot start
%  Note that the vertex-coordinates in DT are given as indices in XGrid and
%  YGrid, not as x- and y-coordinates. You can plot the triangulation DT using
%  <a href="matlab:help triplot">triplot</a>, <a href="matlab:help trimesh">trimesh</a> and <a href="matlab:help trisurf">trisurf</a>. To plot triangulations on grids with undefined 
%  areas use the <a href="matlab:help isInterior ">isInterior</a> method, for example:
%  tf = isInterior(DT); 
%  triplot(DT.ConnectivityList(tf, :), DT.Points(:, 1), DT.Points(:, 2))    
%
%  GRID2TRI(TRI) continues from an existing triangulation in the structure TRI,
%  i.e., does a "hot-start"
%  
%  GRID2TRI(N) does a demo triangulation of MATLAB's <a href="matlab:help peaks">peaks</a> surface, the optional
%  parameter N is the grid size used, default is 49.
%
%  GRID2TRI(...,propertyName,propertyValue) allows additional parameter/value
%  pairs in the call to control the triangulation, plotting etc. See <a href="matlab:help trioptions">trioptions</a>
%  for available settings.
%
%  VERSION 1.0
%  Author: Patrik Forssén 2015
%
%   See also <a href="matlab:help trioptions">trioptions</a>, <a href="matlab:help trival">trival</a>, <a href="matlab:help triplot">triplot</a>, <a href="matlab:help trimesh">trimesh</a>, <a href="matlab:help trisurf">trisurf</a>, <a href="matlab:help isInterior ">isInterior</a>, 
%   <a href="matlab:help dpsimplify">dpsimplify</a>, delaunayTriangulation, ndgrid, meshgrid, peaks.
%

% Parse input
[XGrid, YGrid, ZMat, options] = parseInput(varargin);

% Check output
if (nargout > 1)
  error('Incorrect number of output parameters, should be none or one!')
end

hotStart = 0;
if (isempty(YGrid))
  
  % Hot-start, a TRI structure was supplied
  TRI       = XGrid;
  hotStart  = 1;
  try
    DT     = TRI.DT;
    XGrid  = TRI.XGrid;
    YGrid  = TRI.YGrid;
    ZMat   = TRI.hotStart.ZMat;
    bndMat = TRI.hotStart.bndMat;
    iterStatus.ZMatSize = size(ZMat);
    iterStatus.nIterOffset = TRI.hotStart.nIterOffset;
    iterStatus.tOffset     = TRI.hotStart.tOffset;
    unusedInd              = TRI.hotStart.unusedInd;  
    iterStatus.constFlag   = 0;
    if (~isempty(bndMat))
      iterStatus.constFlag = 1;
    end
    % Check if relative errors can be used
    iterStatus.relErrFlag = 1;
    if (min(ZMat(:)) <= 0 && max(ZMat(:)) >= 0)
      iterStatus.relErrFlag = 0;
    end
  catch
    error('Input is not a valid TRI structure!')
  end   
else
  
  % Cold start
  iterStatus.ZMatSize = size(ZMat);
  % Check if relative errors can be used
  iterStatus.relErrFlag = 1;
  if (min(ZMat(:)) <= 0 && max(ZMat(:)) >= 0)
    iterStatus.relErrFlag = 0;
  end
  
  if (max(isnan(ZMat(:))))
    iterStatus.constFlag = 1;  
    % Need to discretize boundaries
    bndMat = discretizeBoundary(ZMat);
    [constMat, pts0List] = bndMat2const(bndMat, iterStatus.ZMatSize);
  else
    % Add corner points
    iterStatus.constFlag = 0;
    bndMat   = [];
    pts0List = [...
      1            , 1            , ZMat(1  , 1  ); ...
      1            , size(ZMat, 2), ZMat(1  , end); ...
      size(ZMat, 1), 1            , ZMat(end, 1  ); ...
      size(ZMat, 1), size(ZMat, 2), ZMat(end, end)];
  end
  
  % Make an initial triangulation
  if (iterStatus.constFlag)
    DT = delaunayTriangulation(pts0List(:, 1:2), constMat);
  else
    DT = delaunayTriangulation(pts0List(:, 1:2));
  end
  
  % Initial iteration parameters
  iterStatus.nIterOffset = 0;
  iterStatus.tOffset     = 0;
end

% Used to convert row/column indices to x/y coordinates 
iterStatus.row2x = griddedInterpolant(1:length(XGrid), XGrid, ...
  'linear', 'linear');
iterStatus.col2y = griddedInterpolant(1:length(YGrid), YGrid, ...
  'linear', 'linear');

% Use linear indexing
[rowMat, colMat] = ndgrid(1:iterStatus.ZMatSize(1), 1:iterStatus.ZMatSize(2));
rowVec    = rowMat(:);
colVec    = colMat(:);
ZVec      = ZMat(:);
ZTriVec   = ZVec;
if (~hotStart)
  tmp       = sub2ind(iterStatus.ZMatSize, pts0List(:, 1), pts0List(:, 2));
  unusedInd = setdiff(setdiff(1:length(ZVec), find(isnan(ZVec))), tmp)';
end

% Make a convergence plot
if (strcmpi(options.PlotConv, 'on'))
  iterPlot('init', iterStatus, options)
end

% Main iteration
iterStatus.nIter   = 0;
iterStatus.tStart  = tic;
while (1)
    
  % Get current triangulation
  ptsList  = DT.Points;
  triList  = DT.ConnectivityList;
  constMat = DT.Constraints;
  iterStatus.nPts = size(ptsList, 1);
  iterStatus.nTri = size(triList, 1);
  
  % Evaluate triangulation at unused grid points
  if (~isempty(unusedInd))
    gridPts = [rowVec(unusedInd), colVec(unusedInd)];
    ptsInd  = sub2ind(iterStatus.ZMatSize, ptsList(:, 1), ptsList(:, 2));
    [triInd, bCenter] = pointLocation(DT, gridPts);
    triVals  = ZVec(ptsInd(DT(triInd, :)));
    ZTriEval = dot(bCenter', triVals')';
    ZTriVec(unusedInd) = ZTriEval;
  else
    ZTriVec = [];
  end
  
  % Calculate error and check for convergence
  ZErrVec = ZTriVec(unusedInd) - ZVec(unusedInd);
  iterStatus.errors.errVec = ZErrVec;
  iterStatus.errors.absErr = abs(ZErrVec);
  if (iterStatus.relErrFlag)
    iterStatus.errors.relErr = abs(ZErrVec./ZVec(unusedInd));
  end
  iterStatus = chkConvFcn(iterStatus, options);
  if (iterStatus.conv), break, end
        
  % Display 
  if (isnumeric(options.Verbose) && mod(iterStatus.nIter, options.Verbose) == 0)
    displayFcn('iter', iterStatus)
  end
  
  % Plot convergence
  if (strcmpi(options.PlotConv, 'on'))
    iterPlot('iter', iterStatus, options)
  end
  
  % Get vertices to add
  if (strcmpi(options.VertexSelection, 'abs'))
    [~, addInd] = sort(iterStatus.errors.absErr, 'descend');
  else
    [~, addInd] = sort(iterStatus.errors.relErr, 'descend');
  end
  if (options.VerticesAdded < 1)
    nAddPts = round(options.VerticesAdded*iterStatus.nPts);
  else
    nAddPts = options.VerticesAdded;
  end
  nAddPts = max(1, min(nAddPts, length(ZErrVec)));
  if (strcmpi(options.AddSeparateTri, 'on'))
    % The added vertices must be in different triangles
    addTriMat = [triInd(addInd), addInd];
    addInd    = [];
    while (~isempty(addTriMat))
      addInd = [addInd; addTriMat(1, 2)];
      if (length(addInd) == nAddPts)
        break
      end
      addTriMat = addTriMat(addTriMat(:, 1) ~= addTriMat(1, 1), :);
    end
  else
    % The added vertices can be in the same triangle
    addInd  = addInd(1:nAddPts);
  end
  ptInd   = unusedInd(addInd);
  newPts  = [rowVec(ptInd), colVec(ptInd)];
  
  % Add vertices
  if (iterStatus.constFlag)
    % Need to check if any of the new points are boundary points!
    bndInd = ismember(ptInd, bndMat(:, 2));
    bndPts = newPts( bndInd, :);
    intPts = newPts(~bndInd, :);
    if (~isempty(bndPts))
      % Have to add points on the boundary!
      % Current interior points
      oldIntPts = ptsList(size(constMat, 1)+1:end, :);
      % Set boundary points to active in bndMat
      bndMat(ismember(bndMat(:, 2), ptInd), 3) = 1;
      % Update constraints
      [newConstMat, bndPtsList] = bndMat2const(bndMat, iterStatus.ZMatSize);
      newPtsList = [bndPtsList; oldIntPts; intPts];
      % Make new triangulation (cannot be updated if new constraints)
      DT = delaunayTriangulation(newPtsList, newConstMat);
    else
      % Update triangulation
      DT.Points(end+(1:size(intPts, 1)), :) = intPts;
    end
  else
    % Update triangulation
    DT.Points(end+(1:size(newPts, 1)), :) = newPts;
  end
  % Remove added points from the list if unused grid points
  unusedInd = setdiff(unusedInd, ptInd);
   
  % Update number of iterations
  iterStatus.nIter = iterStatus.nIter + 1;
end


if (strcmpi(options.PlotConv, 'on'))
  % Close convergence plot
  iterPlot('final', iterStatus, options)
end

if (iterStatus.constFlag)
  % Get interior triangles
  tf                  = isInterior(DT);
  triList             = DT.ConnectivityList(tf, :);
  iterStatus.nTri     = size(triList, 1);
end

% Set TRI structure
ptsInd    = sub2ind(iterStatus.ZMatSize, ptsList(:, 1), ptsList(:, 2));
TRI.DT    = DT;
TRI.Z     = ZVec(ptsInd);
TRI.XGrid = XGrid;
TRI.YGrid = YGrid;
TRI.interp.x2row     = griddedInterpolant(XGrid, 1:length(XGrid), ...
  'linear', 'linear');
TRI.interp.y2col     = griddedInterpolant(YGrid, 1:length(YGrid), ...
  'linear', 'linear');
TRI.hotStart.ZMat        = ZMat;
TRI.hotStart.unusedInd   = unusedInd;
TRI.hotStart.bndMat      = bndMat;
TRI.hotStart.nIterOffset = iterStatus.nIterOffset + iterStatus.nIter;
TRI.hotStart.tOffset     = iterStatus.tOffset + toc(iterStatus.tStart);

if (isnumeric(options.Verbose) || strcmpi(options.Verbose, 'final'))
  % Display
  displayFcn('final', iterStatus)
end

% Plots
plotFcn(TRI, iterStatus, options)

% Output
if (nargout == 1)
  varargout{1} = TRI;
end



function iterStatus = chkConvFcn(iterStatus, options)
% Check if any convergence criteria was fulfilled

% Default output
iterStatus.conv    = 0;
iterStatus.convMsg = {'Not converged'};

if (isempty(iterStatus.errors.errVec))
  % Complete triangulation of grid
  iterStatus.conv    = 1;
  iterStatus.convMsg = {'Grid completely triangulated'};
  return
end

% Calculate absolute errors
absErr = iterStatus.errors.absErr;
iterStatus.errors.maxAbs    = max(absErr);
iterStatus.errors.meanAbs   = mean(absErr);
iterStatus.errors.medianAbs = median(absErr);
% Percentile
absErrSort = sort(absErr);
vecLen     = length(absErrSort);
level      = 0.95;
if (~isempty(options.AbsPctileTol))
  level = options.AbsPctileTol(2)/100;
end
ind  = max(1, ceil(vecLen*level));
iterStatus.errors.AbsPctile    = absErrSort(ind);
iterStatus.errors.AbsPctileLev = level;

if (iterStatus.relErrFlag)
  % Calculate relative errors
  relErr = iterStatus.errors.relErr;
  iterStatus.errors.maxRel    = max(relErr);
  iterStatus.errors.meanRel   = mean(relErr);
  iterStatus.errors.medianRel = median(relErr);
  % Percentile
  relErrSort = sort(relErr);
  vecLen     = length(relErrSort);
  level      = 0.95;
  if (~isempty(options.RelPctileTol))
    level = options.RelPctileTol(2)/100;
  end
  ind  = max(1, ceil(vecLen*level));
  iterStatus.errors.RelPctile    = relErrSort(ind);
  iterStatus.errors.RelPctileLev = level;
end

% Check iteration parameters for convergence
if (~isempty(options.MaxTriangles) && iterStatus.nTri > options.MaxTriangles)
  % Maximum number of triangles exceeded
  iterStatus.conv    = 1;
  iterStatus.convMsg = {'Maximum number of triangles exceeded'};
  return
end
if (~isempty(options.MaxVertices) && iterStatus.nPts > options.MaxVertices)
  % Maximum number of vertices exceeded
  iterStatus.conv    = 1;
  iterStatus.convMsg = {'Maximum number of vertices exceeded'};
  return
end
if (~isempty(options.MaxIter) && iterStatus.nIterOffset + iterStatus.nIter ...
    > options.MaxIter)
  % Maximum number of iterations exceeded
  iterStatus.conv    = 1;
  iterStatus.convMsg = {'Maximum number of iterations exceeded'};
  return
end
if (~isempty(options.MaxTime) && iterStatus.tOffset + ...
    toc(iterStatus.tStart) > options.MaxTime)
  % Maximum run time exceeded
  iterStatus.conv    = 1;
  iterStatus.convMsg = {'Maximum run time exceeded'};
  return
end

% Check absolute error tolerances for convergence
convChkVec = [];
convChkMsg = {};
if (~isempty(options.MaxAbsTol))
  % Max absolute error
  if (iterStatus.errors.maxAbs <= options.MaxAbsTol)
    convChkMsg = [convChkMsg, {'Max absolute error less than set tolerance'}];
    convChkVec = [convChkVec; 1];
  else
    convChkVec = [convChkVec; 0];
  end
end
if (~isempty(options.MeanAbsTol))
  % Mean absolute error
  if (iterStatus.errors.meanAbs <= options.MeanAbsTol)
    convChkMsg = [convChkMsg, {'Mean absolute error less than set tolerance'}];
    convChkVec = [convChkVec; 1];
  else
    convChkVec = [convChkVec; 0];
  end
end
if (~isempty(options.MedianAbsTol))
  % Median absolute error
  if (iterStatus.errors.medianAbs <= options.MedianAbsTol)
    convChkMsg = [convChkMsg, ...
      {'Median absolute error less than set tolerance'}];
    convChkVec = [convChkVec; 1];
  else
    convChkVec = [convChkVec; 0];
  end
end
if (~isempty(options.AbsPctileTol))
  % Percentile
  if (iterStatus.errors.AbsPctile <= options.AbsPctileTol(1))
    convChkMsg = [convChkMsg, ...
      {'Absolute error percentile limit less than set tolerance'}];
    convChkVec = [convChkVec; 1];
  else
    convChkVec = [convChkVec; 0];
  end
end

if (iterStatus.relErrFlag)
  % Check relative error tolerances for convergence
  if (~isempty(options.MaxRelTol))
    % Max relative error
    if (iterStatus.errors.maxRel <= options.MaxRelTol)
      convChkMsg = [convChkMsg, ...
        {'Max relative error less than set tolerance'}];
      convChkVec = [convChkVec; 1];
    else
      convChkVec = [convChkVec; 0];
    end
  end
  if (~isempty(options.MeanRelTol))
    % Mean relative error
    if (iterStatus.errors.meanRel <= options.MeanRelTol)
      convChkMsg = [convChkMsg, ...
        {'Mean relative error less than set tolerance'}];
      convChkVec = [convChkVec; 1];
    else
      convChkVec = [convChkVec; 0];
    end
  end
  if (~isempty(options.MedianRelTol))
    % Median relative error
    if (iterStatus.errors.medianRel <= options.MedianRelTol)
      convChkMsg = [convChkMsg, ...
        {'Median relative error less than set tolerance'}];
      convChkVec = [convChkVec; 1];
    else
      convChkVec = [convChkVec; 0];
    end
  end
  if (~isempty(options.RelPctileTol))
    % Percentile
    if (iterStatus.errors.RelPctile <= options.RelPctileTol(1))
      convChkMsg = [convChkMsg, ...
        {'Relative error percentile limit less than set tolerance'}];
      convChkVec = [convChkVec; 1];
    else
      convChkVec = [convChkVec; 0];
    end
  end
end

if (strcmpi(options.ConvCriteria, 'all'))
  % All set error tolerances must be fulfilled
  if (min(convChkVec) == 1)
    iterStatus.convMsg = convChkMsg;
    iterStatus.conv    = 1;
  end
else
  % Just one set error tolerances must be fulfilled
  if (max(convChkVec) == 1)
    iterStatus.convMsg = convChkMsg;
    iterStatus.conv    = 1;
  end
end
    


function displayFcn(type, iterStatus)
% Printout

if (strcmpi(type, 'final'))
  % Final printout, reason for convergence
  disp(' ')
  disp('Converged because:')
  for i = 1 : length(iterStatus.convMsg)
    disp(['  ', iterStatus.convMsg{i}])
  end
  disp(' ')
end

disp(['Iterations               : ', num2str(iterStatus.nIterOffset + ...
  iterStatus.nIter)])
disp(['Time                     : ', num2str(iterStatus.tOffset + ...
  toc(iterStatus.tStart)), ' s'])
switch type
  case 'iter'
    disp(['Total triangles          : ', num2str(iterStatus.nTri)])
  case 'final'  
    disp(['Interior triangles       : ', num2str(iterStatus.nTri)])
end
disp(['Vertices                 : ', num2str(iterStatus.nPts)])
disp(['Max Absolute Error       : ', num2str(iterStatus.errors.maxAbs)])
disp(['Mean Absolute Error      : ', num2str(iterStatus.errors.meanAbs)])
disp(['Median Absolute Error    : ', num2str(iterStatus.errors.medianAbs)])
disp(['Absolute Error Percentile: ', num2str(iterStatus.errors.AbsPctile), ...
  ' (', num2str(100*iterStatus.errors.AbsPctileLev), '%)'])
if (iterStatus.relErrFlag)
  disp(['Max Relative Error       : ', ...
    num2str(100*iterStatus.errors.maxRel), '%'])
  disp(['Mean Relative Error      : ', ...
    num2str(100*iterStatus.errors.meanRel), '%'])
  disp(['Median Relative Error    : ', ...
    num2str(100*iterStatus.errors.medianRel), '%'])
  disp(['Relative Error Percentile: ', num2str(...
    100*iterStatus.errors.RelPctile), ' (', num2str(...
    100*iterStatus.errors.RelPctileLev), '%)'])
end

if (strcmpi(type, 'final'))
  % Final printout, storage requirements
  disp(' ')
  gridBytes = 8*prod(iterStatus.ZMatSize);
  disp(['Grid requires ', num2str(prod(iterStatus.ZMatSize)), ' double ', ...
    'floating point numbers (', num2str(gridBytes), ' bytes)'])
  triBytes   = 4*2*iterStatus.nPts + 8*iterStatus.nPts + 4*3*iterStatus.nTri;
  bytesRatio = 100*triBytes/gridBytes;
  bytesRatio = round(bytesRatio*10)/10;
  disp(['Triangulation requires ', num2str(iterStatus.nPts), ' double ', ...
    'floating point numbers and ', num2str(2*iterStatus.nTri), ...
    ' unsigned integers (', num2str(triBytes) , ' bytes, ', ...
    num2str(bytesRatio), '% of grid)'])  
end
disp(' ')



function iterPlot(mode, iterStatus, options)
% Plot convergence, make sure the plot window can be in the background 

persistent cf ca1 ca2 ca3 oldIterStatus

switch mode
  
  case 'init'
    % Initialize figure
    cf = figure;
    set (cf, 'Units', 'normalized', 'Position', [0, 0, 0.35, 1]);
    if (~iterStatus.relErrFlag)
      ca1 = subplot(2, 1, 1);
      box on
      hold on
      title('Triangles')
      xlabel('Iteration')
      ca2 = subplot(2, 1, 2);
      box on
      hold on
      title('Absolute Error')
      xlabel('Iteration')
      ylabel('Error')
    else
      ca1 = subplot(3, 1, 1);
      box on
      hold on
      title('Triangles')
      xlabel('Iteration')
      ca2 = subplot(3, 1, 2);
      box on
      hold on
      title('Absolute Error')
      xlabel('Iteration')
      ylabel('Error')
      ca3 = subplot(3, 1, 3);
      box on
      hold on
      title('Relative Error')
      xlabel('Iteration')
      ylabel('Relative Error [%]')
    end
    
  case 'iter'
    % Iteration update
    if (~ishandle(cf))
      % The user might have closed the window
      return
    end
    
    winSize = 10;
    if (iterStatus.nIter > 0)
      % Plot number of triangles and vertices
      try
        % The user might have closed the window
        set(cf, 'CurrentAxes', ca1)
        
        % Plot number of triangles and vertices
        set(ca1, 'YLimMode', 'auto')
        lp = [];
        lp(1) = line([iterStatus.nIter-1, iterStatus.nIter], ...
          [oldIterStatus.nPts, iterStatus.nPts], 'Color', 'k', ...
          'Linestyle', '-', 'lineWidth', 1.5);
        lp(2) = line([iterStatus.nIter-1, iterStatus.nIter], ...
          [oldIterStatus.nTri, iterStatus.nTri], 'Color', 'm', ...
          'Linestyle', '-', 'lineWidth', 1.5);
        if (~isempty(options.MaxVertices))
          line([iterStatus.nIter-1, iterStatus.nIter], ...
            [options.MaxVertices , options.MaxVertices ], 'Color', 'k', ...
            'Linestyle', ':', 'lineWidth', 1)
        end
        if (~isempty(options.MaxTriangles))
          line([iterStatus.nIter-1, iterStatus.nIter], ....
            [options.MaxTriangles, options.MaxTriangles], 'Color', 'm', ...
            'Linestyle', ':', 'lineWidth', 1)
        end
        if (iterStatus.nIter == 1)
          legend(lp, {'Vertices', 'Triangles'}, 'location', 'northwest')
        end
        xlim([0, iterStatus.nIter])
        yl = ylim;
        ylim([0, yl(2)])
      catch
        return
      end
      
      try
        % The user might have closed the window
        set(cf, 'CurrentAxes', ca2)
        
        % Plot absolute errors
        set(ca2, 'YLimMode', 'auto')
        lp = [];
        lp(1) = line([iterStatus.nIter-1, iterStatus.nIter], ...
          [oldIterStatus.errors.maxAbs   , iterStatus.errors.maxAbs   ], ...
          'Color', 'r', 'Linestyle', '-', 'lineWidth', 1.5);
        lp(2) = line([iterStatus.nIter-1, iterStatus.nIter], ...
          [oldIterStatus.errors.meanAbs  , iterStatus.errors.meanAbs  ], ...
          'Color', 'b', 'Linestyle', '-', 'lineWidth', 1.5);
        lp(3) = line([iterStatus.nIter-1, iterStatus.nIter], ...
          [oldIterStatus.errors.medianAbs, iterStatus.errors.medianAbs], ...
          'Color', 'g', 'Linestyle', '-', 'lineWidth', 1.5);
        lp(4) = line([iterStatus.nIter-1, iterStatus.nIter], ...
          [oldIterStatus.errors.AbsPctile, iterStatus.errors.AbsPctile], ...
          'Color', 'c', 'Linestyle', '-', 'lineWidth', 1.5);
        if (~isempty(options.MaxAbsTol))
          line([iterStatus.nIter-1, iterStatus.nIter], ...
            [options.MaxAbsTol   , options.MaxAbsTol   ], ...
            'Color', 'r', 'Linestyle', ':', 'lineWidth', 1)
        end
        if (~isempty(options.MeanAbsTol))
          line([iterStatus.nIter-1, iterStatus.nIter], ...
            [options.MeanAbsTol  , options.MeanAbsTol  ], ...
            'Color', 'b', 'Linestyle', ':', 'lineWidth', 1)
        end
        if (~isempty(options.MedianAbsTol))
          line([iterStatus.nIter-1, iterStatus.nIter], ...
            [options.MedianAbsTol, options.MedianAbsTol], ...
            'Color', 'g', 'Linestyle', ':', 'lineWidth', 1)
        end
        if (~isempty(options.AbsPctileTol))
          line([iterStatus.nIter-1, iterStatus.nIter], ...
            [options.AbsPctileTol, options.AbsPctileTol], ...
            'Color', 'c', 'Linestyle', ':', 'lineWidth', 1)
        end
        if (iterStatus.nIter == 1)
          legend(lp, {'Max', 'Mean', 'Median', ['Percentile (', ...
            num2str(100*iterStatus.errors.AbsPctileLev), '%)']}, ...
            'location', 'northeast')
        end
        xlim([max(0, iterStatus.nIter - winSize), iterStatus.nIter])
        yl = ylim;
        ylim([0, yl(2)])
      catch
        return
      end
      
      if (iterStatus.relErrFlag)
        try
          % The user might have closed the window
          set(cf, 'CurrentAxes', ca3)
          
          % Plot relative errors
          set(ca3, 'YLimMode', 'auto')
          lp = [];
          lp(1) = line([iterStatus.nIter-1, iterStatus.nIter], ...
            [100*oldIterStatus.errors.maxRel   , ...
            100*iterStatus.errors.maxRel   ], ...
            'Color', 'r', 'Linestyle', '-', 'lineWidth', 1.5);
          lp(2) = line([iterStatus.nIter-1, iterStatus.nIter], ...
            [100*oldIterStatus.errors.meanRel  , ...
            100*iterStatus.errors.meanRel  ],  ...
            'Color', 'b', 'Linestyle', '-', 'lineWidth', 1.5);
          lp(3) = line([iterStatus.nIter-1, iterStatus.nIter], ...
            [100*oldIterStatus.errors.medianRel, ...
            100*iterStatus.errors.medianRel],  ...
            'Color', 'g', 'Linestyle', '-', 'lineWidth', 1.5);
          lp(4) = line([iterStatus.nIter-1, iterStatus.nIter], ...
            [100*oldIterStatus.errors.RelPctile, ...
            100*iterStatus.errors.RelPctile], ...
            'Color', 'c', 'Linestyle', '-', 'lineWidth', 1.5);
          if (~isempty(options.MaxRelTol))
            line([iterStatus.nIter-1, iterStatus.nIter], ...
              [100*options.MaxRelTol   , 100*options.MaxRelTol   ], ...
              'Color', 'r', 'Linestyle', ':', 'lineWidth', 1)
          end
          if (~isempty(options.MeanRelTol))
            line([iterStatus.nIter-1, iterStatus.nIter], ...
              [100*options.MeanRelTol  , 100*options.MeanRelTol  ], ...
              'Color', 'b', 'Linestyle', ':', 'lineWidth', 1)
          end
          if (~isempty(options.MedianRelTol))
            line([iterStatus.nIter-1, iterStatus.nIter], ...
              [100*options.MedianRelTol, 100*options.MedianRelTol], ...
              'Color', 'g', 'Linestyle', ':', 'lineWidth', 1)
          end
          if (~isempty(options.RelPctileTol))
            line([iterStatus.nIter-1, iterStatus.nIter], ...
              [100*options.RelPctileTol, 100*options.RelPctileTol], ...
              'Color', 'c', 'Linestyle', ':', 'lineWidth', 1)
          end
          if (iterStatus.nIter == 1)
            legend(lp, {'Max', 'Mean', 'Median', ['Percentile (', ...
              num2str(100*iterStatus.errors.RelPctileLev), '%)']}, ...
              'location', 'northeast')
          end
          xlim([max(0, iterStatus.nIter - winSize), iterStatus.nIter])
          yl = ylim;
          ylim([0, yl(2)])
        catch
          return
        end
      end
      
    end
    drawnow
    
    % Update
    oldIterStatus = iterStatus;
    
  case 'final'
    % Iterations finished, just close the window
    if (ishandle(cf))
      delete(cf)
    end
    
end



function plotFcn(TRI, iterStatus, options)
% Plots possibly done after convergence

ptsList  = TRI.DT.Points;
if (iterStatus.constFlag)
  tf       = isInterior(TRI.DT);
  triList  = TRI.DT.ConnectivityList(tf, :);
else
  triList  = TRI.DT.ConnectivityList;
end
[X, Y]   = ndgrid(TRI.XGrid, TRI.YGrid);
ZMat     = TRI.hotStart.ZMat;
bndMat   = TRI.hotStart.bndMat;
% Evaluate triangulation at vertices
ptsInd   = sub2ind(size(ZMat), ptsList(:, 1), ptsList(:, 2));
[triInd, bCenter] = pointLocation(TRI.DT, ptsList);
triVals  = ZMat(ptsInd(TRI.DT(triInd, :)));
ZTriEval = dot(bCenter', triVals')';

% Plot surfaces
if (strcmpi(options.PlotSurf, 'on'))
  figure  
  subplot(1, 2, 1)
  title('GRID')
  hold on
  grid on
  sh1 = surf(X, Y, ZMat);
  shading interp
  set(sh1, 'edgecolor', 'k', 'edgealpha', 0.25)
  view(3)
  lightangle(45, 45)
  axis tight
  
  subplot(1, 2, 2)
  title('TRIANGULATION')
  hold on
  grid on
  sh2 = trisurf(triList, iterStatus.row2x(ptsList(:, 1)), ...
    iterStatus.col2y(ptsList(:, 2)), ZTriEval);
  shading interp
  set(sh2, 'edgecolor', 'k', 'edgealpha', 0.25)
  view(3)
  lightangle(45, 45)
  axis tight
end

% Plot triangulation
if (strcmpi(options.PlotTri, 'on'))
  figure
  hold on
  box on
  sh1 = surf(X, Y, ZMat - max(ZMat(:)));
  set(sh1, 'edgecolor', 'none')
  shading interp
  triplot(triList, iterStatus.row2x(ptsList(:, 1)), ...
    iterStatus.col2y(ptsList(:, 2)), 'k-');
  
  % Plot boundaries
  if (iterStatus.constFlag)
    segInd = unique(bndMat(:, 1));
    for segNo = segInd(:)'
      segMat = bndMat(bndMat(:, 1) == segNo, :);
      bndSeg = segMat(segMat(:, 3) == 1, :);
      [plotRow, plotCol] = ind2sub(size(ZMat), bndSeg(:, 2));
      plot(iterStatus.row2x(plotRow), iterStatus.col2y(plotCol), ...
        'r', 'lineWidth', 1)
    end
  end
  
  view(2)
  axis image
end

% Plot error distribution
if (strcmpi(options.PlotErr, 'on'))
  figure
  if (iterStatus.relErrFlag)
    subplot(1, 2, 1)
  else
    subplot(1, 1, 1)
  end
  hh1 = histogram(iterStatus.errors.absErr);
  set(hh1, 'FaceColor', 'b', 'FaceAlpha', 0.25)
  hold on
  title('Grid Points Without Vertices')
  xlabel('Absolute Error')
  yl = ylim;
  lp = zeros(1, 4);
  lp(1) = plot([iterStatus.errors.maxAbs   , iterStatus.errors.maxAbs   ], ...
    [yl(1), yl(2)], 'r', 'lineWidth', 1.5);
  lp(2) = plot([iterStatus.errors.meanAbs  , iterStatus.errors.meanAbs  ], ...
    [yl(1), yl(2)], 'b', 'lineWidth', 1.5);
  lp(3) = plot([iterStatus.errors.medianAbs, iterStatus.errors.medianAbs], ...
    [yl(1), yl(2)], 'g', 'lineWidth', 1.5);
  lp(4) = plot([iterStatus.errors.AbsPctile, iterStatus.errors.AbsPctile], ...
    [yl(1), yl(2)], 'c', 'lineWidth', 1.5);
  legend(lp, {'Max', 'Mean', 'Median', ['Percentile (', num2str(...
    100*iterStatus.errors.AbsPctileLev), '%)']}, 'location', 'northeast')
  
  if (iterStatus.relErrFlag)
    subplot(1, 2, 2)
    hh2 = histogram(100*iterStatus.errors.relErr);
    set(hh2, 'FaceColor', 'r', 'FaceAlpha', 0.25)
    hold on
    title('Grid Points Without Vertices')
    xlabel('Relative Error [%]')
    yl = ylim;
    lp = zeros(1, 4);
    lp(1) = plot([100*iterStatus.errors.maxRel   , ...
      100*iterStatus.errors.maxRel   ], [yl(1), yl(2)], 'r', 'lineWidth', 1.5);
    lp(2) = plot([100*iterStatus.errors.meanRel  , ...
      100*iterStatus.errors.meanRel  ], [yl(1), yl(2)], 'b', 'lineWidth', 1.5);
    lp(3) = plot([100*iterStatus.errors.medianRel, ...
      100*iterStatus.errors.medianRel], [yl(1), yl(2)], 'g', 'lineWidth', 1.5);
    lp(4) = plot([100*iterStatus.errors.RelPctile, ...
      100*iterStatus.errors.RelPctile], [yl(1), yl(2)], 'c', 'lineWidth', 1.5);
    legend(lp, {'Max', 'Mean', 'Median', ['Percentile (', num2str(...
      100*iterStatus.errors.RelPctileLev), '%)']}, 'location', 'northeast')
  end
end



function bndMat = discretizeBoundary(Z)
% Make a boundary matrix if the supplied grid contains NaN's. Contains all
% points on the boundary:
% First column is the boundary segment number
% Second column is the linear point index
% Third column is a flag indicating if the point is used as a boundary point in
% the current triangulation

% Add outer boundary
testMat = NaN*ones(size(Z, 1) + 2, size(Z, 2) + 2);
testMat(2:end-1, 2:end-1) = Z;

% Get contours
testMat(~isnan(testMat)) = 1;
testMat( isnan(testMat)) = 0;
contourMat = contourc(testMat, [0.95, 0.95]);

% Rip out contour segments
boundary   = [];
count      = 1;
contNr     = 0;
totElem    = 0;
while (1)
  numElem   = contourMat(2, count);
  row       = contourMat(2, count+1 : count+numElem) - 1;
  col       = contourMat(1, count+1 : count+numElem) - 1;
  rowColMat = [row(1:end-1)', col(1:end-1)'];
  rowColMat = round(rowColMat);
  % Make all points within index ranges
  rowColMat(rowColMat(:, 1) < 1) = 1;
  rowColMat(rowColMat(:, 1) > size(Z, 1)) = size(Z, 1);
  rowColMat(rowColMat(:, 2) < 1) = 1;
  rowColMat(rowColMat(:, 2) > size(Z, 2)) = size(Z, 2);
  % Remove points outside boundary
  zMat  = Z(sub2ind(size(Z), rowColMat(:, 1), rowColMat(:, 2)));
  rowColMat = rowColMat(~isnan(zMat), :);
  % Remove double points
  rowColMat = unique(rowColMat, 'rows', 'stable');
  totElem = totElem + size(rowColMat, 1);
  % Set boundary data
  contNr  = contNr + 1;
  boundary(contNr).ind = sub2ind(size(Z), rowColMat(:, 1), rowColMat(:, 2));
  boundary(contNr).row = rowColMat(:, 1);
  boundary(contNr).col = rowColMat(:, 2);
    
  % Reduce number of boundary points
  [rowColMatRed, keepInd] = dpsimplify(rowColMat, 0);
  % Set boundary data
  boundary(contNr).indRed = sub2ind(size(Z), rowColMat(keepInd, 1), ...
    rowColMat(keepInd, 2));
  boundary(contNr).rowRed = rowColMatRed(:, 1);
  boundary(contNr).colRed = rowColMatRed(:, 2);
  count = count + numElem + 1;
  if (count > size(contourMat, 2))
    break
  end
end

% Make boundary matrix
rowCount = 1;
bndMat   = zeros(totElem, 3);
for segNo = 1 : length(boundary)
  segInd     = boundary(segNo).ind;
  usedSegInd = boundary(segNo).indRed;
  segStatus  = ismember(segInd, usedSegInd);
  bndMat(rowCount:rowCount + length(segInd)-1, 1) = segNo;
  bndMat(rowCount:rowCount + length(segInd)-1, 2) = segInd(:);
  bndMat(rowCount:rowCount + length(segInd)-1, 3) = segStatus(:);
  rowCount = rowCount + length(segInd);
end



function [constMat, bndPts] = bndMat2const(bndMat, ZMatSize)
% Make constraints for triangulation from a boundary matrix

segRowVec = [];
segColVec = [];
constMat  = [];
segInd = unique(bndMat(:, 1));
for segNo = segInd(:)'
  segMat = bndMat(bndMat(:, 1) == segNo, :);
  bndSeg = segMat(segMat(:, 3) == 1, :);
  [tmpRow, tmpCol] = ind2sub(ZMatSize, bndSeg(:, 2));
  % Constraint matrix
  boundaryLen = length(tmpRow);
  startInd    =  length(segRowVec)+1 : length(segRowVec) + boundaryLen;
  stopInd     = [length(segRowVec)+2 : length(segRowVec) + boundaryLen, ...
    length(segRowVec)+1];
  constMat  = [constMat; startInd(:), stopInd(:)];
  % Boundary points
  segRowVec = [segRowVec; tmpRow];
  segColVec = [segColVec; tmpCol];
end
bndPts = [segRowVec, segColVec];



function [X, Y, Z, options] = parseInput(input)
% Parse & check input

% Default output
X       = [];
Y       = [];
Z       = [];
% Default options
options = trioptions('init', [], [], []);

% Get possible input options
optionsStart = find(cellfun(@ischar, input), 1, 'first');
optionsList  = input(optionsStart:end);
if (~isempty(optionsStart))
  input      = input(1:optionsStart-1);
end
nInput       = length(input);

if (nInput == 0)
  % Default demo mode
  [X, Y, Z] = peaks(49);
  Z = Z';
  X = X(1, :)';
  Y = Y(:, 1);
  options.MaxAbsTol = 0.01*(max(Z(:)) - min(Z(:)));
  options.PlotConv  = 'on';
  options.PlotSurf  = 'on';
  options.PlotTri   = 'on';
  options.PlotErr   = 'on';
  options.Verbose   = 10;

elseif (nInput == 1)
  if (isInt(input{1}, 2))
    % Integer supplied
    [X, Y, Z] = peaks(input{1});
    Z = Z';
    X = X(1, :)';
    Y = Y(:, 1);
    options.MaxAbsTol = 0.01*(max(Z(:)) - min(Z(:)));
    options.PlotConv  = 'on';
    options.PlotSurf  = 'on';
    options.PlotTri   = 'on';
    options.PlotErr   = 'on';
    options.Verbose   = 10;
    
  elseif (isnumeric(input{1}) && ismatrix(input{1}) && ...
      ~isvector(input{1}) && isreal(input{1}))
    % A matrix supplied
    Z = input{1};
    X = (1:size(input{1}, 1))';
    Y = (1:size(input{1}, 2))';
    options.MaxAbsTol = 0.01*(max(Z(:)) - min(Z(:)));
    
  elseif (isstruct(input{1}))
    % Possibly hot-start
    X = input{1};
    try
      Z = X.Z;
      options.MaxAbsTol = 0.01*(max(Z(:)) - min(Z(:)));
    catch
      error('Input is not a valid TRI structure!')
    end
    
  else
    error(['Input should be an integer > 1, a real two-dimensional ', ...
      'matrix or a TRI structure!'])
  end
  
elseif (nInput == 3)
  % X, Y, Z supplied
 XTmp = input{1};
 YTmp = input{2};
 ZTmp = input{3};
 if (~isnumeric(XTmp) || ~isnumeric(YTmp) || ~isnumeric(ZTmp) || ...
     ~isreal(XTmp) || ~isreal(YTmp) || ~isreal(ZTmp))
   error('Input should be real matrices and/or vectors!')
 end
 if (isvector(ZTmp) || ~ismatrix(ZTmp))
   error('Z should be a real two-dimensional matrix!')
 end
 [ZRows, ZCols] = size(ZTmp);
 % X and/or Y matrices
 if (~isvector(XTmp) && ~ismatrix(XTmp))
    error('X should be a real two-dimensional matrix!')
 end
 if (~isvector(YTmp) && ~ismatrix(YTmp))
    error('Y should be a real two-dimensional matrix!')
 end
  if (~isvector(XTmp) && size(XTmp, 1) ~= ZRows && size(XTmp, 2) ~= ZCols)
   error(['X should be a real matrix with ', num2str(ZRows), ' rows and ', ...
     num2str(ZCols), ' columns!'])
  end
 if (~isvector(YTmp) && size(YTmp, 1) ~= ZRows && size(YTmp, 2) ~= ZCols)
   error(['Y should be a real matrix with ', num2str(ZRows), ' rows and ', ...
     num2str(ZCols), ' columns!'])
 end
 
 % Possibly transpose Z
 if (~isvector(XTmp))
   % X is a matrix
   if (max(diff(XTmp(:, 1))) == 0)
     % Meshgrid-format needs to transpose!
     XTmp = XTmp';
     ZTmp = ZTmp';
     if (~isvector(YTmp))
       % Y is a matrix
       YTmp = YTmp';
       YVec = YTmp(1, :);
     else
       % Y is a vector
       YVec = YTmp;
     end
   else
     % ndgrid-format
     if (~isvector(YTmp))
       % Y is a matrix
       YVec = YTmp(1, :);
     else
       % Y is a vector
       YVec = YTmp;
     end
   end
   XVec = XTmp(:, 1);
 else
   % X is a vector
   XVec = XTmp;
   if (~isvector(YTmp))
     % Y is a matrix
     if (max(diff(YTmp(1, :))) == 0)
       % Meshgrid-format needs to transpose!
       YTmp = YTmp';
       ZTmp = ZTmp';
     end
     YVec = YTmp(1, :);
   else
     % Y is a vector
     YVec = YTmp;
   end  
 end
 [ZRows, ZCols] = size(ZTmp);
 
 % X and/or Y vectors, check lengths
 if (length(XVec) ~= ZRows)
   error(['X should be a real vector with ', num2str(ZRows), ' elements!'])
 end
 if (length(YVec) ~= ZCols)
   error(['Y should be a real vector with ', num2str(ZCols), ' elements!'])
 end
 
 % Sorting might be required
 [XVec, rowInd] = sort(XVec);
 [YVec, colInd] = sort(YVec);
 if (min(diff(XVec)) == 0)
   error('X contains duplicate values!')
 end
 if (min(diff(YVec)) == 0)
   error('Y contains duplicate values!')
 end
 ZTmp = ZTmp(rowInd, colInd);
 
 % Output
 X = XVec(:);
 Y = YVec(:);
 Z = ZTmp;
 options.MaxAbsTol = 0.01*(max(Z(:)) - min(Z(:)));
 
else
  % Wrong number of inputs
  error('Wrong number of input parameters!')
  
end

% Parse user options
options = trioptions('set', options, optionsList, Z);




function flag = isInt(value, lowLim, highLim)
% Check if the input is a real integer number, optionally within given bounds

% Default input
if (nargin < 3)
  highLim = inf;
end
if (nargin < 2)
  lowLim  = -inf;
end
if (isempty(lowLim)  ), lowLim   = -inf; end
if (isempty(highLim) ), highLim  =  inf; end

flag = 1;
if (~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
    rem(value, round(value)) ~= 0 || value < lowLim || value > highLim)
  flag = 0;
end