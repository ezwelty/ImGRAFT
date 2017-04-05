function options = trioptions(mode, options, optionsList, Z)

%TRIOPTIONS function for handling options for <a href="matlab:help grid2tri">grid2tri</a>. Not to be called, these
%options should be set in the call to <a href="matlab:help grid2tri">grid2tri</a>. Allowed parameter/value pairs
%are described below,
%
%CONVERGENCE CHECKING PARAMETERS
%
%  Name        : MaxTriangles
%  Value       : integer > 0 or empty (not used)
%  Default     : []
%  Description : Terminates if the supplied maximum number of triangles is
%                exceeded
%
%  Name        : MaxVertices
%  Value       : integer > 0 or empty (not used)
%  Default     : []
%  Description : Terminates if the supplied maximum number of vertices is
%                exceeded
%
%  Name        : MaxIter
%  Value       : integer > 0 or empty (not used)
%  Default     : []
%  Description : Terminates if the supplied maximum number of iterations is
%                exceeded
%
%  Name        : MaxTime
%  Value       : scalar > 0 or empty (not used)
%  Default     : []
%  Description : Terminates if the supplied maximum time spent in seconds is 
%                exceeded
%
%  Name        : MaxAbsTol
%  Value       : scalar > 0 or empty (not used)
%  Default     : 1% of the z-range
%  Description : The maximum absolute error tolerance in the triangulation
% 
%  Name        : MeanAbsTol
%  Value       : scalar > 0 or empty (not used)
%  Default     : []
%  Description : The mean absolute error tolerance in the triangulation
% 
%  Name        : MedianAbsTol
%  Value       : scalar > 0 or empty (not used)
%  Default     : []
%  Description : The median absolute error tolerance in the triangulation
% 
%  Name        : AbsPctileTol
%  Value       : scalar > 0, a 2 element vector where the first element is > 0
%                and the second > 0 and <= 100. If empty not used.
%  Default     : []
%  Description : The triangulation absolute error distribution 95% percentile
%                limit, i.e., 95% of the absolute errors must lie below the
%                supplied limit. Can optionally be a vector where the first
%                element is the percentile limit and the second is the used
%                percentile. For example [0.1, 99] means that 99% of the
%                absolute error must be below 0.1.
% 
%  Name        : MaxRelTol
%  Value       : scalar > 0 or empty (not used)
%  Default     : []
%  Description : The maximum relative error tolerance in the triangulation.
%                Here 0.05 means that relative error tolerance is 5%. Cannot be
%                used if the range of Z contains 0.
% 
%  Name        : MeanRelTol
%  Value       : scalar > 0 or empty (not used)
%  Default     : []
%  Description : The mean relative error tolerance in the triangulation.
%                Here 0.05 means that relative error tolerance is 5%. Cannot be
%                used if the range of Z contains 0.
% 
%  Name        : MedianRelTol
%  Value       : scalar > 0 or empty (not used)
%  Default     : []
%  Description : The median relative error tolerance in the triangulation.
%                Here 0.05 means that relative error tolerance is 5%. Cannot be
%                used if the range of Z contains 0.
%
%  Name        : RelPctileTol
%  Value       : scalar > 0, a 2 element vector where the first element is > 0
%                and the second > 0 and <= 100. If empty not used.
%  Default     : []
%  Description : The triangulation relative error distribution 95% percentile
%                limit, i.e., 95% of the relative errors must lie below the
%                supplied limit. Can optionally be a vector where the first
%                element is the percentile limit and the second is the used
%                percentile. For example [0.05, 99] means that 99% of the
%                relative errors must be below 0.05 (5%). Cannot be
%                used if the range of Z contains 0.
%
%  Name        : ConvCriteria
%  Value       : string 'all' or 'any'
%  Default     : 'all'
%  Description : If 'all' the triangulation terminates if all errors are less
%                than the corresponding supplied tolerances. If 'any' the
%                triangulation terminates if any of the errors is less than the
%                corresponding supplied tolerance.
%
%TRIANGULATION PARAMETERS
%
%  Name        : VertexSelection
%  Value       : a string 'abs' or 'rel'
%  Default     : 'abs'
%  Description : If 'abs' selects new vertices from the points with the greatest
%                absolute error. If 'rel' selects new vertices from the points
%                with the greatest relative error, cannot be used if the range
%                of Z contains 0.
%
%  Name        : VerticesAdded
%  Value       : a scalar > 0 and < 1 or an integer >= 1
%  Default     : 0.05
%  Description : If a scalar < 1 the number of new vertices added in each
%                iteration is the scalar*number of vertices. If an integer this
%                number of new vertices are added. Adding more vertices in each
%                iteration can result in faster convergence, but will yield more
%                triangles.
%
%  Name        : AddSeparateTri
%  Value       : a string 
%  Default     : 'on'
%  Description : If multiple points are added 'on' means that they are added to
%                different triangles, if 'off' (not recommended) they can be 
%                added to the same triangle.
%
%PRINTING AND PLOTTING PARAMETERS
%
%  Name        : Verbose
%  Value       : a string 'off' or 'final' or an integer >= 1
%  Default     : 'off'
%  Description : If 'off' no printing is done, if 'final' a summary is printed
%                after convergence and if an integer N status is printed every
%                N:th iteration.
%
%  Name        : PlotConv
%  Value       : a string 'on' or 'off'
%  Default     : 'off'
%  Description : If 'on' plots convergence during triangulation
%
%  Name        : PlotSurf
%  Value       : a string 'on' or 'off'
%  Default     : 'off'
%  Description : If 'on' plots the grid and the corresponding triangulated 
%                surface in 3D
%
%  Name        : PlotTri
%  Value       : a string 'on' or 'off'
%  Default     : 'off'
%  Description : If 'on' plots the triangulation in 2D
%
%  Name        : PlotErr
%  Value       : a string 'on' or 'off'
%  Default     : 'off'
%  Description : If 'on' plots the absolute and relative error distribution
%
%  VERSION 1.0
%  Author: Patrik Forssén 2015
%
%   See also <a href="matlab:help grid2tri">grid2tri</a>
%

switch mode
  case 'init'
    % Default options
    options.MaxTriangles    = [];
    options.MaxVertices     = [];
    options.MaxIter         = [];
    options.MaxTime         = [];
    options.MaxAbsTol       = [];
    options.MeanAbsTol      = [];
    options.MedianAbsTol    = [];
    options.AbsPctileTol    = [];
    options.MaxRelTol       = [];
    options.MeanRelTol      = [];
    options.MedianRelTol    = [];
    options.RelPctileTol    = [];
    options.ConvCriteria    = 'all';
    options.VertexSelection = 'abs';
    options.VerticesAdded   = 0.05;
    options.AddSeparateTri  = 'on';
    options.Verbose         = 'off';
    options.PlotConv        = 'off';
    options.PlotSurf        = 'off';
    options.PlotTri         = 'off';
    options.PlotErr         = 'off';
    
  case 'set'
    % Parse user options
    if (rem(length(optionsList), 2))
      error('Invalid parameter/value list!')
    end
    optionsList = reshape(optionsList, 2, length(optionsList)/2);
    validParam  = fieldnames(options);
    for i = 1 : size(optionsList, 2)
      paramName  = optionsList{1, i};
      paramValue = optionsList{2, i};
      if (~ischar(paramName))
        error('Invalid parameter/value list!')
      end
      paramInd = find(ismember(upper(validParam), upper(paramName)));
      if (isempty(paramInd))
        error(['Invalid parameter name ''', paramName, '''!'])
      end
      paramName = validParam{paramInd};
      switch paramName
        case {'MaxTriangles', 'MaxVertices', 'MaxIter'}
          if (~isempty(paramValue) && ~isInt(paramValue, 2))
            error(['The parameter ''', paramName, ...
              ''' should be empty or an integer > 1!'])
          end
        case {'MaxTime', 'MaxAbsTol', 'MeanAbsTol', 'MedianAbsTol', ...
            'MaxRelTol', 'MeanRelTol', 'MedianRelTol'}
          if (~isempty(paramValue) && ~isFloat(paramValue, eps))
            error(['The parameter ''', paramName, ...
              ''' should be empty or a scalar > 0!'])
          end
        case {'AbsPctileTol', 'RelPctileTol'}
          if (~isempty(paramValue))
            if (isnumeric(paramValue))
              if ((numel(paramValue) == 1 && isFloat(paramValue, eps)) || ...
                  (numel(paramValue) == 2 && isFloat(paramValue(1), eps) && ...
                  isFloat(paramValue(2), eps, 100)))
              else
                error(['The parameter ''', paramName, ''' should be empty', ...
                  ', a scalar > 0 or a 2 element vector where the first', ...
                  ' element is > 0 and the second > 0 and <= 100!'])
              end
            else
              error(['The parameter ''', paramName, ''' should be empty', ...
                ', a scalar > 0 or a 2 element vector where the first', ...
                ' element is > 0 and the second > 0 and <= 100!'])
            end
            if (numel(paramValue) == 1)
              paramValue = [paramValue, 95];
            end
          end
        case 'ConvCriteria'
          if (~ischar(paramValue) || (~strcmpi(paramValue, 'all') && ...
              ~strcmpi(paramValue, 'any')))
            error(['The parameter ''ConvCriteria'' should be a ', ...
              'string ''all'' or ''any''!'])
          end
        case 'VertexSelection'
          if (~ischar(paramValue) || (~strcmpi(paramValue, 'abs') && ...
              ~strcmpi(paramValue, 'rel')))
            error(['The parameter ''VertexSelection'' should be a ', ...
              'string ''abs'' or ''rel''!'])
          end
        case 'VerticesAdded'
          if (~isFloat(paramValue, eps, 1-eps) && ~isInt(paramValue, 1))
            error(['The parameter ''VerticesAdded'' should be a scalar', ...
              ' > 0 and < 1 or an integer >= 1'])
          end
        case 'Verbose'
          if ((~ischar(paramValue) && ~isInt(paramValue, 1)) || ...
              (ischar(paramValue) && ~strcmpi(paramValue, 'off') && ...
              ~strcmpi(paramValue, 'final')))
            error(['The parameter ''Verbose'' should be a string ''off''', ...
              ' or ''final'' or an integer >= 1!'])
          end
        case {'AddSeparateTri', 'PlotConv', 'PlotSurf', 'PlotTri'}
          if (~ischar(paramValue) || (~strcmpi(paramValue, 'on') && ...
              ~strcmpi(paramValue, 'off')))
            error(['The parameter ''', paramName, ...
              ''' should be a string ''on'' or ''off''!'])
          end
      end
      % Set parameter
      options.(paramName) = paramValue;
    end
    
    % Check if relative error is allowed
    if (min(Z(:)) <= 0 && max(Z(:)) >= 0)
      if (~isempty(options.MaxRelTol) || ~isempty(options.MeanRelTol) || ...
          ~isempty(options.MedianRelTol) || ...
          ~isempty(options.RelPctileTol) || ...
          strcmpi(options.VertexSelection, 'rel'))
        error([...
          'Z range contains zero, setting relative error tolerances or ', ...
          'using vertex selection based on relative error is not possible!'])
      end
    end
    
end



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



function flag = isFloat(value, lowLim, highLim, cplxFlag)
% Check if the input is a floating point number, optionally within given bounds
% and optionally allow complex value

% Default input
if (nargin < 4)
  cplxFlag = 0;
end
if (nargin < 3)
  highLim = inf;
end
if (nargin < 2)
  lowLim  = -inf;
end
if (isempty(lowLim)  ), lowLim   = -inf; end
if (isempty(highLim) ), highLim  =  inf; end
if (isempty(cplxFlag)), cplxFlag = 0   ; end

flag = 1;
if (~isnumeric(value) || ~isscalar(value) || value < lowLim || value > highLim)
  flag = 0;
end
if (~cplxFlag && ~isreal(value))
  flag = 0;
end