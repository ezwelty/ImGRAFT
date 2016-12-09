function tridemos(demoNo)

%TRIDEMOS some demos for <a href="matlab:help grid2tri">grid2tri</a>
%
%  TRIDEMOS(N) runs demo number N (default 1)
%  
%  Demo 1 shows how to evaluate a triangulation using <a href="matlab:help trival">trival</a>
%  Demo 2 shows how to use hot-start to refine a previously calculated
%         triangulation
%  Demo 3 is an example of a triangulation of a grid with undefined areas
%
%  VERSION 1.0
%  Author: Patrik Forssén 2015
%
%   See also <a href="matlab:help grid2tri">grid2tri</a>, <a href="matlab:help trival">trival</a>
%

if (nargin == 0)
  demoNo = 1;
end


switch demoNo
  
  case 1
    % Make a triangulation and evaluate it on a grid
    [X, Y, Z] = peaks(300);
    tri = grid2tri(X, Y, Z, ...
      'Verbose' , 'final' , ...
      'PlotConv', 'on'    , ...
      'PlotSurf', 'on');
    
    [XI, YI] = meshgrid(linspace(-2, 2), linspace(-2, 2));
    ZI = trival(tri, XI, YI);
    % Plot
    figure
    hold on
    grid on
    title('Triangulation Evaluated on a Grid')
    sh1 = surf(XI, YI, ZI);
    shading interp
    set(sh1, 'edgecolor', 'k', 'edgealpha', 0.25)
    view(3)
    lightangle(45, 45)
    
  case 2
    % Make a triangulation and use hot-start to refine it
    [X, Y, Z] = peaks(300);
    tri    = grid2tri(X, Y, Z, ...
      'Verbose' , 'final', ...
      'PlotConv', 'on'   , ...
      'PlotSurf', 'on');
    
    grid2tri(tri, ...
      'MaxAbsTol', 0.05   , ...
      'Verbose'  , 'final', ...
      'PlotConv' , 'on'   , ...
      'PlotSurf' , 'on');
    
  case 3
    % Make a triangulation where there is undefined areas in the supplied grid
    [X, Y, Z] = peaks(300);
    AN        = ones(size(X));
    R         = sqrt(X.^2 + Y.^2);
    AN(R > 3) = NaN;
    AN(R < 1) = NaN;
    ZAN       = Z.*AN;
    
    grid2tri(X, Y, ZAN    , ...
      'MaxAbsTol', 0.05   , ...
      'Verbose'  , 'final', ...
      'PlotConv' , 'on'   , ...
      'PlotTri'  , 'on'   , ...
      'PlotSurf' , 'on');
    
end