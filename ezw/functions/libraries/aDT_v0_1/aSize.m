%Image and Terrain Modeling using Incremental Delaunay Triangulation
%
%Version : 0.1
%Date : 07.02.2010

function [m n] = aSize(dt)

n = max(dt.X(:,1))+0.5;
m = max(dt.X(:,2))+0.5;
