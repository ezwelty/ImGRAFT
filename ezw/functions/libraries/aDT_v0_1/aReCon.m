%Image and Terrain Modeling using Incremental Delaunay Triangulation
%
%Version : 0.1
%Date : 07.02.2010

function output_im = aReCon(dt,z)

[m n] = aSize(dt);
[y x] = ndgrid(1:m,1:n);
x = x - 0.5;y = y - 0.5;
F = TriScatteredInterp(dt,z);
output_im = F(x,y);