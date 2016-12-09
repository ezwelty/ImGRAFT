%Image and Terrain Modeling using Incremental Delaunay Triangulation
%
%Version : 0.1
%Date : 07.02.2010

function output_dt = aResize(input_dt,ratio)

vertex = input_dt.X;
vertex = round((vertex - 0.5) * ratio)+0.5;
output_dt = DelaunayTri(vertex);
