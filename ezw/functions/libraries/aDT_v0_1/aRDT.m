%Image and Terrain Modeling using Incremental Delaunay Triangulation
%
%Version : 0.1
%Date : 07.02.2010


%Image Approximation by Random Delaunay Triangulation
function [dt z] = aRDT(input_im,number_of_vertex)

[m n] = size(input_im);
input_im = double(input_im);
vertex = [n-0.5 m-0.5;0.5 m-0.5;n-0.5   0.5;0.5   0.5;floor(rand(number_of_vertex-4,1)*n)+0.5 floor(rand(number_of_vertex-4,1)*m)+0.5];
dt = DelaunayTri(vertex(:,1),vertex(:,2));
for i=1:size(dt.X,1)
    z(i,1) = input_im(dt.X(i,2)+0.5,dt.X(i,1)+0.5);
end

