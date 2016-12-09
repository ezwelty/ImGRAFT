%Image and Terrain Modeling using Incremental Delaunay Triangulation
%
%Version : 0.1
%Date : 07.02.2010


%Image Approxumation by Incremental Delaunay Triangulation
function [dt z maxerr] = aIDT(input_im,max_allowed_error,max_number_of_iterations)
% max_allowed_error and max_number_of_iterations 
% are stopping condition parameters

[m n] = size(input_im);
input_im = double(input_im);
[y x] = ndgrid(1:m,1:n);
x = x - 0.5;y = y - 0.5;
dt = DelaunayTri();
dt.X = [n-0.5 m-0.5;0.5 m-0.5;n-0.5 0.5;0.5 0.5];

maxerr = Inf;
iter = 0;
while(maxerr > max_allowed_error && iter < max_number_of_iterations)
    iter = iter + 1;    
    IntpFunc = TriScatteredInterp(dt, input_im((dt.X(:,1)-0.5)*m+dt.X(:,2)+0.5));
    interpolation = IntpFunc(x,y);    %this line is not efficient
    error = abs(input_im-interpolation);
    [maxrow rindex] = max(error);
    [maxerr cindex] = max(maxrow);    
    dt.X(end+1,:) = [cindex-0.5 rindex(cindex)-0.5];      
    fprintf('%f\n',maxerr);
end

for i=1:size(dt.X,1)
    z(i,1) = input_im(dt.X(i,2)+0.5,dt.X(i,1)+0.5);
end