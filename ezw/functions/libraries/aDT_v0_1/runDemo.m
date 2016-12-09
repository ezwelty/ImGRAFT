%Image and Terrain Modeling using Incremental Delaunay Triangulation
%
%Version : 0.1
%Date : 07.02.2010

clc

I = imread('peak.256.png');
[m n] = size(I);
try
    I = rgb2gray(I);
end
[dt z maxerr] = aIDT(I,15,5000);
A = aReCon(dt,z);
figure;imshow(I,[]);xlabel('Original Image');
figure;imshow(A,[]);ylabel('Approximated Image');xlabel(['Absolute Error: ',num2str(maxerr)]);
figure;imshow(A,[]);xlabel('Approximated Image and Delaunay Triangulation');
hold on;triplot(dt,'r');xlabel('Delaunay Triangulation');axis([0 n 0 m]);

dt2 = aResize(dt,2.0);
A2 = aReCon(dt2,z);
figure;imshow(A2,[]);xlabel('Reconstructed Image x2');