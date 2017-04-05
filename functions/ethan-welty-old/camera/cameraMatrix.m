% CAMERAMATRIX  Create a camera matrix (K) from linear camera parameters.
%
%   K = cameraMatrix(fpx, c)
% 
% Constructs camera matrix (K) according to the standard convention
%
%   K = [ fpx(1) 0      c(1) 
%         0      fpx(2) c(2) 
%         0      0      1   ]
%
% where fpx is the focal length in pixels and c is the principal point in pixels. 
%
% 
% Input:    f           focal length in pixels [fx (fy)]
%           c           principal point in pixels [cx cy]
%           
% Output:   K           camera matrix
%
% See also fmm2fpx, imgcenter.

function K = cameraMatrix(fpx, c)

% initialize matrix
K = zeros(3);
K(3,3) = 1;

% insert focal length (assume equal in x and y if length(fpx) == 1)
K(1,1) = fpx(1);
if length(fpx)>1
    K(2,2) = fpx(2);
else
    K(2,2) = fpx(1);
end

% insert principal point
K(1,3) = c(1);
K(2,3) = c(2);