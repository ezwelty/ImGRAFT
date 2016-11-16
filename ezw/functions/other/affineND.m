% AFFINEND  Computes the affine transformation between two sets of ND
% points.
%
%   [s R T tA err] = affineND(A, B, doScale = 1)
%
% Computes the orientation and position (and optionally the uniform scale 
% factor) of the affine transformation between two corresponding point sets 
% Ai and Bi such that they are related by
%
%     Bi = s*R*Ai + T.
%
% Based on "A procedure for determining rigid body transformation
% parameters" by J.H. Challis (1995).
%
%
% Input:    A           DxN matrix of N points of dimension D
%           B           DxN matrix of N points of dimension D
%           doScale     boolean whether to allow scaling [default=1]
%
% Output:   s           uniform scale factor
%           R           DxD rotation matrix
%           T           Dx1 translation vector
%           tA          transformed A
%           err         1xN residual error:  ||Bi - tAi||

function [s, R, T, tA, err] = affineND(A, B, doScale)

% Check number of arguments
if (nargin < 2)
    error('Specify at least 2 input arguments.')
end

% Test size of point sets
[d1, n1] = size(A);
[d2, n2] = size(B);
if (d1 ~= d2)
    error('Point sets need to have the same number of dimensions.')
end
if (n1 ~= n2)
    error('Point sets need to have the same number of points.')
end
if (n1 < d1)
    error('Computation requires at least as many points as dimensions.')
end

% Enforce defaults
if (nargin < 3)
    doScale = 1;
end

% Put centroid of data points at origin (cA, cB)
normA = bsxfun(@minus,A,mean(A,2));
normB = bsxfun(@minus,B,mean(B,2));

% Compute rotation matrix (R)
c = (normB * normA') / size(normA,2); % correlation matrix
[u,~,v] = svd(c);   % singular value decomposition
R = u * v';
R = u * diag([ones(1,d1-1) det(R)],0) * v';    % enforce det(R) = 1

% Compute uniform scale factor (s)
if (doScale && n1 > 1)    
    a =0; b=0;
    for i=1:size(normA,2)
        a = a + normB(:,i)'*R*normA(:,i);
        b = b + normB(:,i)'*normB(:,i);
    end 
    s = b/a;    
else
    s = 1;
end

% Compute translation vector (T)
T = mean(B,2) - s*R*mean(A,2);

% Apply transformation (tA)
tA = bsxfun(@plus, s*R*A, T); 

% Return distance errors
err = sqrt(sum((B - tA).^2,1));