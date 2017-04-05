% COLORMASK  Filter image pixels by color.
% 
%   mask = colormask(img, col, tol = 0)
%
% Returns a binary mask of pixels within tolerance tol of the colors col.
%
%
% Input:    img         D-channel image
%           col         color vectors [NxD]
%           tol         color tolerance [scalar or NxD]
%
% Output:   mask        single-channel binary image
%
% See also mask2pts, svg2struct.

function mask = colormask(img, col, tol)

% Enforce defaults
if (nargin < 3)
    tol = 0;
end

% Dimensions
D = size(img,3);
N = size(col,1);
if length(tol) == 1
    tol = repmat(tol, N, D);
end

% Build mask
mask = zeros(size(img(:,:,1)));
for n = 1:N % color
    temp = ones(size(img(:,:,1)));
    for d = 1:D % channel
        temp = temp & (img(:,:,d) >= col(n,d)-tol(n,d)) & (img(:,:,d) <= col(n,d) + tol(n,d));
    end
    mask = mask | temp;
end