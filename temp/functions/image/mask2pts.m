% MASK2PTS  Extract image coordinates of points in a binary image mask.
%
%   impts = mask2pts(mask, method = 'skeleton')
%
% By default returns the image coordinates of the topological skeleton of
% the mask, useful for filtering noise from points and lines marked 
% manually on an image. Otherwise, returns the image coordinates of all the
% mask pixels. Points are returned sorted in increasing x, then increasing
% y.
%
% Input:    mask        binary image
%           method      either 'skeleton' or 'all'
%
% Output:   impts       NX2 matrix of image coordinates [x y]
%
% See also colormask, svg2struct.

function impts = mask2pts(mask, method)

% Enforce defaults
if (nargin < 2)
    method = 'skeleton';
end

% Extract point coordinates
switch method
    case 'all'
        [i,j] = find(mask);
    case 'skeleton'
        mask = bwmorph(mask, 'clean');
        skel = bwmorph(mask, 'thin', Inf);
        [i,j] = find(skel);
    otherwise
        error(['Unsupported method encountered: ' method])
end

% Set origin to [0 0] rather than Matlab's [0.5 0.5]
impts = bsxfun(@minus, [j i], [0.5 0.5]);