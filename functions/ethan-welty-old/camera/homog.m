% HOMOG  Return the homogeneous equivalent of a matrix.
%
%   Xh = homog(X, pts)
%
% If X represents points (pts = 1, the default), a row of ones is appended 
% to the bottom of X:
%   
%   Xh = [ x11 x12 x13 ... 
%          x21 x22 x23 ...
%          x31 x32 x33 ...
%          ... ... ... ...
%           1   1   1  ... ]
%
% If X represents a transform (pts = 0), two outcomes are possible. If X is
% square, a row and a column of zeros are appended to the bottom and left 
% of X, with a one on the main diagonal:
%
%   Xh = [ x11 x12 x13 ...  0
%          x21 x22 x23 ...  0
%          x31 x32 x33 ...  0
%          ... ... ... ... ...
%           0   0   0  ...  1  ]
%
% If X is N-1xN, a row of zeros is appended to the bottom of X, with a one 
% on the main diagonal:
%
%   Xh = [ x11 x12 x13 ... x1N
%          x21 x22 x23 ... x2N
%          x31 x32 x33 ... x3N
%          ... ... ... ... ...
%           0   0   0  ...  1  ]
%
%
% Input:    X           DxN points matrix, NxN or N-1xN transform matrix
%           pts         boolean: points(1) or transform(0) [default=1]
%
% Output:   Xh          homogeneous equivalent of X
%
% See also inhomog.

function Xh = homog(X, pts)

% Check number of arguments
if (nargin < 1)
    error('Specify whether a points or transform matrix.')
end

% Enforce defaults
if (nargin < 2)
    pts = 1;
end

% Return homogeneous equivalent
[~, ncol] = size(X);

if (pts)
    
    % if points, append row [1 1 1 ... ]
    Xh = [X ; ones(1,ncol)];
    
else 
    
    % if transform, check dimensions of transform matrix
    if (nrow == ncol)
        
        % if NxN
        % append row [0 0 0 ... 1] and column [0 0 0 ... 1]'
        Xh = [X zeros(nrow,1) ; zeros(1,ncol) 1];
    
    elseif (nrow == ncol-1)
        
        % if N-1xN
        % append row [0 0 0 ... 1]
        Xh = [X ; zeros(1,ncol-1) 1];
        
    else
        error('Transform matrix must be either NxN or N-1xN.')
    end
    
end