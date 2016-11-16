% INHOMOG  Return the inhomogeneous equivalent of a matrix.
%
%   X = inhomog(Xh, pts)
%
% If X represents points (as columns, pts = 1), a row of ones is stripped
% from the bottom of X.
%
% If X represents a transform (pts = 0), a row of zeros, with a one on the
% main diagonal, is stripped from the bottom of X. If the rightmost column
% of the resulting matrix is only zeros, the rightmost column is also
% stripped.
%
%
% Input:    Xh          points or transform matrix
%           pts         boolean: points(1) or transform(0) [default=1]
%
% Output:   X           inhomogeneous equivalent of X
%
% See also homog.

function X = inhomog(Xh, pts)

% Check number of arguments
if (nargin < 1)
    error('Specify either a points or transform matrix.')
end

% Enforce defaults
if (nargin < 2)
    pts = 1;
end

% Return inhomogeneous equivalent
[~, ncol] = size(Xh);

if (pts)

    % if points, strip bottom row if [1 1 1 ... ]
    matches = Xh(end,:) == ones(1,ncol);
    if (sum(matches) == ncol)
        X = Xh(1:end-1,:);
    else
        error('Matrix is not in a homogeneous format.')
    end
    
else 
    
    % if transform, strip bottom row if [0 0 0 ... 1]
    matches = Xh(end,:) == [zeros(1,ncol-1) 1];
    if (sum(matches) == ncol)
        X = Xh(1:end-1,:);
        
        % and strip rightmost column if [0 0 0 ... ]'
        matches = X(:,end) == zeros(nrow-1,1);
        if (sum(matches) == nrow-1)
            X = X(:,1:end-1);
        end
    else
        error('Matrix is not in a homogeneous format.')
    end
    
end