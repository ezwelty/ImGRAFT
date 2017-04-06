function out = dimfunc(func, data, dim, out_length, handles_matrix)
%dimfunc: apply a function along a specified dimension of an n-D array.
%
% The function need only handle vectors (or 2d matrices, if faster) and may
% return vector output. The size of the output matches that of the input,
% except that the chosen dimension changes to match the length of the
% output of the function (often, but not necessarily, 1).
%
% Usage: out = dimfunc(func, data, dim, out_length, handles_matrix)
% Input:
%  func - function handle (or function name string) to apply
%  data - input array, of arbitrary dimensionality and size
%  dim  - dimension to operate along, defaulting to first non-singleton dim
% Optional inputs:
%  out_length     - length of vector output from func (default: 1, scalar)
%  handles_matrix - Boolean flag, if true, assume func will operate on
%                   columns of a matrix input (default: false)
%
% Example:
%  v = randn(10,10,3);        % a 10x10 array of random 3-vectors
%  n = sqrt(sum(v.^2, 3));    % a 10x10 matrix of their Euclidean norms
%  N = dimfunc('norm', v, 3); % the same as n, but perhaps clearer
%
% Note that norm does not allow the "dim behaviour" that sum does; in this
% case, it is easy enough to use sum directly, but for other functions
% there might be no simple alternatives. I hope that dimfunc will make it
% slightly easier for others to write functions by removing the need to
% worry about n-D arrays, instead just handling either 1 or 2-d inputs.
%
% See also: shiftdata, shiftdim, reshape2matrix
%
% reshape2matrix by Carlos Adrian Vargas Aguilera is available here:
% http://www.mathworks.com/matlabcentral/fileexchange/25472
%
% Copyright 2010 Ged Ridgway
% http://www.mathworks.com/matlabcentral/fileexchange/authors/27434

%% Optional inputs
if nargin < 5 || isempty(handles_matrix)
    handles_matrix = false; % safer (though slower) to assume this
end
if nargin < 4 || isempty(out_length)
    out_length = 1;
end

%% Dimension (also optional -- but more commonly specified)
sz = size(data);
if nargin < 3 || isempty(dim)
    dim = find(sz > 1, 1, 'first'); % first non-singleton dimension
elseif dim < 1 || dim > length(sz) || dim ~= round(dim)
    warning('dimfunc:nodim', ...
        'Requested dim=%g but only 1 to %d available; returning data.', ...
        dim, length(sz));
    out = data;
    return
end

%% Convert function name string (like 'norm') to function handle (@norm)
if ischar(func)
    func = str2func(func);
end

%% Permute desired dimension into first place
perm = [dim 1:dim-1 dim+1:ndims(data)]; % (NB perm = 1:ndims if dim is 1)
if dim > 1
    data = permute(data, perm);
end

%% Reshape to matrix version
rows = sz(dim);
cols = numel(data) / rows;
if length(sz) > 2
    data = reshape(data, [rows cols]);
end

%% Apply function
if handles_matrix
    out = func(data);
else
    % Loop over column vectors
    out = nan(out_length, cols);
    for c = 1:cols
        out(:, c) = func(data(:, c));
    end
end

%% Revert output to correspond with original data
if length(sz) > 2
    out = reshape(out, [out_length sz(perm(2:end))]);
end
if dim > 1
    out = ipermute(out, perm);
end