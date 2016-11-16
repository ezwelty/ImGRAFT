% INTERPXN  Nonparametric interpolation in N dimensions.
%
%   Xi = interpXn(X, t, ti, type = 'linear', param = [])
%
% Interpolates to find Xi, the values of the variables (columns) in X
% specified at points t, at the query points ti. Currently assumes constant
% "velocity" between points t.
%
% NOTE: If method is not 'linear' or a kernel, but one from interp1, a path 
% is fitted to the points according to the method, resampled at density
% specified by param, then linearly interpolated at ti.
%
% Input:    X          NxD matrix of D variables
%           t          N vector of points corresponding to each row in X
%           ti         M vector of points at which to reinterpolate X 
%           type       'linear' - linear interpolation
%                      'unifk' - uniform kernel
%                      'bisqk' - bisquared kernel
%                      Any method from interp1
%           param      ignored for linear, kernel radius if kernel, and
%                      path density if method from interp1
%
% Output:   Xi         MXD matrix of variables X calculated at points ti

function Xi = interpXn(X, t, ti, type, param)

% Enforce defaults
if nargin < 4
    type = 'linear';
end
if ~strcmp(type,'linear') && nargin < 5
    error(['Specify the parameter value for method ' type '.'])
end

% Enforce t,ti are vertical
t = reshape(t,length(t),1);
ti = reshape(ti,length(ti),1);

% Sort reference positions in time
[t,ind] = sort(t);
X = X(ind,:);

% Initialize results matrix
Xi = nan(length(ti),size(X,2));

% Bisquared kernel
if strcmp(type, 'bisqk')
    radius = param;
    for n=1:size(X,2)
        for i=1:length(ti)
            [ival,weights] = bisqkernel(t,ti(i),radius);  % weights for measurement
            Xi(i,n)= sum(weights.*X(ival,n)) / sum(weights);
        end
    end

% Uniform kernel
elseif strcmp(type, 'unifk')
    radius = param;
    for n=1:size(X,2)
        for i=1:length(ti)
            [ival,weights] = unifkernel(t,ti(i),radius);  % weights for measurement
            Xi(i,n)= sum(weights.*X(ival,n)) / sum(weights);
        end
    end
   
% Linear interpolation    
elseif strcmp(type, 'linear')
    for i=1:length(ti)
        p1 = find(t <= ti(i),1,'last');
        if t(p1) == ti(i)
            Xi(i,:) = X(p1,:);
        elseif (~isempty(p1) && p1 < length(t))
            p2 = p1 + 1;
            vel = (X(p2,:) - X(p1,:)) / (t(p2) - t(p1));
            Xi(i,:) = X(p1,:) + vel * (ti(i) - t(p1));
        elseif (ti(i) < t(1))
            vel = (X(2,:) - X(1,:)) / (t(2) - t(1));
            Xi(i,:) = X(1,:) + vel * (ti(i) - t(1));
        elseif (ti(i) > t(end))
            vel = (X(end,:) - X(end-1,:)) / (t(end) - t(end-1));
            Xi(i,:) = X(end,:) + vel * (ti(i) - t(end));
        end 
    end

% Fit path to points, linear interpolation
else
    
    % densify points
    density = param; %(for now, just making this constant)
    [N, D] = size(X);
    Xhd = nan((N-1)*density+1,D);
    for d = 1:D
        Xhd(:,d) = interp1(1:N, X(:,d), linspace(1,N,(N-1)*density+1), type);
    end
    
    % compute cumulative chordal arclength
    L = [0 ; cumsum(sqrt(sum(diff(Xhd,1,1).^2,2)))];
    
    % recalculate t
    Lx = L(1:density:end); Lxm = repmat(Lx',density,1); L0 = Lxm(:); L0(end-density+1:end) = [];
    dL = diff(Lx); dLv = repmat(dL',density,1); dLv = dLv(:);
    t0 = repmat(t',density,1); t0 = t0(:); t0(end-density+1:end) = [];
    dt = diff(t); dtv = repmat(dt',density,1); dtv = dtv(:);
    thd = t0 + ((L(1:end-1) - L0) ./ dLv) .* (dtv); thd = [thd ; t(end)];
    
    % linear interpolation
    Xi = interpXn(Xhd,thd,ti,'linear');
end    

% Compute bisquared kernel of weights
function [ival,weights] = bisqkernel(t,ti,radius)
dist = sqrt((t-ti).^2);
ival = dist < radius;
weights = 15/16*(1-(dist(ival)/radius).^2).^2;

% Compute uniform kernel of weights 
% (=1 for all points within range)
function [ival,weights] = unifkernel(t,ti,radius)
dist = sqrt((t-ti).^2);
ival = dist < radius;
weights = ones(sum(ival),1);