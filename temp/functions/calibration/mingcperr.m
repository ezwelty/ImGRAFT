% MINGCPERR  Minimizes the mean reprojection error of ground control points 
% over multiple images by adjusting the specified camera parameters.
%
%   images = mingcperr(images, params, [options])
%
% Uses an optimization routine to minize the mean reprojection error of
% ground control points (paired image-world coordinates in required field 
% 'gcp') over multiple images by adjusting the camera parameters of the 
% images as initialized in params.
%
% Input:    images      N-camera structure with field 'gcp' [imx imy x y z]
%           params      Structure (length 1 or N) with any of following 
%                       fields (oblang, fpx, k, c, P) specifying the
%                       parameters to be adjusted and their initial
%                       values. If only the first element of a field is
%                       non-empty, it is applied to all images. Any NaN
%                       prevents that value from being adjusted. For
%                       example, params(1).P = [NA NA 2] allows only 
%                       images(1).P(3) to be adjusted.
%
%           options     Optimization options (see optimset)
%
% Output:   images      Updated images structure
%                       - updated fields (as specified by params)
%                       - reprojection errors: reprojerr [dx dy]
%
% FIXME: Instead return new params, and write new function to update images
% from a params structure.
% FIXME: Write helper functions for flattening and expanding structures.
%
% See also reprojerr

function images = mingcperr(images, params, options)

% Remove invalid fields from params
validfields = {'oblang', 'fpx', 'k', 'c', 'P'};
paramfields = fieldnames(params);
params = rmfield(params, paramfields(~ismember(paramfields,validfields)));
f = fieldnames(params);
if isempty(f)
    error('Params contains no valid fields (oblang, fpx, k, c, P).');
end

% As needed, add param fields to images
addfields = f(~ismember(f, fieldnames(images)));
if ~isempty(addfields)
    for i = addfields
        eval(['images(end).' char(i) ' = [];']);
    end
end

% Dimensions
Ni = length(images);
Np = length(params);

% As needed, populate images.oblang from images.R
if ismember('oblang', f) && ismember('R', fieldnames(images))
    for n = 1:Ni
        if ~isempty(images(n).R)
            images(n).oblang = rot2oblang(images(n).R);
        end
    end
end

% Flatten params to vector
% (all because fmin functions in MATLAB only optimize on double input)
X0 = [];
for i = 1:length(f)
  for n = 1:Np
    param = params(n).(f{i});
    image = images(n).(f{i});
    % Unless param empty
    if ~isempty(param)
        % if image not empty and different size, error
        if ~isempty(image) && ~isequal(size(param), size(image))
            error(['Unequal, non-empty dimensions at params(' n ').' f{i} ' => images(' n ').' f{i}]);
         % else, add to input
        else
            X0 = [X0 ; reshape(param, [], 1)];
        end
    end
  end
end
% Remove NaN from input
X0 = X0(~isnan(X0));

% Optimize
if (nargin < 3)
    options = optimset('Display', 'on', 'LargeScale', 'off');
end
[X, fval] = fminunc(@reprojerr_fmin, X0, options);
disp(['MEAN ERROR'])
disp(fval)

% Remove oblang field
if ismember('oblang', fieldnames(images))
    images = rmfield(images, 'oblang');
end

% (nested function)
function meanerr = reprojerr_fmin(X)
    
    % Expand vector
    for i = 1:length(f)
        applyall = true;
        for n = 1:Np
            param = params(n).(f{i});
            image = params(n).(f{i});
            % if empty, skip
            if isempty(param)
                continue;
            % otherwise, update image structure
            else
              nans = isnan(param);
              L = numel(param) - sum(nans(:));
              % if empty image, update all values
              if isempty(image)
                  % throw error if param contains nan
                  if any(nans)
                      error(['Param with NaN at empty image: params(' n ').' f{i} ' => images(' n ').' f{i}]);
                  else
                      new = reshape(X(1:L), size(param));
                  end
              % otherwise, update values skipping nans
              else
                  new = image;
                  new(~nans(:)) = X(1:L);
              end
              % apply changes
              images(n).(f{i}) = new;
              % Apply all only if 1st param non-empty
              if applyall && n > 1
                applyall = false;
              end
            end
            X = X((L+1):end);
        end
        % As needed, apply change to all images
        if applyall
          [images.(f{i})] = deal(new);
        end
    end
    
    % As needed, update images.R from images.oblang
    if ismember('oblang', f)
        for n = 1:Ni
            images(n).R = oblang2rot(images(n).oblang);
        end
    end
  
    % For each image
    Ni = length(images);
    err = nan(Ni,1);
    for n=1:Ni
        % reprojection errors
        images(n).gcperr = (reprojerr(images(n).gcp(:,1:2)', images(n).gcp(:,3:5)', images(n)))';
        % calculate mean error
        err(n) = mean(sqrt(sum(images(n).gcperr.^2, 2)));    
    end

    % Return mean per-image reprojection error
    meanerr = mean(err);
end

end