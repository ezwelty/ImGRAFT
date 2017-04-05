% REPROJERR  Calculate the reprojection errors of a set of points.
%
%   err = reprojerr(Xi, X, camera)
%
% Returns the rerojection errors (in pixels) of a set of points given their
% world coordinates and nominal image coordinates.
%
% Input:    Xi          2XN image coordinates
%           X           3xN world coordinates
%           camera      camera structure                  
%           
% Output:   err         2xN pixel reprojection errors
%           
% See also mingcperr.

function err = reprojerr(Xi, X, camera)

% Project world coordinates onto the image
p = world2image(X, camera);

% Return reprojection errors
err = p - Xi;