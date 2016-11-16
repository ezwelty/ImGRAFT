% UNDISTORT  Remove radial and tangential lens distortion from normalized
% image coordinates.
%
%   X = undistort(Xd, k)
%
% Uses the algorithm presented in "A four-step camera calibration procedure
% with implicit image correction" (Heikkila and Silven 1997).
%
% Input:    Xd          2xN distorted normalized image coordinates
%           k           distortion parameters [k1 k2 p1 p2 k3]
%           
% Output:   X           2XN undistorted normalized image coordinates
%
% See also image2world, image2camera.

% Adapted from code in the MATLAB Calibration Toolbox 
% (http://www.vision.caltech.edu/bouguetj/calib_doc/).

function X = undistort(Xd, k)

% If no distortion, nothing to do
if sum(abs(k)) == 0
    
    X = Xd;

% If only k1 is nonzero, use closed form solution
% Cubic roots solution from Numerical Recipes in C 2nd Ed. pages 184-185.
% (c) Jean-Yves Bouguet - April 27th, 1998
elseif sum(k ~= 0) == 1 && k(1) ~= 0
    
    ph = atan2(Xd(2,:), Xd(1,:));
    Q = -1 / (3 * k(1));
    R = -Xd(1,:) ./ (2 * k(1) * cos(ph));

    % Works in all practical situations (fails for very large k1)
    if k(1) < 0
        if k(1) < -0.2
            warning('Large, negative k1. Undistort may fail.')
        end
        th = acos(R ./ sqrt(Q^3));
        r = -2 * sqrt(Q) * cos((th - 2 * pi) / 3);
    % Always works, even for very large k1
    else
        A = (sqrt(R.^2 - Q^3) - R).^(1/3); 
        B = Q * (1./A);
        r = (A + B);
    end

    X = [r .* cos(ph) ; r .* sin(ph)];
    X = real(X);

% Otherwise, use iterative solution
else
    
    % Iterate
    X = Xd; % initial guess
    for n = 1:20
        % radial lens distortion
        % r = sqrt( x^2 + y^2 )
        % dr = (1 + k1 * r^2 + k2 * r^4 + k5 * r^6)
        r = sqrt(sum(X.^2, 1));
        dr = 1 + ones(2,1) * (k(1)*r.^2 + k(2)*r.^4 + k(5)*r.^6);

        % tangential lens distortion
        % dx = 2xy * k3 + k4 * (r^2 + 2x^2)
        % dy = k3 * (r^2 + 2y^2) + 2xy * k4
        dtx = 2*X(1,:).*X(2,:)*k(3) + k(4)*(r.^2 + 2*X(1,:).^2);
        dty = k(3)*(r.^2 + 2*X(2,:).^2) + 2*X(1,:).*X(2,:)*k(4);

        % undo distortion
        X = (Xd - [dtx;dty]) ./ dr;         
    end
end