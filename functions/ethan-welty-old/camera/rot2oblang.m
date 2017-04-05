% ROT2OBLANG  Decompose a rotation matrix into heading, elevation, and bank
% angles.
%
%   angles = rot2oblang(R, verbose = false)
%
% Returns heading, elevation and bank angles for a rotation matrix (R). See
% oblang2rot.m for more angle definitions and the equivalence with the 
% w,p,k convention.
%
%
% Input:    R           3x3 rotation matrix
%           verbose     boolean whether to display warning messages
%                       [default=0]
%
% Output:   angles [heading elevation bank]:
%             heading     degrees from y-axis (+CW, -CCW)
%             elevation   degrees from xy-plane (+lens up, -lens down)
%             bank        degrees about optical axis (+CW, -CCW, from behind)
%           
% See also oblang2rot, rot2ang, ang2rot.

function angles = rot2oblang(R, verbose)

% Check number of arguments
if (nargin == 0)
    error('Specify a 3x3 rotation matrix.')
end

% Test size of R
[nr, nc] = size(R);
if (nr ~= 3 || nc ~= 3)
    error('Rotation matrix must be 3x3.')
end

% Enforce defaults
if (nargin < 2)
    verbose = false;
end

% Calculate angles (in degrees)
% cos(elevation) != 0
if (abs(R(3,3)) ~= 1)
    w1 = asin(R(3,3)); 
    p1 = atan2(R(3,1)/cos(w1), R(3,2)/cos(w1)); 
    k1 = atan2(-R(1,3)/cos(w1), -R(2,3)/cos(w1));
    if (verbose)
        w2 = pi - w1;
        p2 = atan2(R(3,1)/cos(w2), R(3,2)/cos(w2));
        k2 = atan2(-R(1,3)/cos(w2), -R(2,3)/cos(w2));
        disp('Two valid solutions exist. Only the first was returned.')
        disp(rad2deg([p1 w1 k1 ; p2 w2 k2]))
    end
    angles = rad2deg([p1 w1 k1]);

% cos(elevation) == 0
else 
    k = 0; % (unconstrained)
    if (R(3,3) == 1)
        w = pi/2;
        p = -k + atan2(-R(1,2), R(1,1));
    else
        w = -pi/2;
        p = k + atan2(-R(1,2), R(1,1));
    end
    
    if (verbose)
        disp('Infinite solutions exist. bank(z) was set to zero.')
    end
    angles = rad2deg([p w k]);
end