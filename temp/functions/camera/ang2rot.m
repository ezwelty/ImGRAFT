% ANG2ROT  Calculate a rotation matrix from right-handed angles of rotation. 
%
%   R = ang2rot(angles)
%
% Computes rotation matrix (R) from angles omega (w), phi (p), kappa (k),
% which correspond to right-handed rotations (in degrees) about the x, y
% and z axes of a 3D euclidean coordinate system. By convention, the
% rotations are applied in the order w - p - k:
%
%   R = Rk * Rp * Rw
%
%
% Input:    angles      the x,y,z angles of rotation (in degrees) [w p k]
%
% Output:   R           3x3 rotation matrix
%
% See also rot2ang, rot2oblang, oblang2rot.

function R = ang2rot(angles)

% Check number of arguments
if (length(angles) == 3)
    k = angles(3); p = angles(2); w = angles(1);
else
    error('Specify three angles omega, phi and kappa (in degrees).')
end

% Prepare matrix element values
% (use of sind,cosd fixes floating point error of sin(pi), cos(pi/2))
sw = sind(w); cw = cosd(w);
sp = sind(p); cp = cosd(p);
sk = sind(k); ck = cosd(k);

% Build individual rotation matrices
rw = [ 1   0  0;
       0  cw sw;
       0 -sw cw ];
    
rp = [ cp  0 -sp;
        0  1   0;
       sp  0  cp ];
   
rk = [ ck  sk  0;
      -sk  ck  0;
        0   0  1 ];

% Concatenate rotation matrices
R = rk*rp*rw;
