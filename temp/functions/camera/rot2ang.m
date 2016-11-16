% ROT2ANG  Decompose a rotation matrix into right-handed angles of rotation. 
%
%   angles = rot2ang(R)
%
% Computes angles omega (w), phi (p), kappa (k), corresponding to right-
% handed rotations (in degrees) about the x, y and z axes of a 3D euclidean
% coordinate system, from a rotation matrix (R). By convention, the
% rotations are applied in the order w - p - k:
%
%   R = Rk * Rp * Rw
%
%
% Input:    R           3x3 rotation matrix
%
% Output:   angles      the x,y,z angles of rotation (in degrees) [w p k]
%
% See also ang2rot, rot2oblang, oblang2rot.

function angles = rot2ang(R)

% Check number of arguments
if (nargin == 0)
    error('Specify a 3x3 rotation matrix.')
end

% Test size of R
[nr, nc] = size(R);
if (nr ~= 3 || nc ~= 3)
    error('Rotation matrix must be 3x3.')
end

% Calculate w,p,k (in degrees)
w = atand(-R(3,2) / R(3,3));
p = asind(R(3,1));
k = atand(-R(2,1) / R(1,1));
angles = [w p k];