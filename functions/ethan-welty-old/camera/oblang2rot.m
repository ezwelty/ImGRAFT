% OBLANG2ROT  Calculate a rotation matrix from heading, elevation, and bank
% angles.
%
%   R = oblang2rot(angles)
%
% Computes rotation matrix (R) from heading, elevation and bank angles by
% expressing them in turns of corresponding right-handed axis rotations
% about the z, x and y axes. By convention, the rotations are applied in
% the order
%
%   heading - elevation - bank,
%
% starting from a camera pointing in the +y axis (true north) flush with
% the horizon.
%
% Since the convention for R in use here is to be relative to a camera
% pointing in the +z axis, an initial -90 degree roll about the x-axis is
% needed. The final rotation matrix is given by
%
%   R = Rb * Re * Rh * Ri, or equivalently
%   R = Rk * Rw * Rp * Ri
%
% in terms of w,p,k, right-handed rotations about x,y,z (see ang2rot.m).
%
%
% Input:    angles [heading elevation bank]:
%           heading     degrees from y-axis (+CW, -CCW)
%           elevation   degrees from xy-plane (+lens up, -lens down)
%           bank        degrees about optical axis (+CW, -CCW, from behind)
%
% Output:   R           3x3 rotation matrix
%
% See also rot2oblang, rot2ang, ang2rot.

function R = oblang2rot(angles)

% Check number of arguments
if (length(angles) == 3)
    bank = angles(3); elevation = angles(2); heading = angles(1);
else
    error('Specify three angles heading, elevation, and bank (in degrees).')
end

% Initial rotation
% (+z to +y camera pointing direction, -90 deg roll about x-axis)
% (+z now pointing north, +x east, +y down)
ri =  [ 1   0  0;
        0   0 -1;
        0   1  0 ];

% Relate to w,p,k, right-handed rotations about x,y,z (see ang2rot.m)
p = heading; % (y-axis rotation)
w = elevation; % (x-axis rotation)
k = bank; % (z-axis rotation)

% Prepare matrix element values
% note: use of sind,cosd fixes floating point error of sin(pi), cos(pi/2)
sp = sind(p); cp = cosd(p);
sw = sind(w); cw = cosd(w);
sk = sind(k); ck = cosd(k);

% Build individual rotation matrices
rp = [ cp  0 -sp;
    0  1   0;
   sp  0  cp ];

rw = [ 1   0  0;
       0  cw sw;
       0  -sw cw ];

rk = [ ck   sk  0;
      -sk   ck  0;
        0    0  1 ];

% Concatenate rotation matrices
R = rk*rw*rp*ri;
