function f = homogenous2traditional(H)

% HOMOGENOUS2TRADITIONAL estimates the traditional translation, rotation
% and scaling parameters from a homogenous transformation matrix. It will
% give an error if the homogenous matrix also describes a perspective
% transformation.
%
% Use as
%   f = homogenous2traditional(H)
% where H is a 4x4 homogenous transformation matrix and f is a vector with
% nine elements describing
%   x-shift
%	  y-shift
%	  z-shift
% followed by the
%	  pitch (rotation around x-axis)
%	  roll  (rotation around y-axis)
%	  yaw   (rotation around z-axis)
% followed by the
%   x-rescaling factor
%   y-rescaling factor
%	  z-rescaling factor
%
% The order in which the transformations would be done is exactly opposite
% as the list above, i.e. first z-rescale ... and finally x-shift.

% Copyright (C) 2005, Robert Oostenveld
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% remember the input homogenous transformation matrix
Horg = H;

% The homogenous transformation matrix is built up according to
%   H = T * R * S
% where
%   R = Rx * Ry * Rz

% estimate the translation
tx = H(1,4);
ty = H(2,4);
tz = H(3,4);
T = [
  1 0 0 tx
  0 1 0 ty
  0 0 1 tz
  0 0 0 1
  ];
% recompute the homogenous matrix excluding the translation
H = inv(T) * H;

% estimate the scaling
sx = norm(H(1:3,1));
sy = norm(H(1:3,2));
sz = norm(H(1:3,3));
S = [
  sx 0 0 0
  0 sy 0 0
  0 0 sz 0
  0 0 0  1
  ];
% recompute the homogenous matrix excluding the scaling
H = H * inv(S);

% the difficult part is to determine the rotations
% the order of the rotations matters

% compute the rotation using a probe point on the z-axis
p = H * [0 0 1 0]';
% the rotation around the y-axis is resulting in an offset in the positive x-direction
ry = asin(p(1));

% the rotation around the x-axis can be estimated by the projection on the yz-plane
if abs(p(2))<eps && abs(p(2))<eps
  % the rotation around y was pi/2 or -pi/2, therefore I cannot estimate the rotation around x any more
  error('need another estimate, not implemented yet');
elseif abs(p(3))<eps
  % this is an unstable situation for using atan, but the rotation around x is either pi/2 or -pi/2
  if p(2)<0
    rx = pi/2
  else
    rx = -pi/2;
  end
else
  % this is the default equation for determining the rotation
  rx = -atan(p(2)/p(3));
end

% recompute the individual rotation matrices
Rx = rotate([rx 0 0]);
Ry = rotate([0 ry 0]);
Rz = inv(Ry) * inv(Rx) * H; % use left side multiplication

% compute the remaining rotation using a probe point on the x-axis
p = Rz * [1 0 0 0]';
rz = asin(p(2));

% the complete rotation matrix was
R = rotate([rx ry rz]);

% compare the original translation with the one that was estimated
H = T * R * S;
%fprintf('remaining difference\n');
%disp(Horg - H);

f = [tx ty tz rx ry rz sx sy sz];


function [output] = rotate(R, input);

% ROTATE performs a 3D rotation on the input coordinates
% around the x, y and z-axis. The direction of the rotation 
% is according to the right-hand rule. The rotation is first
% done around the x-, then the y- and finally the z-axis.
% 
% Use as
%   [output] = rotate(R, input)
% where
%   R		[rx, ry, rz] rotations around each of the axes in degrees
%   input	Nx3 matrix with the points before rotation
%   output	Nx3 matrix with the points after rotation
%
% Or as
%   [Tr] = rotate(R)
% where
%   R		[rx, ry, rz] in degrees
%   Tr		corresponding homogenous transformation matrix

% Copyright (C) 2000-2004, Robert Oostenveld
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

rotx = eye(3);
roty = eye(3);
rotz = eye(3);

rx = pi*R(1) / 180;
ry = pi*R(2) / 180;
rz = pi*R(3) / 180;

if rx~=0
  % rotation around x-axis
  rotx(2,:) = [   0    cos(rx)  -sin(rx) ];
  rotx(3,:) = [   0    sin(rx)   cos(rx) ];
end

if ry~=0
  % rotation around y-axis
  roty(1,:) = [  cos(ry)   0     sin(ry) ];
  roty(3,:) = [ -sin(ry)   0     cos(ry) ];
end

if rz~=0
  % rotation around z-axis
  rotz(1,:) = [  cos(rz)  -sin(rz)   0 ];
  rotz(2,:) = [  sin(rz)   cos(rz)   0 ];
end

if nargin==1
  % compute and return the homogenous transformation matrix
  rotx(4,4) = 1;
  roty(4,4) = 1;
  rotz(4,4) = 1;
  output = rotz * roty * rotx;
else
  % apply the transformation on the input points
  output = ((rotz * roty * rotx) * input')';
end

