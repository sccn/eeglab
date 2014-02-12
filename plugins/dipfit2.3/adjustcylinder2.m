% adjustcylinder() - Adjust 3d object coordinates to match a pair of points
%
% Usage:
%   >> [x y z] = adjustcylinder( x, y, z, pos1, pos2);
%
% Inputs:
%  x,y,z      - 3-D point coordinates
%  pos1       - position of first point [x y z]
%  pos2       - position of second point [x y z]
%
% Outputs:
%  x,y,z      - updated 3-D point coordinates
%
% Author: Arnaud Delorme, CNL / Salk Institute, 30 Mai 2003

% Copyright (C) 2003 Arnaud Delorme
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

function [x, y, z] = adjustcylinder2( h, pos1, pos2);
    
    % figure; plot3(x(2,:),y(2,:),z(2,:)); [ x(2,:)' y(2,:)' z(2,:)']
    
    % stretch z coordinates to match for vector length
    % ------------------------------------------------
    dist = sqrt(sum((pos1-pos2).^2));
    z = get(h, 'zdata');
    zrange = max(z(:)) - min(z(:)); 
    set(h, 'zdata', get(h, 'zdata') /zrange*dist);
    
    % rotate in 3-D to match vector angle [0 0 1] -> vector angle)
    % only have to rotate in the x-z and y-z plane
    % --------------------------------------------
    vectrot = [ pos2(1)-pos1(1) pos2(2)-pos1(2) pos2(3)-pos1(3)];
    [thvect phivect] = cart2sph( vectrot(1), vectrot(2), vectrot(3) ); 
    
    rotate(h, [0 0 1], thvect/pi*180, [0 0 0]);
    rotate(h, [thvect+pi/2 0]/pi*180, (pi/2-phivect)/pi*180, [0 0 0]);    

    x = get(h, 'xdata') + pos1(1);
    y = get(h, 'ydata') + pos1(2);
    z = get(h, 'zdata') + pos1(3);
    
    set(h, 'xdata', x);
    set(h, 'ydata', y);
    set(h, 'zdata', z);
    return;
    
