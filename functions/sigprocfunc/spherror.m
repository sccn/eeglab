% spherror() - chancenter() sub function to compute minimum distance
%               of Cartesian coordinates to a sphere
%
% Author: Scott Makeig, CNL / Salk Institute, 2000

% Copyright (C) Scott Makeig, CNL / Salk Institute, 2000
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

function wobble = spherror(center,x,y,z)
  
% 03/14/02 corrected wobble calculation -lf
  
x = x - center(1);  % center the points at (0,0,0)
y = y - center(2);
z = z - center(3);
radius = (sqrt(x.^2+y.^2+z.^2)); % distances from the center
wobble = std(radius-mean(radius)); % test if xyz values are on a sphere
