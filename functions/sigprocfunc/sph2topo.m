% sph2topo() - Convert from a 3-column headplot file in spherical coordinates
%             to a 3-column topoplot file in polar (not cylindrical) coords.
%             Used for topoplot() and other 2-d topographic plotting programs.
%             Assumes a spherical coordinate system in which horizontal angles 
%             has a range [-180,180] with zero pointing to the right ear. 
%             In the output polar coordinate system, zero points to the nose.
% Usage:
%          >> [chan_num,angle,radius] = sph2topo(input,shrink_factor,method);
%
% Inputs:
%   input         = [channo,az,horiz] = chan_number, azumith (deg), horiz. angle (deg)
%                   When az>0, horiz=0 -> right ear, 90 -> nose 
%                   When az<0, horiz=0 -> left ear, -90 -> nose
%   shrink_factor = radial scaling factor>=1 (Note: 1 -> plot edge 90 deg az
%                   1.5 -> plot edge is +/-135 deg az {default 1}
%   method        = [1|2], optional. 1 is standard, 2 is for
%                   compatibility with cart2sph.
%
% Outputs:
%   channo  = channel number (as in input)
%   angle   = horizontal angle (0 -> nose; 90 -> right ear; -90 -> left ear)
%   radius  = arc radius from vertex (Note: 90 deg az -> 0.5/shrink_factor);
%             By convention, radius=0.5 is the outer edge of topoplot().
%             Use shrink_factor>1 to plot chans with abs(az)>90.
%
% Author: Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 6/12/98 
%
% See also: cart2topo(), topo2sph()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 6/12/98 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% $Log: not supported by cvs2svn $

% corrected left/right orientation mismatch, Blair Hicks 6/20/98
% changed name sph2pol() -> sph2topo() for compatibility -sm
% 01-25-02 reformated help & license -ad 
% 01-25-02 changed computation so that it works with sph2topo -ad 

function [channo,angle,radius] = sph2topo(input,factor, method)

chans = size(input,1);
angle = zeros(chans,1);
radius = zeros(chans,1);

if nargin < 1
   help sph2topo
   return
end
   
if nargin< 2
  factor = 0;
end
if factor==0
  factor = 1;
end
if factor < 1
  help sph2topo
  return
end

if size(input,2) ~= 3
   help sph2topo
   return
end

channo = input(:,1);
az = input(:,2);
horiz = input(:,3);

if exist('method')== 1 & method == 2
  angle  = -horiz;
  radius = 0.5 - abs(az/180)/factor;
else
  radius = abs(az/180)/factor;
  i = find(az>=0);
  angle(i) = 90-horiz(i);
  i = find(az<0);
  angle(i) = -90-horiz(i);
end;
