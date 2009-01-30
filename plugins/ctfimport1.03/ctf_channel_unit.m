function [sens_unit] = ctf_channel_unit(ctf)

% ctf_channel_unit - project CTF channel locations onto a unit sphere
% 
% [sens_unit] = ctf_channel_unit(ctf)
%
% sens_unit output is a projection of all the CTF sensor locations onto a
% unit sphere, assuming the origin is at (0,0,0).  The output vector is Nx3
% Cartesian coordinates.
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%      <                                                      >
%      <                    DISCLAIMER:                       >
%      <                                                      >
%      < THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. >
%      < THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   >
%      <                    OFFICIAL USE.                     >
%      <                                                      >
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>
%


% $Revision: 1.1 $ $Date: 2009-01-30 03:49:26 $

% Copyright (C) 2005  Darren L. Weber
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Modified: 01/2005, Darren.Weber_at_radiology.ucsf.edu
%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ver = '$Revision: 1.1 $ $Date: 2009-01-30 03:49:26 $';
fprintf('\nCTF_CHANNEL_UNIT [v %s]\n',ver(11:15));

% translate the CTF sensors into a unit sphere, 
% assuming the origin is at (0,0,0)
sens = ctf.sensor.location';
distance = sqrt( sens(:,1).^2 + sens(:,2).^2 + sens(:,3).^2 );
distance = repmat(distance,1,3);
sens_unit = sens ./ distance;

return
