% convertelocs() - Convert between frames for electrode locations
%
% Usage: >> newchans = convertelocs( EEG, 'command');
%
% Input:
%   chanlocs  - EEGLAB dataset or EEG.chanlocs structure
%   'command' - Can be any of the following
%               ['cart2topo'|'sph2topo'|'sphbesa2topo'|
%               'sph2cart'|'topo2cart'|'sphbesa2cart'|
%               'cart2sph'|'sphbesa2sph'|'topo2sph'|
%               'cart2sphbesa'|'sph2sphbesa'|'topo2sphbesa'|
%               'cart2all'|'sph2all'|'sphbesa2all'|'topo2all']
%               to convert between the 4 coordinate frames
%               3-D Carthesian (cart), Matlab spherical (sph)
%               Besa spherical (sphbesa) and 2-D polar (topo)
%               The command can also be 'auto'. Then the function
%               find the most complex coordinate frame and constrain
%               all the others to this one. By order, it first searches
% 					 for carthesian ccordinates, then spherical and finaly
%               polar. Default is 'auto'.
%
% Outputs:
%   newchans      - new EEGLAB channel locations structure
%
% Ex:  CHANSTRUCT = convertelocs( CHANSTRUCT, 'cart2topo'
%        % Convert Cathesian coordinates to topographical ones
%
% Author: Arnaud Delorme, CNL / Salk Institute, 22 Dec 2002
%
% See also: readlocs()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 22 Dec 2002, arno@salk.edu
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

function chans = convertelocs(chans, command);

if nargin < 1
   help convertelocs;
   return;
end;

if nargin < 2
   command = 'auto';
end;

% test if value exists for default
% --------------------------------
if strcmp(command, 'auto')
   if isfield(chans, 'X') & ~isempty(chans(1).X)
      command = 'cart2all';
      disp('Uniformize all coordinate frames using Carthesian coords');
   else
	   if isfield(chans, 'sph_theta') & ~isempty(chans(1).sph_theta)
   	   command = 'sph2all';
	      disp('Uniformize all coordinate frames using spherical coords');
      else
		   if isfield(chans, 'sph_thetabesa') & ~isempty(chans(1).sph_thetabesa)
   		   command = 'sphbesa2all';
      		disp('Uniformize all coordinate frames using BESA spherical coords');
         else
            command = 'topo2all';
		      disp('Uniformize all coordinate frames using polar coords');
         end;
      end;
   end;
end;

% convert
% -------         
switch command
case 'topo2sph',
   [sph_phi sph_theta] = topo2sph( [cell2mat({chans.theta})' cell2mat({chans.radius})'] );
   for index = 1:length(chans)
      chans(index).sph_theta  = sph_theta(index);
      chans(index).sph_phi    = sph_phi  (index);
      chans(index).sph_radius = 1;
   end;
case 'topo2sphbesa',
   chans = convertelocs(chans, 'topo2sph'); % search for spherical coords
   chans = convertelocs(chans, 'sph2sphbesa'); % search for spherical coords
case 'topo2cart'
   chans = convertelocs(chans, 'topo2sph'); % search for spherical coords
   disp('Warning: spherical coordinates automatically updated');
   chans = convertelocs(chans, 'sph2cart'); % search for spherical coords
case 'topo2all',
   chans = convertelocs(chans, 'topo2sph'); % search for spherical coords
   chans = convertelocs(chans, 'sph2sphbesa'); % search for spherical coords
   chans = convertelocs(chans, 'sph2cart'); % search for spherical coords
case 'sph2cart',
   [x y z] = sph2cart(cell2mat({chans.sph_theta})'/180*pi, cell2mat({chans.sph_phi})'/180*pi, cell2mat({chans.sph_radius})');
   for index = 1:length(chans)
      chans(index).X = x(index);
      chans(index).Y = y(index);
      chans(index).Z = z(index);
   end;
case 'sph2topo',
   disp('Warning: all radii considered to one for this transformation');
   try, [chan_num,angle,radius] = sph2topo([ones(length(chans),1)  cell2mat({chans.sph_phi})' cell2mat({chans.sph_theta})'], 1, 2); % using method 2
   catch, error('Cannot process empty values'); end;
   for index = 1:length(chans)
      chans(index).theta  = angle(index);
      chans(index).radius = radius(index);
   end;
case 'sph2sphbesa',
   % using polar coordinates
   [chan_num,angle,radius] = sph2topo([ones(length(chans),1)  cell2mat({chans.sph_phi})' cell2mat({chans.sph_theta})'], 1, 2);
   [sph_thetabesa sph_phibesa] = topo2sph([angle radius], 1, 1);
   for index = 1:length(chans)
      chans(index).sph_thetabesa  = sph_thetabesa(index);
      chans(index).sph_phibesa    = sph_phibesa  (index);
   end;   
case 'sph2all',
   chans = convertelocs(chans, 'sph2topo'); % search for spherical coords
   chans = convertelocs(chans, 'sph2sphbesa'); % search for spherical coords
   chans = convertelocs(chans, 'sph2cart'); % search for spherical coords
case 'sphbesa2sph',
   % using polar coordinates
   [chan_num,angle,radius] = sph2topo([ones(length(chans),1)  cell2mat({chans.sph_phi})' cell2mat({chans.sph_theta})'], 1, 1);
   [sph_theta sph_theta] = topo2sph([angle radius], 2);
   for index = 1:length(chans)
      chans(index).sph_theta  = sph_theta(index);
      chans(index).sph_phi    = sph_phi  (index);
   end;
case 'sphbesa2topo',
   chans = convertelocs(chans, 'sphbesa2sph'); % search for spherical coords
   chans = convertelocs(chans, 'sph2topo'); % search for spherical coords
case 'sphbesa2cart',
   chans = convertelocs(chans, 'sphbesa2sph'); % search for spherical coords
   chans = convertelocs(chans, 'sph2cart'); % search for spherical coords   
case 'sphbesa2all',
   chans = convertelocs(chans, 'sphbesa2sph'); % search for spherical coords
   chans = convertelocs(chans, 'sph2all'); % search for spherical coords
case 'cart2topo',
   chans = convertelocs(chans, 'cart2sph'); % search for spherical coords
   chans = convertelocs(chans, 'sph2topo'); % search for spherical coords
case 'cart2sphbesa',
   chans = convertelocs(chans, 'cart2sph'); % search for spherical coords
   chans = convertelocs(chans, 'sph2sphbesa'); % search for spherical coords
case 'cart2sph',
	[th phi radius] = cart2sph(cell2mat({chans.X}), cell2mat({chans.Y}), cell2mat({chans.Z}));
	for index = 1:length(chans)
		 chans(index).sph_theta     = th(index)/pi*180;
		 chans(index).sph_phi       = phi(index)/pi*180;
		 chans(index).sph_radius    = radius(index);
	end;
case 'cart2all',
   chans = convertelocs(chans, 'cart2sph'); % search for spherical coords
   chans = convertelocs(chans, 'sph2all'); % search for spherical coords
end;