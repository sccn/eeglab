% topo2sph() - convert a topoplot() style 2-D polar-coordinates 
%              channel locations file to a 3-D spherical-angle
%              file for use with headplot()
%
% Usage: 
%   >> [c h] = topo2sph('eloc_file','eloc_angles_file');
%   >> [c h] = topo2sph( topoarray );
%
% Inputs:
%   'eloc_file' = filename of polar 2-d electrode locations file used by topoplot()
%                 See >> topoplot example or cart2topo()
%   'eloc_angles_file' = output file of electrode locations in spherical angle coords.
%                        for use in headplot().
%   topoarray = polar array of 2-d electrode locations, with polar angle in the
%               first column and radius in the second one.
%
% Outputs:
%   c = coronal rotation
%   h = horizontal rotation
%
% Author: Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 1999 
%
% See also: sph2topo(), cart2topo()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 1999 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 3-16-00 changed name to topo2sph() for compatibility with cart2topo() -sm
% 01-25-02 reformated help & license -ad 
% 03-22-02 complete remodeling for returning arguments and taking arrays -ad 

function [c, h] = topo2sph(eloc_locs,eloc_angles)

MAXCHANS = 1024;

if nargin < 1
    help topo2sph;
    return;
end;
    
if isstr(eloc_locs)
	fid = fopen(eloc_locs);
	if fid<1,
	    fprintf('topo2sph()^G: cannot open eloc_loc file (%s)\n',eloc_locs)
	    return
	end
	E = fscanf(fid,'%d %f %f  %s',[7 MAXCHANS]);
	E = E';
	fclose(fid);
else
    E = eloc_locs;
    E = [ ones(size(E,1),1) E ];
end;
    
if nargin > 1
	if exist(eloc_angles)==2,
	   fprintf('plo2sph(): eloc_angles file (%s) already exists.\n',eloc_angles);
	   return
	end

	fid = fopen(eloc_angles,'a');
	if fid<1,
	    fprintf('topo2sph()^G: cannot open eloc_angles file (%s)\n',eloc_angles)
	    return
	end
end;

t = E(:,2); % theta
r = E(:,3); % radius
h = -t;  % horizontal rotation
c = (0.5-r)*180;

for e=1:size(E,1)
   % (t,r) -> (c,h)

   %if t>=0
   %   h(e) = 90-t; % horizontal rotation
   %else
   %   h(e) = -(90+t);
   %end
   %if t~=0
   %   c(e) = sign(t)*180*r; % coronal rotation
   %else
   %   c(e) = 180*r;
   %end

   if nargin > 1
        chan = E(e,4:7);
        fprintf('%d	%g	%g	%s\n',E(e,1),c(e),h(e),chan);
        fprintf(fid,'%d	%g	%g	%s\n',E(e,1),c(e),h(e),chan);
   end;     
end

