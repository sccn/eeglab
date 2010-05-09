% readelp() - read electrode locations from an .elp (electrode positions)
%             file as generated, for example, by a Polhemus tracking device 
% Usage:
%   >> [eloc, elocnames, X, Y, Z] = readelp(filename);
%
% Inputs:
%   filename - name of the .elp file containing cartesian (XYZ) electrode
%              locations
% Outputs:
%   eloc      - structure containing the names and locations of the channels
%   elocnames - cell array containing the names of the electrodes
%   X,Y,Z     - vectors containing the xyz locations of the electrodes
%
% Author: Arnaud Delorme, CNL / Salk Institute, 28 Feb 2002
%
% Note: ignores comments and the sensor type field
% Note: convert output XYZ locations to polar coordinates using cart2pol()
%
% See also: readlocs(), snapread(), floatread(), cart2pol()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 28 Feb 2002
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

function [eloc, names, x, y, z] = readelp( filename ); 

if nargin < 1
	help readelp;
	return;
end;

% open file
% ---------
fid = fopen(filename, 'r');
if fid == -1
  disp('Cannot open file'); return;
end;
for index=1:8	fgetl(fid); end;

% scan file
% ---------
index = 1;
tmpstr = fgetl(fid);
noteof = 1;
countfid  = 1;
fidlabels = { 'Nz' 'LPA' 'RPA' };
while noteof
    if ~isempty(deblank(tmpstr))
        if ~((tmpstr(1) == '/') & (tmpstr(2) == '/'))
            
            if tmpstr(1) == '%'
                if tmpstr(2) == 'F' % fiducial
                    tmp = sscanf(tmpstr(3:end), '%f');
                    
                    eloc(index).labels = fidlabels{countfid};
                    eloc(index).X  = tmp(1); x(index) = tmp(1);
                    eloc(index).Y  = tmp(2); y(index) = tmp(2);
                    eloc(index).Z  = tmp(3); z(index) = tmp(3);
                    eloc(index).type = 'FID';
                    index     = index    + 1;
                    countfid  = countfid + 1;
                    
                elseif tmpstr(2) == 'N' % regular channel
                    
                    eloc(index).labels =  strtok( tmpstr(3:end) );
                    tmpstr = fgetl(fid);
                    tmp = sscanf(tmpstr, '%f');
                    if isempty(tmpstr)
                        tmpstr = fgetl(fid);
                        tmp = sscanf(tmpstr, '%f');
                    end;
                    
                    eloc(index).X  = tmp(1); x(index) = tmp(1);
                    eloc(index).Y  = tmp(2); y(index) = tmp(2);
                    eloc(index).Z  = tmp(3); z(index) = tmp(3);
                    eloc(index).type = 'EEG';
                    index = index + 1;
                end;
            end;
        end;
    end;
    tmpstr = fgetl(fid);
    if ~isstr(tmpstr) & tmpstr == -1
        noteof = 0;
    end;  
end;  

names = { eloc.labels };
