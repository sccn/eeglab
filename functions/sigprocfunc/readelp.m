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

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.6  2003/11/27 02:44:15  arno
% header typo
%
% Revision 1.5  2002/12/27 23:25:45  scott
% edit header msg -sm
%
% Revision 1.4  2002/12/27 22:48:43  arno
% fg
% name -> labels
%
% Revision 1.3  2002/11/15 01:39:58  scott
% Can not -> cannot
%
% Revision 1.2  2002/04/21 23:11:17  scott
% edited help msg -sm
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

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
for index=1:12	fgetl(fid); end;

% scan file
% ---------
index = 1;
tmpstr = fgetl(fid);
noteof = 1;
while noteof
    if ~isempty(deblank(tmpstr))
        if ~((tmpstr(1) == '/') & (tmpstr(2) == '/'))
            if (tmpstr(1) == '%') & (tmpstr(2) == 'N')
                eloc(index).labels =  strtok( tmpstr(3:end) );
                tmpstr = fgetl(fid);
                tmp = sscanf(tmpstr, '%f');
                eloc(index).X  = tmp(1); x(index) = tmp(1);
                eloc(index).Y  = tmp(2); y(index) = tmp(2);
                eloc(index).Z  = tmp(3); z(index) = tmp(3);
                index = index + 1;
            end;
        end;
    end;
    tmpstr = fgetl(fid);
    if ~isstr(tmpstr) & tmpstr == -1
        noteof = 0;
    end;  
end;  

names = { eloc.labels };
