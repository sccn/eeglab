% pop_saveh() - save the EEGLAB session command history stored in ALLCOM
%
% Usage:
%   >> pop_saveh( ALLCOM, filename, filepath);
%
% Inputs:
%   ALLCOM     - cell array of strings containing the EEGLAB command history 
%   filename   - name of the file to save to (optional, default "eeglabhist.m"
%   filepath   - path of the file to save to (optional, default pwd)
%
% Author: Arnaud Delorme, CNL / Salk Institute, 22 March 2002
%
% See also: h(), eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 22 March 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.3  2002/08/13 16:40:23  scott
% help msg and text edits
%
% Revision 1.2  2002/04/23 23:33:14  arno
% changing last call
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-29-02 added update menu -ad 

function com = pop_saveh( allcoms, curfilename, curfilepath);

com = '';
if nargin < 1
	help pop_saveh;
	return;
end;
	
if nargin < 3
	[curfilename, curfilepath] = uiputfile('eeglabhist.m', 'Save the EEGLAB session command history with .m extension -- pop_saveh()');
    drawnow;
	if curfilename == 0 return; end;
end;

disp('Saving the EEGLAB session command history...');

fid = fopen( [ curfilepath curfilename ], 'w');
if fid == -1
    error('pop_saveh(): Cannot open named file');
end;    
fprintf(fid, '%% EEGLAB history file generated on the %s\n', date);
fprintf(fid, '%% ------------------------------------------------\n');
for index = length(allcoms):-1:1
    fprintf(fid, '%s\n', allcoms{index});
end;
fprintf(fid, 'eeglab redraw;\n');
fclose(fid);
    
com = sprintf('pop_saveset( %s, ''%s'', ''%s'');', inputname(1), curfilename, curfilepath);

return;
