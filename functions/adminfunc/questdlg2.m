% questdlg2() - questdlg function clone with coloring and help for 
%               eeglab().
%
% Usage: same as questdlg()
%
% Warning: 
% Case of button text and result might be changed by the function
%         
%
% Author: Arnaud Delorme, CNL / Salk Institute, La Jolla, 11 August 2002
%
% See also: inputdlg2(), errordlg2(), supergui(), inputgui()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, arno@salk.edu
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
% Revision 1.9  2002/10/15 14:45:12  arno
% display warning
%
% Revision 1.8  2002/10/15 14:35:52  arno
% default case for buttons
%
% Revision 1.7  2002/08/14 18:17:00  arno
% new supergui arg
%
% Revision 1.6  2002/08/13 22:36:43  arno
% debug for error
%
% Revision 1.5  2002/08/13 17:29:35  arno
% new superguicall
%
% Revision 1.4  2002/08/12 18:49:20  arno
% header
%
% Revision 1.3  2002/08/12 18:24:29  arno
% debug
%
% Revision 1.2  2002/08/12 18:02:47  arno
% debug
%
% Revision 1.1  2002/08/12 18:01:34  arno
% Initial revision
%

function [result] = questdlg2(Prompt,Title,varargin);

result = varargin{end};
if nargin < 2
   help questdlg2;
   return;
end;
if Prompt(end) == 10, Prompt(end) = []; end;
if Prompt(end) == 10, Prompt(end) = []; end;
if Prompt(end) == 10, Prompt(end) = []; end;
if Prompt(end) == 10, Prompt(end) = []; end;

fig = figure('visible', 'off');
set(gcf, 'name', Title);

geometry = {[1]};
listui = {{ 'Style', 'text', 'string' Prompt }};

geometry = { geometry{:} ones(1,length(varargin)-1) };
for index = 1:length(varargin)-1 % ignoring default val
	listui = {listui{:} { 'Style', 'pushbutton', 'string', varargin{index}, 'callback', ['set(gcbf, ''userdata'', ''' varargin{index} ''');'] }  };
	if strcmp(varargin{index}, varargin{end})
		listui{end}{end+1} = 'fontweight';
		listui{end}{end+1} = 'bold';
	end;
end;

cr = length(find(Prompt == char(10)))+1;
if cr == 1
	cr = size(Prompt,1);
end;

[tmp tmp2 allobj] = supergui( fig, geometry, [cr 1], listui{:} );

waitfor( fig, 'userdata');
drawnow;
try,
	result = get(fig, 'userdata');
	close(fig);
end;