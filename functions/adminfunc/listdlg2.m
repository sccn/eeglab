% listdlg2() - listdlg function clone with coloring and help for 
%               eeglab().
%
% Usage: same as listdlg()
%
% Author: Arnaud Delorme, CNL / Salk Institute, La Jolla, 16 August 2002
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
% Revision 1.16  2006/02/23 00:42:57  arno
% nothing
%
% Revision 1.15  2006/02/23 00:40:32  arno
% multiple selection
%
% Revision 1.14  2006/02/23 00:35:02  arno
% empty value case
%
% Revision 1.13  2006/02/23 00:33:52  arno
% adding value argument
%
% Revision 1.12  2006/02/08 23:15:23  arno
% allowing non-cell array input
%
% Revision 1.11  2006/01/10 00:40:28  arno
% fix multiline
%
% Revision 1.10  2004/11/10 17:06:39  arno
% initval -> initialvalue
%
% Revision 1.9  2004/11/10 16:51:25  arno
% remove dbug message
%
% Revision 1.8  2004/11/10 16:48:20  arno
% nothing
%
% Revision 1.7  2004/11/10 16:33:18  arno
% maximum height
%
% Revision 1.6  2004/06/28 15:44:43  arno
% redrawing figure
%
% Revision 1.5  2004/06/16 22:00:20  arno
% debug for integer
%
% Revision 1.4  2002/10/15 17:24:13  arno
% same
%
% Revision 1.3  2002/10/15 16:56:45  arno
% drawnow for windows
%
% Revision 1.2  2002/08/21 18:13:37  arno
% debug for one element only
%
% Revision 1.1  2002/08/17 03:06:48  arno
% Initial revision
%

function [vals, okornot, strval] = listdlg2(varargin);

if nargin < 2
   help listdlg2;
   return;
end;
for index = 1:length(varargin)
	if iscell(varargin{index}), varargin{index} = { varargin{index} }; end;
	if isstr(varargin{index}), varargin{index} = lower(varargin{index}); end;
end;
g = struct(varargin{:});

try,  g.promptstring;  catch, g.promptstring = ''; end;
try,  g.liststring;    catch, error('''liststring'' must be defined'); end;
try,  g.selectionmode; catch, g.selectionmode = 'multiple'; end;
try,  g.listsize;      catch, g.listsize = []; end;
try,  g.initialvalue;  catch, g.initialvalue = []; end;
try,  g.name;          catch, g.name = ''; end;
if isempty(g.value), g.value = 1; end;

fig = figure('visible', 'off');
set(gcf, 'name', g.name);
if isstr(g.liststring)
    allstr =  g.liststring;
else
    allstr = '';
    for index = 1:length(g.liststring)
        allstr = [ allstr '|' g.liststring{index} ];
    end;
    allstr = allstr(2:end);
end;

geometry = {[1] [1 1]};
geomvert = [min(length(g.liststring), 10) 1];
if ~strcmpi(g.selectionmode, 'multiple') | ...
        (iscell(g.liststring) & length(g.liststring) == 1) | ...
        (isstr (g.liststring) & size  (g.liststring,1) == 1 & isempty(find(g.liststring == '|')))
	if isempty(g.initialvalue), g.initialvalue = 1; end;
    minval = 1;
	maxval = 1;
else
    minval = 0;
    maxval = 2;
end;
listui = {{ 'Style', 'listbox', 'tag', 'listboxvals', 'string', allstr, 'max', maxval, 'min', minval } ...
		  { 'Style', 'pushbutton', 'string', 'Cancel', 'callback', ['set(gcbf, ''userdata'', ''cancel'');'] }  ...
		  { 'Style', 'pushbutton', 'string', 'Ok'    , 'callback', ['set(gcbf, ''userdata'', ''ok'');'] } };

if ~isempty(g.promptstring)
	geometry = {[1] geometry{:}};
	geomvert = [1 geomvert];
	listui = { { 'Style', 'text', 'string', g.promptstring } listui{:}};
end;
[tmp tmp2 allobj] = supergui( fig, geometry, geomvert, listui{:} );

% assign value to listbox
% must be done after creating it 
% ------------------------------
lstbox = findobj(fig, 'tag', 'listboxvals');
set(lstbox, 'value', g.initialvalue);

if ~isempty(g.listsize)
	pos = get(gcf, 'position');
	set(gcf, 'position', [ pos(1:2) g.listsize]);
end;
h = findobj( 'parent', fig, 'tag', 'listboxvals');
	
okornot = 0;
strval = '';
vals = [];
figure(fig);
drawnow;
waitfor( fig, 'userdata');
try,
	vals = get(h, 'value');
	strval = '';
    if iscell(g.liststring)
        for index = vals
            strval = [ strval ' ' g.liststring{index} ];
        end;
    else
        for index = vals
            strval = [ strval ' ' g.liststring(index,:) ];
        end;
    end;        
	strval = strval(2:end);
	if strcmp(get(fig, 'userdata'), 'cancel')
		okornot = 0;
	else
		okornot = 1;
	end;
	close(fig);
    drawnow;
end;