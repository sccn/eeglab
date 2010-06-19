% listdlg2() - listdlg function clone with coloring and help for 
%               eeglab().
%
% Usage: same as listdlg()
%
% Author: Arnaud Delorme, CNL / Salk Institute, La Jolla, 16 August 2002
%
% See also: inputdlg2(), errordlg2(), supergui(), inputgui()

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
