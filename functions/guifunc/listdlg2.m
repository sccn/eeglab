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
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [vals, okornot, strval] = listdlg2(varargin);

if nargin < 2
   help listdlg2;
   return;
end
for index = 1:length(varargin)
	if iscell(varargin{index}), varargin{index} = { varargin{index} }; end
	if ischar(varargin{index}), varargin{index} = lower(varargin{index}); end
end
g = struct(varargin{:});

try,  g.promptstring;  catch, g.promptstring = ''; end
try,  g.liststring;    catch, error('''liststring'' must be defined'); end
try,  g.selectionmode; catch, g.selectionmode = 'multiple'; end
try,  g.listsize;      catch, g.listsize = []; end
try,  g.initialvalue;  catch, g.initialvalue = []; end
try,  g.name;          catch, g.name = ''; end

fig = figure('visible', 'off');
set(gcf, 'name', g.name);
if ischar(g.liststring)
    allstr =  g.liststring;
else
    allstr = '';
    for index = 1:length(g.liststring)
        allstr = [ allstr '|' g.liststring{index} ];
    end
    allstr = allstr(2:end);
end

geometry = {[1] [1 1]};
geomvert = [min(length(g.liststring), 10) 1];
if ~strcmpi(g.selectionmode, 'multiple') || ...
        (iscell(g.liststring) && length(g.liststring) == 1) || ...
        (isstr (g.liststring) && size  (g.liststring,1) == 1 && isempty(find(g.liststring == '|')))
	if isempty(g.initialvalue), g.initialvalue = 1; end
    minval = 1;
	maxval = 1;
else
    minval = 0;
    maxval = 2;
end
listui = {{ 'Style', 'listbox', 'tag', 'listboxvals', 'string', allstr, 'max', maxval, 'min', minval } ...
		  { 'Style', 'pushbutton', 'string', 'Cancel', 'callback', ['set(gcbf, ''userdata'', ''cancel'');'] }  ...
		  { 'Style', 'pushbutton', 'string', 'Ok'    , 'callback', ['set(gcbf, ''userdata'', ''ok'');'] } };

if ~isempty(g.promptstring)
	geometry = {[1] geometry{:}};
	geomvert = [2 geomvert];
	listui = { { 'Style', 'text', 'string', g.promptstring } listui{:}};
end
[tmp, tmp2, allobj] = supergui( fig, geometry, geomvert, listui{:} );

% assign value to listbox
% must be done after creating it 
% ------------------------------
lstbox = findobj(fig, 'tag', 'listboxvals');
set(lstbox, 'value', g.initialvalue);

if ~isempty(g.listsize)
	pos = get(gcf, 'position');
	set(gcf, 'position', [ pos(1:2) g.listsize]);
end
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
        end
    else
        for index = vals
            strval = [ strval ' ' g.liststring(index,:) ];
        end
    end;        
	strval = strval(2:end);
	if strcmp(get(fig, 'userdata'), 'cancel')
		okornot = 0;
	else
		okornot = 1;
	end
	close(fig);
    drawnow;
end
