% pluguin_uiupdate() - support function for plugin_menu

% Copyright (C) 2019 Arnaud Delorme
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

function plugin_uiupdate(fig)

tmpobj = get(fig, 'userdata');
listboxVal = get(findobj(fig, 'tag', 'pluginlist'), 'value');

allPlugins = tmpobj.allplugins;
selectedPlugins = tmpobj.selectedplugins;
selection       = tmpobj.selection;
if length(listboxVal) == 1
    lastSelect = listboxVal;
else
    lastSelect = setdiff(listboxVal, selection);
    
    % check if selection is consistent
    if ~all([selectedPlugins([ listboxVal ]).installed]) && ~all(~[selectedPlugins([ listboxVal ]).installed] )
        listboxVal = lastSelect(1);
        set(findobj(fig, 'tag', 'pluginlist'), 'value', listboxVal);
    end
end
if isempty(lastSelect), lastSelect = listboxVal(1); end
if length(lastSelect) > 1, lastSelect = lastSelect(1); end

tmpstr = selectedPlugins(lastSelect).description; 

% set text
hDescription = findobj(fig, 'tag', 'description');
hTags        = findobj(fig, 'tag', 'tags');
hStatus      = findobj(fig, 'tag', 'status');
tmpstr = textwrap(hDescription, {tmpstr}); 
set(hDescription, 'string', tmpstr);
set(hTags       , 'string', selectedPlugins(lastSelect).rawtags );

% installed status
if selectedPlugins(lastSelect).installed
    if selectedPlugins(lastSelect).installorupdate
        set(hStatus, 'string', ['version ' selectedPlugins(lastSelect).currentversion ' installed - an update is available'] );
    else
        set(hStatus , 'string', 'installed' );
    end
else
    set(hStatus, 'string', 'not installed');
end

% set buttons
set(findobj(fig, 'tag', 'install'), 'enable', 'on');
if selectedPlugins(lastSelect).installed
    set(findobj(fig, 'tag', 'remove'), 'enable', 'on');
else
    set(findobj(fig, 'tag', 'remove'), 'enable', 'off');
end
if ~isempty(selectedPlugins(lastSelect).webdoc)
    set(findobj(fig, 'tag', 'documentation'), 'enable', 'on', 'callback', ['web(''' selectedPlugins(lastSelect).webdoc ''', ''-browser'');' ]);
else
    set(findobj(fig, 'tag', 'documentation'), 'enable', 'off')
end
set(findobj(fig, 'tag', 'rating'),  'callback', ['web(''' selectedPlugins(lastSelect).webrating ''', ''-browser'');' ]);

tmpobj.selection = listboxVal;
set(fig, 'userdata', tmpobj);