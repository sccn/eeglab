% plugin_uifilter() - support function for plugin_menu

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

function plugin_uifilter(fig)

tmpobj = get(fig, 'userdata');
allPlugins = tmpobj.allplugins;

filterInstall = get(findobj(fig, 'tag', 'filter1'), 'value');

filterVal = get(findobj(fig, 'tag', 'filter2'), 'value');
filterStr = get(findobj(fig, 'tag', 'filter2'), 'string');
filterStr = filterStr{filterVal};

if filterVal == 1
    selectedPlugins = allPlugins;
else
    spaceInd = find(filterStr == ' ');
    filterStr = filterStr(spaceInd(2)+1:end);
    
    % find plugins
    pluginList = [];
    for iPlugin = 1:length(allPlugins)
        if ~isempty(strmatch(filterStr, allPlugins(iPlugin).tags, 'exact'))
            pluginList = [ pluginList iPlugin ];
        end
    end
    selectedPlugins = allPlugins(pluginList);
end

if filterInstall == 2
    selectedPlugins( [selectedPlugins.installed] == 0) = [];
elseif filterInstall == 3
    selectedPlugins( [selectedPlugins.installed] == 1) = [];
end

% update GUI
set(findobj(fig, 'tag', 'pluginlist'), 'string', { selectedPlugins.text }, 'value', []);
tmpobj.selection = [];
tmpobj.selectedplugins = selectedPlugins;

set(fig, 'userdata', tmpobj);