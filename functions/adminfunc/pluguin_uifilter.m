function pluguin_uifilter(fig)

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