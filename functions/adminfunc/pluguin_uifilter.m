function pluguin_uifilter(fig)

tmpobj = get(fig, 'userdata');
allPlugins = tmpobj.allplugins;

filterVal = get(findobj(fig, 'tag', 'filter'), 'value');
filterStr = get(findobj(fig, 'tag', 'filter'), 'string');
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

% update GUI
set(findobj(fig, 'tag', 'pluginlist'), 'string', { selectedPlugins.text }, 'value', []);
tmpobj.selection = [];
tmpobj.selectedplugins = selectedPlugins;

set(fig, 'userdata', tmpobj);