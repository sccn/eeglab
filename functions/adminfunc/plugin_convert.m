function plugin = plugin_convert(pluginOri)

for iRow = 1:length(pluginOri)
    plugin(iRow).currentversion = pluginOri(iRow).version;
    plugin(iRow).foldername     = pluginOri(iRow).foldername;
    plugin(iRow).status         = pluginOri(iRow).status;
    plugin(iRow).name           = pluginOri(iRow).plugin;
    plugin(iRow).installed      = 1;
end;