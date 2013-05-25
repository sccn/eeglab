function restartEeglabFlag = plugin_managedeactivated(pluginlist)

% type may be 'import' or 'process'
restartEeglabFlag = false;
plugin = plugin_getweb('import' , pluginlist);
plugin = plugin_getweb('process', plugin);

for iRow = 1:length(plugin)
    plugin(iRow).remove       = 0;
    plugin(iRow).reactivate   = 0;
end;

if isempty(strmatch('deactivated', { plugin.status }, 'exact'))
    warndlg2('There are no deactivated plugins');
    return;
end;

maxchar = 60;
uilist =  { {} { 'style' 'text' 'string' '             List of deactivated plugins                                 ' 'fontweight' 'bold' 'fontsize' 18 'tag' 'title' } };
uilist =  { uilist{:} ...
            { 'style' 'text' 'string' 'I' 'tag' 'reactivate' } ...
            { 'style' 'text' 'string' 'I' 'tag' 'remove1' } ...
            { 'style' 'text' 'string' 'Plugin' 'fontweight' 'bold' } ...
            { 'style' 'text' 'string' 'Version'    'tag' 'verweb'     'fontweight' 'bold' } ...
            { 'style' 'text' 'string' 'Description' 'fontweight' 'bold' } {}};
lineGeom = [ 0.28 0.28 0.95 0.8 3 0.35 ];
geom = { [1 5.5] lineGeom };
geomvert = [1 1];
pluginIndices = [];
for iRow = 1:length(plugin)
    if strcmpi(plugin(iRow).status, 'deactivated') 
        % text for description
        description = plugin(iRow).description;
        if length(description) > maxchar+2
             description = [ description(1:min(maxchar,length(description))) '...' ];
        end;

        cb = [ 'tmptag = get(gcbo, ''tag'');' ...
               'if tmptag(3) == ''1'', tmptag(3) = ''2''; else tmptag(3) = ''1''; end;' ...
               'set(findobj(gcbf, ''tag'', tmptag), ''value'', ~get(gcbo, ''value''));' ];
        userdata = '';
        enableWebDoc = fastif(isempty(plugin(iRow).webdoc), 'off', 'on');
        uilist = { uilist{:}, ...
                  { 'style' 'checkbox' 'string' '' 'tag' [ 'cb1' int2str(iRow) ] 'callback' cb }, ...
                  { 'style' 'checkbox' 'string' '' 'tag' [ 'cb2' int2str(iRow) ] 'callback' cb }, ...
                  { 'style' 'text' 'string' plugin(iRow).name }, ...
                  { 'style' 'text' 'string' plugin(iRow).version 'tag' 'latestversion' }, ...
                  { 'style' 'text' 'string' description }, ...
                  { 'style' 'pushbutton' 'string' 'Doc' 'enable' enableWebDoc 'callback' [ 'web(''' plugin(iRow).webdoc ''');' ] } };              
        geom = { geom{:}, lineGeom };
        geomvert = [ geomvert 1];
        pluginIndices = [ pluginIndices iRow ];
    end;
end;

evalStr = [ 'uisettxt(gcf, ''reactivate''     , ''Reactivate'' , ''rotation'', 90, ''fontsize'', 14);' ...
            'uisettxt(gcf, ''remove1''        , ''Remove''     , ''rotation'', 90, ''fontsize'', 14);' ...
            'tmppos = get(gcf, ''position''); set(gcf, ''position'', [tmppos(1:2) max(700, tmppos(3)) tmppos(4)]);' ...
            'set(findobj(gcf, ''tag'', ''title''), ''fontsize'', 16);' ...
            'tmpobj = findobj(gcf, ''userdata'', ''colortored'');' ...
            'set(tmpobj, ''Foregroundcolor'', [1 0 0]);' ...
            'clear tmppos tmpobj;'];

res = inputgui('uilist', uilist, 'geometry', geom, 'geomvert', geomvert, 'eval', evalStr);
if isempty(res), return; end;

% decode inputs
% -------------
for iRow = 1:length(pluginIndices)
    plugin(pluginIndices(iRow)).reactivate  = res{(iRow-1)*2+1};
    plugin(pluginIndices(iRow)).remove      = res{(iRow-1)*2+2};
end;

% install plugins
% ---------------
firstPlugin = 1;
for iRow = 1:length(plugin)
    if plugin(iRow).remove
        if ~firstPlugin, disp('---------------------------------'); end; firstPlugin = 0;
        
        fprintf('Removing plugin %s\n', plugin(iRow).name);
        plugin_remove(plugin(iRow).foldername);
    elseif plugin(iRow).reactivate
        if ~firstPlugin, disp('---------------------------------'); end; firstPlugin = 0; 
        restartEeglabFlag = true;
        fprintf('Reactivating plugin %s\n', plugin(iRow).name);
        plugin_reactivate(plugin(iRow).foldername);
    end;
end;
