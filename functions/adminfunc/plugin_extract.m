function restartEeglabFlag = plugin_extract(type, pluginlist, page)

if nargin < 3, page = 1; end;
% type may be 'import' or 'process'
restartEeglabFlag = false;
pluginsPerPage    = 15;

% check the presence of unzip
%str = evalc('!unzip');
%if length(str) < 200
%    error([ '"unzip" could not be found. Instal unzip and make sure' 10 'it is accessible under Matlab by adding the program to' 10 'the path and typing "!unzip"' ]);
%end;

if ~isstruct(type)
    plugin = plugin_getweb(type, pluginlist, 'newlist');
    % sort plugins by download score
    [tmp scoreOrder] = sort([ plugin.downloads ], 2, 'descend');
    plugin = plugin(scoreOrder);
else
    plugin    = type;
end;

% select page
allPlugins = plugin;
numPlugin  = length(plugin);
moreThanOnePage = 0;
if numPlugin > pluginsPerPage
    plugin = plugin(pluginsPerPage*(page-1)+1:min(length(plugin),pluginsPerPage*page));
    moreThanOnePage = 1;
end;

% find which menu to show
newInstallFlag  = false;
installedFlag   = false;
deactivatedFlag = false;
for iPlugin = length(plugin):-1:1
    if ~plugin(iPlugin).installed, newInstallFlag = true; end;
    if plugin(iPlugin).installed && ~strcmpi(plugin(iPlugin).status, 'deactivated'),  installedFlag  = true; end;
    if strcmpi(plugin(iPlugin).status, 'deactivated'), deactivatedFlag = true; end;
end;

uilist   = {};
geom     = {};
geomvert = [];
pluginIndices = [];
callback = [ 'tmptag = get(gcbo, ''tag'');' ...
             'if tmptag(3) == ''1'', tmptag(3) = ''2''; else tmptag(3) = ''1''; end;' ...
             'if get(gcbo, ''value''), set(findobj(gcbf, ''tag'', tmptag), ''value'', 0); end; clear tmptag;' ];

% ------------------
% plugins to install
% ------------------
maxchar = 60;
geom    = {};
lineGeom = [ 0.28 0.28 0.95 0.6 0.6 3 0.35 ];
if newInstallFlag
    uilist =  { {} { 'style' 'text' 'string' 'Extensions available for install on the internet' 'fontweight' 'bold' 'fontsize' 18 'tag' 'title' } };
    uilist =  { uilist{:} { 'style' 'text' 'string' 'I' 'tag' 'install' } { } ...
        { 'style' 'text' 'string' 'Plugin' 'fontweight' 'bold' } ...
        { 'style' 'text' 'string' 'Vers.'     'tag' 'verweb'     'fontweight' 'bold' } ...
        { 'style' 'text' 'string' 'Score' 'fontweight' 'bold' } ...
        { 'style' 'text' 'string' 'Description' 'fontweight' 'bold' } {}};
    geom = { [1 5.5] lineGeom };
    geomvert = [1 1];
    for iRow = 1:length(plugin)
        if ~plugin(iRow).installed && ~strcmpi(plugin(iRow).status, 'deactivated')
            % text for description
            description = plugin(iRow).description;
            if length(description) > maxchar+2
                description = [ description(1:min(maxchar,length(description))) '...' ];
            end;
            
            enableWebDoc = fastif(isempty(plugin(iRow).webdoc), 'off', 'on');
            
            userdata = '';
            if plugin(iRow).installed && plugin(iRow).installorupdate, userdata = 'colortored'; end;
            uilist = { uilist{:}, ...
                { 'style' 'checkbox' 'string' '' 'value'  0 'enable' 'on' }, ...
                { 'style' 'checkbox' 'string' '' 'visible' 'off' }, ...
                { 'style' 'text' 'string' plugin(iRow).name }, ...
                { 'style' 'text' 'string' plugin(iRow).version 'tag' 'latestversion' 'userdata' userdata }, ...
                { 'style' 'text' 'string' int2str(plugin(iRow).downloads) }, ...
                { 'style' 'text' 'string' description }, ...
                { 'style' 'pushbutton' 'string' 'Doc' 'enable' enableWebDoc 'callback' myweb(plugin(iRow).webdoc) } };
            geom = { geom{:}, lineGeom };
            geomvert = [ geomvert 1];
            pluginIndices = [ pluginIndices iRow ];
        end;
    end;
end;

% -----------------
% installed plugins
% -----------------
if installedFlag
    uilist =  { uilist{:} {} {} { 'style' 'text' 'string' 'Installed extensions' 'fontweight' 'bold' 'fontsize' 18 'tag' 'title' } };
    uilist =  { uilist{:} { 'style' 'text' 'string' 'I' 'tag' 'update' } ...
        { 'style' 'text' 'string' 'I' 'tag' 'deactivate' } ...
        { 'style' 'text' 'string' 'Plugin' 'fontweight' 'bold' } ...
        { 'style' 'text' 'string' 'Vers.'     'tag' 'verweb'     'fontweight' 'bold' } ...
        { 'style' 'text' 'string' 'Score' 'fontweight' 'bold' } ...
        { 'style' 'text' 'string' 'Description' 'fontweight' 'bold' } {}};
    
    geom = { geom{:} 1 [1 5.5] lineGeom };
    geomvert = [geomvert 1 1 1];
    for iRow = 1:length(plugin)
        if plugin(iRow).installed && ~strcmpi(plugin(iRow).status, 'deactivated')
            % text for description
            description = plugin(iRow).description;
            if length(description) > maxchar+2
                description = [ description(1:min(maxchar,length(description))) '...' ];
            end;
            
            enableWebDoc = fastif(isempty(plugin(iRow).webdoc), 'off', 'on');
            
            userdata = '';
            if plugin(iRow).installorupdate,
                textnew = [ 'Click update to install version ' plugin(iRow).version ' now available on the web' ];
                userdata = 'colortored';
            else
                textnew = description;
            end;
            uilist = { uilist{:}, ...
                { 'style' 'checkbox' 'string' '' 'value' 0 'enable' fastif(plugin(iRow).installorupdate, 'on', 'off') 'tag' [ 'cb1' int2str(iRow) ] 'callback' callback }, ...
                { 'style' 'checkbox' 'string' '' 'enable' 'on' 'tag' [ 'cb2' int2str(iRow) ] 'callback' callback }, ...
                { 'style' 'text' 'string' plugin(iRow).name }, ...
                { 'style' 'text' 'string' plugin(iRow).currentversion 'tag' 'latestversion'}, ...
                { 'style' 'text' 'string' int2str(plugin(iRow).downloads) }, ...
                { 'style' 'text' 'string' textnew 'userdata' userdata }, ...
                { 'style' 'pushbutton' 'string' 'Doc' 'enable' enableWebDoc 'callback' myweb(plugin(iRow).webdoc) } };
            geom = { geom{:}, lineGeom };
            geomvert = [ geomvert 1];
            pluginIndices = [ pluginIndices iRow ];
        end;
    end;
end;

% -------------------
% deactivated plugins
% -------------------
%geom = { geom{:} 1 1 };
%geomvert = [geomvert 0.5 1];
%uilist = { uilist{:} {} { 'style' 'text' 'string' 'To manage deactivated plugins, use menu item File > Manage plugins > Manage deactivated plugins' } };              
if deactivatedFlag
    uilist =  { uilist{:} {} {} { 'style' 'text' 'string' 'List of deactivated extensions                         ' 'fontweight' 'bold' 'fontsize' 18 'tag' 'title' } };
    uilist =  { uilist{:} ...
        { 'style' 'text' 'string' 'I' 'tag' 'reactivate' } ...
        { 'style' 'text' 'string' 'I' 'tag' 'remove1' } ...
        { 'style' 'text' 'string' 'Plugin' 'fontweight' 'bold' } ...
        { 'style' 'text' 'string' 'Vers.'     'tag' 'verweb'     'fontweight' 'bold' } ...
        { 'style' 'text' 'string' 'Score' 'fontweight' 'bold' } ...
        { 'style' 'text' 'string' 'Description' 'fontweight' 'bold' } {}};
    geom = { geom{:} 1 [1 5.5] lineGeom };
    geomvert = [geomvert 1 1 1];
    for iRow = 1:length(plugin)
        if strcmpi(plugin(iRow).status, 'deactivated')
            % text for description
            description = plugin(iRow).description;
            if length(description) > maxchar+2
                description = [ description(1:min(maxchar,length(description))) '...' ];
            end;
            
            userdata = '';
            enableWebDoc = fastif(isempty(plugin(iRow).webdoc), 'off', 'on');
            uilist = { uilist{:}, ...
                { 'style' 'checkbox' 'string' '' 'tag' [ 'cb1' int2str(iRow) ] 'callback' callback }, ...
                { 'style' 'checkbox' 'string' '' 'tag' [ 'cb2' int2str(iRow) ] 'callback' callback }, ...
                { 'style' 'text' 'string' plugin(iRow).name }, ...
                { 'style' 'text' 'string' plugin(iRow).version 'tag' 'latestversion' }, ...
                { 'style' 'text' 'string' int2str(plugin(iRow).downloads) }, ...                
                { 'style' 'text' 'string' description }, ...
                { 'style' 'pushbutton' 'string' 'Doc' 'enable' enableWebDoc 'callback' myweb(plugin(iRow).webdoc) } };
            geom = { geom{:}, lineGeom };
            geomvert = [ geomvert 1];
            pluginIndices = [ pluginIndices iRow ];
        end;
    end;
end;

evalStr = [ 'uisettxt(gcf, ''update''         , ''Update''     , ''rotation'', 90, ''fontsize'', 14);' ...
            'uisettxt(gcf, ''deactivate''     , ''Deactivate'' , ''rotation'', 90, ''fontsize'', 14);' ...
            'uisettxt(gcf, ''install''        , ''Install''    , ''rotation'', 90, ''fontsize'', 14);' ...
            'uisettxt(gcf, ''reactivate''     , ''Reactivate'' , ''rotation'', 90, ''fontsize'', 14);' ...
            'uisettxt(gcf, ''remove1''        , ''Remove''     , ''rotation'', 90, ''fontsize'', 14);' ...
            'set(findobj(gcf, ''tag'', ''title''), ''fontsize'', 16);' ...
            'tmpobj = findobj(gcf, ''userdata'', ''colortored'');' ...
            'set(tmpobj, ''Foregroundcolor'', [1 0 0]);' ...
            'tmppos = get(gcf, ''position'');' ...
            'set(gcf, ''position'', [tmppos(1:2) 800 tmppos(4)]);' ...
            'clear tmpobj tmppos;' ...
            ];
        
if 1
    % version with button
    if page == 1,                       enablePpage = 'off'; else enablePpage = 'on'; end;
    if page*pluginsPerPage > numPlugin, enableNpage = 'off'; else enableNpage = 'on'; end;
    callBackPpage = [ 'tmpobj = get(gcbf, ''userdata''); close gcbf; restartEeglabFlag = plugin_extract(tmpobj, [], ' int2str(page-1) '); clear tmpobj;' ];
    callBackNpage = [ 'tmpobj = get(gcbf, ''userdata''); close gcbf; restartEeglabFlag = plugin_extract(tmpobj, [], ' int2str(page+1) '); clear tmpobj;' ];
    
    uilist = { uilist{:}, {} { 'width' 80 'align' 'left'  'Style', 'pushbutton', 'string', '< Prev. page', 'tag' 'ppage' 'callback', callBackPpage 'enable' enablePpage } };
    uilist = { uilist{:},    { 'width' 80 'align' 'left'  'stickto' 'on', 'Style', 'pushbutton', 'string', 'Next page >', 'tag' 'npage' 'callback', callBackNpage 'enable' enableNpage } };
    uilist = { uilist{:},    { 'width' 80 'align' 'right' 'Style', 'pushbutton', 'string', 'Cancel', 'tag' 'cancel' 'callback', 'close gcbf' } };
    uilist = { uilist{:},    { 'width' 80 'align' 'right' 'stickto' 'on' 'Style', 'pushbutton', 'tag', 'ok', 'string', 'OK', 'callback', 'set(gcbo, ''userdata'', ''retuninginputui'');' } };
    geom     = { geom{:} [1] [1 1 1 1] };
    geomvert = [ geomvert 1 1];    
    res = inputgui('uilist', uilist, 'geometry', geom, 'geomvert', geomvert, 'eval', evalStr, 'addbuttons', 'off', 'skipline', 'off', 'userdata', allPlugins);
    
    try, restartEeglabFlag = evalin('base', 'restartEeglabFlag;'); catch, end;
    evalin('base', 'clear restartEeglabFlag;');
else
    % no buttons
    res = inputgui('uilist', uilist, 'geometry', geom, 'geomvert', geomvert, 'eval', evalStr);
end;
if isempty(res), return; end;

% decode inputs
% -------------
for iRow = 1:length(pluginIndices)
    plugin(pluginIndices(iRow)).install = res{(iRow-1)*2+1};
    plugin(pluginIndices(iRow)).remove  = res{(iRow-1)*2+2};
end;

% install plugins
% ---------------
firstPlugin = 1;
for iRow = 1:length(plugin)
    if plugin(iRow).install
        restartEeglabFlag = true;
        if ~firstPlugin, disp('---------------------------------'); end; firstPlugin = 0;
        
        if strcmpi(plugin(iRow).status, 'deactivated')
            fprintf('Reactivating extension %s\n', plugin(iRow).name);
            plugin_reactivate(plugin(iRow).foldername);
            if plugin(iRow).installorupdate
                res = questdlg2([ 'Extension ' plugin(iRow).foldername ' has been reactivated but' 10 'a new version is available. Do you want to install it?' ], 'Warning', 'No', 'Yes', 'Yes');
                if strcmpi(res, 'yes')
                    plugin_deactivate(plugin(iRow).foldername);
                    plugin_install(plugin(iRow).zip, plugin(iRow).name, plugin(iRow).version);
                    plugin_remove(plugin(iRow).foldername);
                end;
            end;
        else
            if plugin(iRow).installed
                fprintf('Updating extension %s\n', plugin(iRow).name);
                plugin_deactivate(plugin(iRow).foldername);
                plugin_install(plugin(iRow).zip, plugin(iRow).name, plugin(iRow).version);
                plugin_remove(plugin(iRow).foldername);
            else
                fprintf('Installing extension %s\n', plugin(iRow).name);
                plugin_install(plugin(iRow).zip, plugin(iRow).name, plugin(iRow).version);
            end;
        end;
    elseif plugin(iRow).remove
        if ~firstPlugin, disp('---------------------------------'); end; firstPlugin = 0; 
        restartEeglabFlag = true;
        
        if strcmpi(plugin(iRow).status, 'deactivated')
            fprintf('Removing extension %s\n', plugin(iRow).name);
            plugin_remove(plugin(iRow).foldername);
        else
            fprintf('Deactivating extension %s\n', plugin(iRow).name);
            plugin_deactivate(plugin(iRow).foldername);
        end;
    end;
end;

function str = myweb(url);

    %if isempty(strfind(url, 'wiki'))
    %     str = [ 'web(''' url ''');' ];
    %else
        str = [ 'web(''' url ''', ''-browser'');' ];
    %end;
    
