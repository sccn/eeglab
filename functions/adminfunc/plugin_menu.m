function restartEeglabFlag = plugin_menu(pluginlist)

% type may be 'import' or 'process'
restartEeglabFlag = false;
pluginsPerPage    = 15;

% check the presence of unzip
%str = evalc('!unzip');
%if length(str) < 200
%    error([ '"unzip" could not be found. Instal unzip and make sure' 10 'it is accessible under Matlab by adding the program to' 10 'the path and typing "!unzip"' ]);
%end

plugin = plugin_getweb('', pluginlist, 'newlist');

% sort plugins by download score
[~,scoreOrder] = sort([ plugin.downloads ], 2, 'descend');
plugin = plugin(scoreOrder);

% sort plugins by download score
[~, scoreOrder] = sort([ plugin.downloads ], 2, 'descend');
plugin = plugin(scoreOrder);

uilist   = {};
geom     = {};
geomvert = [];
pluginIndices = [];
         
% ------------------
% plugins to install
% ------------------
maxchar = 60;
geom    = {};

lsitboxtext = {};
for iRow = 1:length(plugin)
    plugin(iRow).text = [ '<html><font size=+0> ' htmlrating(plugin(iRow).rating, plugin(iRow).numrating) ];
    if plugin(iRow).installed
        if plugin(iRow).installorupdate
            plugin(iRow).text = [ plugin(iRow).text '<b><font color=red>' ];
        else
            plugin(iRow).text = [ plugin(iRow).text '<b>' ];
        end
    end
    plugin(iRow).text =  [ plugin(iRow).text plugin(iRow).name ' v'  plugin(iRow).version ' (' int2str(plugin(iRow).downloads) ' downloads' ];
    if plugin(iRow).numrating
        plugin(iRow).text =  [ plugin(iRow).text '; ' int2str(plugin(iRow).numrating) ' rating' ];
    end
    plugin(iRow).text =  [ plugin(iRow).text ')</b></font></font></html>' ];
end

%cb_select = 'tmpobj = get(gcbf, ''userdata''); tmpstr = tmpobj(get(gcbo, ''value'')).longdescription; tmpstr = textwrap(findobj(gcbf, ''tag'', ''description''), {tmpstr}); set(findobj(gcbf, ''tag'', ''description''), ''string'', tmpstr); clear tmpobj tmpstr;';
filterList = { 'No filter' ...
               'Filter by import' ...
               'Filter by export' ...
               'Filter by artifact' ...
               'Filter by ica'  ...
               'Filter by preprocessing' ...
               'Filter by erp' ...
               'Filter by source' ...
               'Filter by study' ...
               'Filter by time-freq' ...
               'Filter by other' };

uilist =  {
    { 'style', 'text', 'string', 'List of plugins (bolded plugins are installed)' 'fontweight' 'bold' } ...
    { 'style', 'popupmenu', 'string', filterList  'callback' 'pluguin_uifilter(gcbf);' 'tag' 'filter' } ...
    { 'style', 'listbox', 'string', { plugin.text } 'callback' 'pluguin_uiupdate(gcbf);' 'Min', 0, 'Max', 2, 'value' [] 'tag', 'pluginlist' 'fontsize', 16 } ...
    { 'style', 'pushbutton', 'string', [ 'Rate plugin' ] 'tag' 'rating' } ...
    { 'style', 'pushbutton', 'string', [ 'Web documentation' ] 'tag' 'documentation' } ...
    { 'style', 'text', 'string', 'Tags:' 'fontweight' 'bold' } ...
    { 'style', 'text', 'string', 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 'tag' 'tags' } ...
    { 'style', 'text', 'string', 'Status:' 'fontweight' 'bold' } ...
    { 'style', 'text', 'string', 'installed' 'tag' 'status'} ...
    { 'style', 'text', 'string', 'Description:' 'fontweight' 'bold' } ...
    { 'style', 'text', 'string', [ 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 10 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 10 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 10 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' ] 'tag' 'description' } ...
    {} ...
    {} ...
    { 'Style', 'pushbutton', 'string', 'Cancel', 'tag' 'cc' 'callback', 'close gcbf' } ...
    { 'Style', 'pushbutton', 'tag', 'remove', 'string', 'Remove', 'callback', 'set(gcbo, ''userdata'', ''retuninginputui'');' } ...
    { 'Style', 'pushbutton', 'tag', 'install', 'string', 'Install/Update', 'callback', 'set(gcbo, ''userdata'', ''retuninginputui'');' } ...
    };

usrDat.allplugins = plugin;
usrDat.selectedplugins = plugin;
usrDat.selection = [];
fig = figure('visible', 'off');
supergui('fig', fig, 'uilist', uilist, 'geomhoriz', {[1 0.5] 1 [1 1] [0.2 1] [0.2 1] 1 1 1 [0.43 0.37 0.4 0.5]}, 'geomvert', [1 10 1 1 1 1 2.5 1 1], 'userdata', usrDat);
pos = get(fig, 'position');
set(fig, 'position', [pos(1) pos(2) pos(3)/841*200 pos(4) ]);
set(fig, 'visible', 'on');
waitfor( findobj('parent', fig, 'tag', 'install'), 'userdata');

% Cancel
if ~(ishandle(fig)), return; end % Check if figure still exist

usrDat = get(fig, 'userdata');
close(fig);

plugin = usrDat.selectedplugins;
for iSelect = 1:length(usrDat.selection)
    plugin(iSelect).install = true;
end

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
                    if plugin_install(plugin(iRow).zip, plugin(iRow).name, plugin(iRow).version) == -1
                        plugin_reactivate(plugin(iRow).foldername);
                    else
                        plugin_remove(plugin(iRow).foldername);
                    end
                end
            end
        else
            if plugin(iRow).installed
                fprintf('Updating extension %s\n', plugin(iRow).name);
                plugin_deactivate(plugin(iRow).foldername);
                if plugin_install(plugin(iRow).zip, plugin(iRow).name, plugin(iRow).version) == -1
                    plugin_reactivate(plugin(iRow).foldername);
                else
                    plugin_remove(plugin(iRow).foldername);
                end
            else
                fprintf('Installing extension %s\n', plugin(iRow).name);
                plugin_install(plugin(iRow).zip, plugin(iRow).name, plugin(iRow).version);
            end
        end
    elseif plugin(iRow).remove
        if ~firstPlugin, disp('---------------------------------'); end; firstPlugin = 0; 
        restartEeglabFlag = true;
        
        if strcmpi(plugin(iRow).status, 'deactivated')
            fprintf('Removing extension %s\n', plugin(iRow).name);
            plugin_remove(plugin(iRow).foldername);
        else
            fprintf('Deactivating extension %s\n', plugin(iRow).name);
            plugin_deactivate(plugin(iRow).foldername);
        end
    end
end

function str = htmlrating(rating, numRating)
    roundRating = round(rating);
    if isempty(roundRating) || isnan(roundRating) || numRating == 0
        str = 'no rating &nbsp;- ';
    else
        str = '<font size=-1><font color="#CCCC00">';
        for i=1:roundRating
            str = [ str '<span>&#9733;</span>' ];
        end
        str = [ str '</font><font color="gray">' ];
        for i=roundRating+1:5
            str = [ str '<span>&#9733;</span>' ];
        end
        str = [ str '</font></font> - ' ];
    end
    
function str = myweb(url);

    %if isempty(strfind(url, 'wiki'))
    %     str = [ 'web(''' url ''');' ];
    %else
        str = [ 'web(''' url ''', ''-browser'');' ];
    %end
