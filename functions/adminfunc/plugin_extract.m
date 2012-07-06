function res = plugin_extract(type);

% type may be 'import' or 'process'
res = 0;
global PLUGINLIST;

% check the presence of unzip
%str = evalc('!unzip');
%if length(str) < 200
%    error([ '"unzip" could not be found. Instal unzip and make sure' 10 'it is accessible under Matlab by adding the program to' 10 'the path and typing "!unzip"' ]);
%end;

try
    disp('Retreiving URL with plugins...');
    if strcmpi(type, 'import')
        [tmp status] = urlread('http://sccn.ucsd.edu/wiki/Plugin_list_import');
    else
        [tmp status] = urlread('http://sccn.ucsd.edu/wiki/Plugin_list_process');
    end;
catch,
    error('Cannot connect to the Internet to retrieve plugin list');
end;
if status == 0
    error('Cannot connect to the Internet to retrieve plugin list');
end;    
    
% parse the web page
% ------------------
try
    plugin = parseTable(tmp);
catch
    disp('PLUGIN PAGE PARSING ERROR, USING BACKUP PLUGIN LIST');
    [tmp status] = urlread('http://sccn.ucsd.edu/wiki/Plugin_backup_list');
    plugin = parseTable(tmp);
end;

% find correspondance with plugin list
% ------------------------------------
currentNames = lower({ PLUGINLIST.plugin });
for iRow = 1:length(plugin)
    indMatch = strmatch(lower(plugin(iRow).name), currentNames, 'exact');
    if isempty(indMatch)
        plugin(iRow).currentversion  = '-';
        plugin(iRow).installed       = 0;
        plugin(iRow).installorupdate = 1;
    else
        plugin(iRow).currentversion = PLUGINLIST(indMatch).version;
        plugin(iRow).foldername     = PLUGINLIST(indMatch).foldername;
        plugin(iRow).installed      = 1;
        if strcmpi(plugin(iRow).currentversion, plugin(iRow).version)
            plugin(iRow).installorupdate = 0;
        else
            plugin(iRow).installorupdate = 1;
        end;
    end;
end;

% put all the installed plugins first
% -----------------------------------
[tmp reorder] = sort([plugin.installed], 'descend');
plugin = plugin(reorder);

% name
% description
% uninstall
% install
% ignore
% update available
%% create GUI
uilist =  { { 'style' 'text' 'string' 'I' 'tag' 'install' } ...
            { 'style' 'text' 'string' 'I' 'tag' 'problemo' } ...
            { 'style' 'text' 'string' 'I' 'tag' 'ignore' } ...
            { 'style' 'text' 'string' 'Plugin' 'fontweight' 'bold' } ...
            { 'style' 'text' 'string' 'Vers.' 'tag' 'versheader' 'fontweight' 'bold' } ...
            { 'style' 'text' 'string' 'Vers.' 'tag' 'verweb'     'fontweight' 'bold' } ...
            { 'style' 'text' 'string' 'Description' 'fontweight' 'bold' } {}};
lineGeom = [ 0.25 0.25 0.4 0.8 0.4 0.4 3 0.35 ];
geom = { lineGeom };
geomvert = [2];
maxchar = 60;
for iRow = 1:length(plugin)
    % text for description
    description = plugin(iRow).description;
    if length(description) > maxchar+2
         description = [ description(1:min(maxchar,length(description))) '...' ];
    end;
    
    enableWebDoc = fastif(isempty(plugin(iRow).webdoc), 'off', 'on');

    userdata = '';
    if plugin(iRow).installed && plugin(iRow).installorupdate, userdata = 'colortored'; end;
    uilist = { uilist{:}, ...
              { 'style' 'checkbox' 'string' '' 'value' plugin(iRow).installed 'enable' 'off' }, ...
              { 'style' 'checkbox' 'string' '' 'enable' fastif(plugin(iRow).installed, 'on', 'off') }, ...
              { 'style' 'checkbox' 'string' '' 'enable' fastif(plugin(iRow).installorupdate, 'on', 'off') }, ...}, ...
              { 'style' 'text' 'string' plugin(iRow).name }, ...
              { 'style' 'text' 'string' plugin(iRow).currentversion 'tag' 'currentversion' }, ...
              { 'style' 'text' 'string' plugin(iRow).version 'tag' 'latestversion' 'userdata' userdata }, ...
              { 'style' 'text' 'string' description }, ...
              { 'style' 'pushbutton' 'string' 'Doc' 'enable' enableWebDoc 'callback' [ 'web(''' plugin(iRow).webdoc ''');' ] } };              
    geom = { geom{:}, lineGeom };
    geomvert = [ geomvert 1];
end;

% evalStr = [ 'tmpobj = findobj(gcf, ''tag'', ''version'');' ...
%             'set(tmpobj, ''Foregroundcolor'', [1 0 0]);' ...
%             'tmpobj = findobj(gcf, ''tag'', ''testi'');' ...
%             'tmppos = get(tmpobj, ''position'');' ...
%             'delete(tmpobj);' ...
%             'axes(''position'', tmppos);' ...
%             'axis(''off'');' ...
%             'tmp = text(5,0, ''Installed'');' ... 
%             'xlim([0 10]);' ...
%             'set(tmp, ''rotation'', 90, ''fontweight'', ''bold'');' ];

evalStr = [ 'tmpobj = findobj(gcf, ''userdata'', ''colortored'');' ...
            'set(tmpobj, ''Foregroundcolor'', [1 0 0]);' ...
            'putverticaltext(gcf, ''versheader''     , {''Installed'' ''version'' });' ...
            'putverticaltext(gcf, ''verweb''         , {''Latest''    ''version'' });' ...
            'putverticaltext(gcf, ''problemo''     , ''Remove'');' ...
            'putverticaltext(gcf, ''install'', ''Installed'');' ...
            'putverticaltext(gcf, ''ignore'',  { ''Install or'' ''Update'' });' ...
            ];

res = inputgui('uilist', uilist, 'geometry', geom, 'geomvert', geomvert, 'eval', evalStr);
if isempty(res), return; end;

% decode inputs
% -------------
for iRow = 1:length(plugin)
    plugin(iRow).install = res{(iRow-1)*3+3};
    plugin(iRow).remove  = res{(iRow-1)*3+2};
end;

% install plugins
% ---------------
firstPlugin = 1;
for iRow = 1:length(plugin)
    if plugin(iRow).install
        res = 1;
        if ~firstPlugin
            disp('---------------------------------');
        end;
        firstPlugin = 0;
        res = plugin_install(plugin(iRow).zip, plugin(iRow).name, plugin(iRow).version);
        if plugin(iRow).installed && res == 1
            plugin_remove(plugin(iRow).foldername);
        end;
    elseif plugin(iRow).remove
        res = 1;
        plugin_remove(plugin(iRow).foldername);
    end;
end;

% parse the web table
% ===================
function plugin = parseTable(tmp);

% get table content
% -----------------
tableBeg = findstr('Plugin name', tmp);
tableEnd = findstr('</table>', tmp(tableBeg:end));
tableContent = tmp(tableBeg:tableBeg+tableEnd-2);
endFirstLine = findstr('</tr>', tableContent);
tableContent = tableContent(endFirstLine(1)+5:end);

% parse table entries
% -------------------
posBegRow = findstr('<tr>' , tableContent);
posEndRow = findstr('</tr>', tableContent);
if length(posBegRow) ~= length(posEndRow) || isempty(posBegRow)
    error('Cannot connect to the Internet to retrieve plugin list');
end;
for iRow = 1:length(posBegRow)
    rowContent = tableContent(posBegRow(iRow)+4:posEndRow(iRow)-1);
    posBegCol = findstr('<td>' , rowContent);
    posEndCol = findstr('</td>', rowContent);
    for iCol = 1:length(posBegCol)
        table{iRow,iCol} = rowContent(posBegCol(iCol)+4:posEndCol(iCol)-1);
    end;
end;

%% extract zip link and plugin name from first column
% --------------------------------------------------
%href="http://www.unicog.org/pm/uploads/MEG/ADJUST_PLUGIN.zip" class="external text" title="http://www.unicog.org/pm/uploads/MEG/ADJUST_PLUGIN.zip" rel="nofollow">ADJUST PLUGIN</a></td
for iRow = 1:size(table,1)
    
    % get link
    [plugin(iRow).name plugin(iRow).webdoc] = parsehttplink(table{iRow,1});
    plugin(iRow).version = table{iRow,2};
    plugin(iRow).description = deblank(table{iRow,3});
    [tmp plugin(iRow).zip] = parsehttplink(table{iRow,4});
    
end;

function [txt link] = parsehttplink(currentRow)
    openTag  = find(currentRow == '<');
    closeTag = find(currentRow == '>');
    if isempty(openTag)
        link = '';
        txt = currentRow;
    else
        % parse link
        link = currentRow(openTag(1)+1:closeTag(1)-1);
        hrefpos = findstr('href', link);
        link = link(hrefpos:end);
        quoteInd = find(link == '"');
        link = link(quoteInd(1)+1:quoteInd(2)-1);
        
        for iTag = length(openTag):-1:1
            currentRow(openTag(iTag):closeTag(iTag)) = [];
        end;
        txt = currentRow;
    end;
