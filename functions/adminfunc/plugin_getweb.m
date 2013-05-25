function plugin = plugin_getweb(type, pluginOri, mode)

if nargin < 1, help plugin_getweb; return; end;
if nargin < 2, pluginOri = []; end;
if nargin < 3, mode = 'merge'; end; % 'merge' or 'newlist'

% convert plugin list format if necessary
if isfield(pluginOri, 'plugin'), pluginOri = plugin_convert(pluginOri); end;

try
    disp( [ 'Retreiving URL with ' type ' plugins...' ] );
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
if ~isempty(pluginOri)
     currentNames = lower({ pluginOri.name });
else currentNames = {};
end;
allMatch = [];
for iRow = 1:length(plugin)
    indMatch = strmatch(lower(plugin(iRow).name), currentNames, 'exact');
    if isempty(indMatch)
        plugin(iRow).currentversion  = '-';
        plugin(iRow).installed       = 0;
        plugin(iRow).installorupdate = 1;
        plugin(iRow).status          = 'notinstalled';
    else
        plugin(iRow).currentversion = pluginOri(indMatch).currentversion;
        plugin(iRow).foldername     = pluginOri(indMatch).foldername;
        plugin(iRow).status         = pluginOri(indMatch).status;
        plugin(iRow).installed      = 1;
        if strcmpi(plugin(iRow).currentversion, plugin(iRow).version)
            plugin(iRow).installorupdate = 0;
        else
            plugin(iRow).installorupdate = 1;
        end;
        allMatch = [ allMatch indMatch ];
    end;
end;

% put all the installed plugins first
% -----------------------------------
if ~isempty(plugin)
    [tmp reorder] = sort([plugin.installed], 'descend');
    plugin = plugin(reorder);
%     plugin(1).currentversion  = '0.9';
%     plugin(1).version         = '1';
%     plugin(1).foldername      = 'test';
%     plugin(1).installed       = 1;
%     plugin(1).installorupdate = 1;
%     plugin(1).description     = 'test';
%     plugin(1).webdoc          = 'test';
%     plugin(1).name            = 'test';
end;

if strcmpi(mode, 'merge')
    indices = setdiff([1:length(pluginOri)], allMatch);
    fields  = fieldnames(pluginOri);
    lenPlugin = length(plugin);
    
    for indPlugin = 1:length(indices)
        for indField = 1:length(fields)
            value  = getfield(pluginOri, { indices(indPlugin)  }, fields{ indField });
            plugin = setfield(plugin   , { lenPlugin+indPlugin }, fields{ indField }, value);
        end;
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