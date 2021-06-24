% plugin_menu() - main function to install EEGLAB plugins
%
% To install plugins from the command line, type in
%
% plugin_askinstall('xxxxxx', [], true); % with xxxx being the name of the plugin
%
% Usage: plugin_menu(PLUGINLIST); % pop up gui

% Copyright (C) 2019 Arnaud Delorme
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function restartEeglabFlag = plugin_menu(pluginlist)

FONTSIZE = 0; % set to 4 for high-res screen

% type may be 'import' or 'process'
restartEeglabFlag = false;
plugin = plugin_getweb('plugin_check', pluginlist);
if isempty(plugin)           
    errordlg2(['Either you are offline, a firewall is blocking EEGLAB from accessing its' char(10) ...
        'plugin server or there is a problem with Java. For Java problems, refer to' char(10) ...
        'https://github.com/sccn/eeglab/issues/20']);
    return;
end

% sort plugins by download score
[~,scoreOrder] = sort([ plugin.downloads ], 2, 'descend');
plugin = plugin(scoreOrder);

% sort plugins by download score
[~, scoreOrder] = sort([ plugin.downloads ], 2, 'descend');
plugin = plugin(scoreOrder);
        
% ------------------
% plugins to install
% ------------------
for iRow = 1:length(plugin)
    plugin(iRow).text = [ '<html><font size=+' int2str(FONTSIZE) '> ' htmlrating(plugin(iRow).rating, plugin(iRow).numrating) ];
    if plugin(iRow).installed
        if plugin(iRow).installorupdate
            plugin(iRow).text = [ plugin(iRow).text '<b><font color=red>' ];
        else
            plugin(iRow).text = [ plugin(iRow).text '<b>' ];
        end
    end
    plugin(iRow).text =  [ plugin(iRow).text plugin(iRow).name ' v'  plugin(iRow).version ];
    if plugin(iRow).installed && plugin(iRow).installorupdate
        plugin(iRow).text =  [ plugin(iRow).text ' update available ' ];
    end
    plugin(iRow).text =  [ plugin(iRow).text ' (' int2str(plugin(iRow).downloads) ' downloads' ];
    if ~isnan(plugin(iRow).numrating) && plugin(iRow).numrating
        plugin(iRow).text =  [ plugin(iRow).text '; ' int2str(plugin(iRow).numrating) ' rating' ];
    end
    plugin(iRow).text =  [ plugin(iRow).text ')</b></font></font></html>' ];
    plugin(iRow).strsearch = lower([ plugin(iRow).name plugin(iRow).rawtags plugin(iRow).description ]);  
end

%cb_select = 'tmpobj = get(gcbf, ''userdata''); tmpstr = tmpobj(get(gcbo, ''value'')).longdescription; tmpstr = textwrap(findobj(gcbf, ''tag'', ''description''), {tmpstr}); set(findobj(gcbf, ''tag'', ''description''), ''string'', tmpstr); clear tmpobj tmpstr;';
filterList1 = { 'No install status filter' ...
                'Show installed only' ...
                'Show non installed only' };
filterList2 = { 'No topic filter' ...
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
search_icon_path = ['file://' fullfile(fileparts(which('plugin_menu.m')),'search-icon.png')];           
uilist =  {
    { 'style', 'text', 'string', 'List of plugins (bolded means installed)' 'fontweight' 'bold' } ...
    { 'style', 'popupmenu', 'string', filterList1 'callback' 'plugin_uifilter(gcbf);' 'tag' 'filter1' } ...
    { 'style', 'popupmenu', 'string', filterList2 'callback' 'plugin_uifilter(gcbf);' 'tag' 'filter2' } ...
    { 'style', 'pushbutton', 'string', ['<html><img width=17 height=16 src="' search_icon_path '"> &nbsp; Search</html>'] 'callback' 'plugin_search(gcbf);' 'tag' 'search' 'tooltipstring' 'Enter search term' } ...
    { 'style', 'listbox', 'string', { plugin.text } 'callback' 'plugin_uiupdate(gcbf);' 'Min', 0, 'Max', 2, 'value' [] 'tag', 'pluginlist' 'fontsize', 16, 'tooltipstring', [ 'Bold plugins are installed.' 10 'Red plugins need updating.' 10 '(Wong font size? Change it in plugin_menu.m)' ] } ...
    { 'style', 'pushbutton', 'string', [ 'Rate this plugin' ] 'tag' 'rating' } ...
    { 'style', 'pushbutton', 'string', [ 'Web documentation' ] 'tag' 'documentation' } ...
    { 'style', 'pushbutton', 'string', [ 'Upload new plugin' ] 'tag' 'upload' 'callback' [ 'web(''http://sccn.ucsd.edu/eeglab/plugin_uploader/upload_form.php'', ''-browser'');' ]} ...
    { 'style', 'text', 'string', 'Tags:' 'fontweight' 'bold' } ...
    { 'style', 'text', 'string', 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 'tag' 'tags' } ...
    { 'style', 'text', 'string', 'Status:' 'fontweight' 'bold' } ...
    { 'style', 'text', 'string', 'Installed' 'tag' 'status'} ...
    { 'style', 'text', 'string', 'Size:' 'fontweight' 'bold' } ...
    { 'style', 'text', 'string', 'Size is not large' 'tag' 'size'} ...
    { 'style', 'text', 'string', 'Description of the plugin:' 'fontweight' 'bold' } ...
    { 'style', 'text', 'string', [ 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 10 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 10 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 10 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' ] 'tag' 'description' } ...
    {} ...
    {} ...
    { 'Style', 'pushbutton', 'string', 'Cancel', 'tag' 'cc' 'callback', 'close gcbf' } ...
    { 'Style', 'pushbutton', 'tag', 'rmbut', 'string', 'Remove', 'callback', 'set(findobj(gcbf, ''tag'', ''install''), ''userdata'', ''remove'');' 'enable' 'off ' } ...
    { 'Style', 'pushbutton', 'tag', 'install', 'string', 'Install/Update', 'callback', 'set(gcbo, ''userdata'', ''install'');' 'userdata' 'test' 'enable' 'off' } ...
    };

usrDat.allplugins = plugin;
usrDat.selectedplugins = plugin;
usrDat.selection = [];
fig = figure('visible', 'off');
supergui('fig', fig, 'uilist', uilist, 'geomhoriz', {[0.9 0.5 0.5 0.5] 1 [1 1 1] [0.2 1] [0.2 1] [0.2 1] 1 1 1 [0.43 0.37 0.4 0.5]}, 'geomvert', [1 10 1 1 1 1 1 2.5 1 1], 'userdata', usrDat);
%pos = get(fig, 'position');
%set(fig, 'position', [pos(1) pos(2) pos(3)/841*200 pos(4) ]);

% Remove text
set(findobj(fig, 'tag', 'tags'), 'string', '');
set(findobj(fig, 'tag', 'status'), 'string', '');
set(findobj(fig, 'tag', 'size'), 'string', '');
set(findobj(fig, 'tag', 'description'), 'string', 'Click on a plugin to show its description');

set(fig, 'visible', 'on');
waitfor( findobj('parent', fig, 'tag', 'install'), 'userdata');

% Cancel
if ~(ishandle(fig)), return; end % Check if figure still exist

% Check if install or remove
removeOrInstall = get(findobj('parent', fig, 'tag', 'install'), 'userdata');
if isempty(removeOrInstall)
    removeOrInstall = 'remove';
    ButtonName = questdlg2( [ 'Are you sure you want to remove the selected plugins?' 10 'All modification you might have made to the code will be lost.'], ...
        'Removal confirmation',  'No', 'Yes', 'No' );
    if strcmpi(ButtonName, 'No')
        close(fig);
        return
    end
end 
    
usrDat = get(fig, 'userdata');
close(fig);

% set flag in plugins
plugin = usrDat.selectedplugins;
plugin(1).install = [];
plugin(1).remove  = [];
for iSelect = 1:length(usrDat.selection)
    pInd = usrDat.selection(iSelect);
    if ~plugin(pInd).installed
        plugin(pInd).install = true;
    else
        % install or update based on which button was pressed
        if strcmpi(removeOrInstall, 'install')
            if plugin(pInd).installorupdate
                plugin(pInd).install = true;
            else
                fprintf('Skipping install for plugin %s as there is no update available\n', plugin(pInd).name);
            end
        else
            plugin(pInd).remove = true;
        end
    end
end

% install plugins
% ---------------
firstPlugin = 1;
for iRow = 1:length(plugin)
    if ~isempty(plugin(iRow).install) && plugin(iRow).install
        restartEeglabFlag = true;
        if ~firstPlugin, disp('---------------------------------'); end; firstPlugin = 0;
        
        if plugin(iRow).installed
            fprintf('Updating extension %s\n', plugin(iRow).name);
            if plugin_install(plugin(iRow).zip, plugin(iRow).name, plugin(iRow).version, plugin(iRow).size) ~= -1
                plugin_remove(plugin(iRow).foldername);
            end
        else
            fprintf('Installing extension %s\n', plugin(iRow).name);
            plugin_install(plugin(iRow).zip, plugin(iRow).name, plugin(iRow).version, plugin(iRow).size);
        end
    elseif ~isempty(plugin(iRow).remove) && plugin(iRow).remove
        if ~firstPlugin, disp('---------------------------------'); end; firstPlugin = 0; 
        restartEeglabFlag = true;
        
        fprintf('Removing extension %s\n', plugin(iRow).name);
        plugin_remove(plugin(iRow).foldername);
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
