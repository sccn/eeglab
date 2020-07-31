% eeglab_update() - assess if EEGLAB new updating and ask user to install
%                   new version
%
% Usage:
%     >> eeglab_update
%
% Authors: Arnaud Delorme, SCCN/INC/UCSD, 2019

% Copyright (C) 2019 Arnaud Delorme, SCCN/INC/UCSD
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

function eeglab_update(varargin)

% return at once if users said not to check version
% -------------------------------------------------
disp('EEGLAB update menu disabled');
return;
eeglab_options;
if ~option_checkversion && nargin == 0
    return
end

%% automatic updater
try
    [~, eeglabVersionNumber, currentReleaseDateString] = eeg_getversion;
    if isempty(eeglabVersionNumber)
        eeglabVersionNumber = 'dev';
    end
    eeglabUpdater = up.updater(eeglabVersionNumber, 'http://sccn.ucsd.edu/eeglab/updater/latest_version.php', 'EEGLAB', currentReleaseDateString);
    
    % place it in the base workspace.
    assignin('base', 'eeglabUpdater', eeglabUpdater);
    
    eeglabUpdater.checkForNewVersion({'eeglab_event' 'setup'});
    if strcmpi(eeglabVersionNumber, 'dev')
        return;
    end
    newMajorRevision = 0;
    
    % eeglab new version
    eeglabv = num2str(eeglabUpdater.latestVersionNumber);
    posperiod = find(eeglabv == '.');
    if isempty(posperiod), posperiod = length(eeglabv)+1; eeglabv = [ eeglabv '.0' ]; end
    if length(eeglabv(posperiod+1:end)) < 2, eeglabv = [ eeglabv '0' ]; end
    %if length(eeglabv(posperiod+1:end)) < 3, eeglabv = [ eeglabv '0' ]; end
    eeglabv = [ eeglabv(1:posperiod+1) '.' eeglabv(posperiod+2) ]; %'.' eeglabv(posperiod+3) ];
    if strcmpi(eeglabv(end-1:end), '.0'), eeglabv(end-1:end) = []; end
    
    msg = '';
    if ~isempty(eeglabUpdater.newMajorRevision)
        msg = sprintf('A new major version of EEGLAB (EEGLAB%s) is now <a href="http://sccn.ucsd.edu/eeglab/">available</a>.', eeglabUpdater.newMajorRevision);
        fprintf('\n%s\n', msg);
        newMajorRevision = 1;
    end
    if eeglabUpdater.newerVersionIsAvailable
        
        stateWarning = warning('query', 'backtrace');
        warning('off', 'backtrace');
        if newMajorRevision
            fprintf('\n');
            msg = sprintf(['\nA critical revision of EEGLAB%d (%s) is also available <a href="%s">here</a>\n' ...
                eeglabUpdater.releaseNotes ' See <a href="matlab: web(''%s'', ''-browser'')">Release notes</a> for more informations\n' ...
                'You may disable this message in the File>Preferences menu but will miss critical updates.\n' ], ...
                floor(eeglabVersionNumber), eeglabv, eeglabUpdater.downloadUrl, eeglabUpdater.releaseNotesUrl);
            if nargin == 0, warning( msg ); end
        else
            msg =  sprintf(['\nA newer version of EEGLAB (%s) is available <a href="%s">here</a>\n' ...
                eeglabUpdater.releaseNotes ' See <a href="matlab: web(''%s'', ''-browser'')">Release notes</a> for more informations.\n' ...
                'You may disable this message in the File>Preferences menu but will miss critical updates.\n' ], ...
                eeglabv, eeglabUpdater.downloadUrl, eeglabUpdater.releaseNotesUrl);
            if nargin == 0, warning( msg ); end
        end
        warning(stateWarning.state, 'backtrace');
        
    elseif isempty(eeglabUpdater.lastTimeChecked)
        msg = [ 'Could not check for the latest EEGLAB version (internet may be disconnected).' 10 ...
                'To prevent long startup time, disable checking for new EEGLAB version (FIle>Preferences).' ];
        fprintf('%s\n', msg);
    else
        if ~newMajorRevision
            msg = 'You are using the latest version of EEGLAB.';
            fprintf('%s\n', msg);
        else
            msg = sprintf('You are currently using the latest revision of EEGLAB%d (no critical update available).', floor(eeglabVersionNumber));
            fprintf('%s\n', msg);
        end
    end
catch
    msg = 'Updater could not be initialized';
    fprintf('%s\n', msg);
end

% return at once if users said to not show this interface again or no new version available
% -----------------------------------------------------------------------------------------
eeglab_options;
if ~option_updateeeglab || ~exist('eeglabUpdater', 'var') || isempty(eeglabUpdater.newerVersionIsAvailable) || ~eeglabUpdater.newerVersionIsAvailable
    if nargin > 0
        warndlg2(msg);
    end
    return
end

% try saving path
% ---------------
cb_notes = [ 'web(''' eeglabUpdater.releaseNotesUrl ''', ''-browser'');' ];
cb_install = 'set(gcbf, ''userdata'', ''install'');';
cb_custom  = 'set(gcbf, ''userdata'', ''custom'');';
cb_ignore  = 'set(gcbf, ''userdata'', ''ignore'');';
% cb_install = 'set(findobj(gcbf, ''tag'', ''eeglab''), ''userdata'', ''install'');';
% cb_custom  = 'set(findobj(gcbf, ''tag'', ''eeglab''), ''userdata'', ''custom'');';
% cb_ignore  = 'set(findobj(gcbf, ''tag'', ''eeglab''), ''userdata'', ''ignore'');';
uilist   = { ...
    { 'style' 'text' 'String' [ 'EEGLAB ' eeglabv ' now available' ] 'fontweight' 'bold' } ...
    { 'style' 'pushbutton' 'string' 'See release notes' 'callback' cb_notes 'tag' 'eeglab' 'userdata' []} ...
    {} ...
    {} { 'style' 'pushbutton' 'String' 'Automatic install' 'callback' cb_install } {} ...
    {} { 'style' 'pushbutton' 'String' 'Custom install'    'callback' cb_custom  } {} ...
    {} { 'style' 'pushbutton' 'String' 'Ignore for now'    'callback' cb_ignore  } {} ...
    {} ...
    { 'style' 'text' 'String' [ 'Note: requires 120MB of free space plus the space ' 10 'for copying plugins from the current version' ] } ...
    { 'style' 'checkbox' 'String' 'Do not show this message again' 'tag' 'hidemsg' } ...
    };
geom     = { [1.3 1] [1] [0.5 1 0.5] [0.5 1 0.5] [0.5 1 0.5] [1] [1] [1] };
geomvert = [ 1     1   1           1           1           1     1.5 1   ];
res = supergui( 'geomhoriz', geom, 'geomvert', geomvert, 'uilist', uilist, 'title', 'Update EEGLAB -- eeglab_update()');
set(gcf, 'userdata', 'wait');
waitfor( gcf, 'userdata');

% hide interface next time
% ------------------------
doNotShow = get(findobj(gcf, 'tag', 'hidemsg'), 'value');
if doNotShow
    pop_editoptions( 'option_updateeeglab', 0);
end

% process GUI output
% ------------------
response = get(gcf, 'userdata');
close(gcf);
if strcmpi(response, 'ignore'), return; end

% try saving path
% ---------------
try
    savepath;
catch
    questdlg2( [ 'EEGLAB cannot modify and save the Matlab path file.' 10 ...
        'Although EEGLAB could still be updated, EEGLAB will not' 10 ...
        'be able to set paths in a way that is persistent after' 10 ...
        'you close Matlab. We therefore recommend that you abord and' 10 ...
        'update EEGLAB manually by downloading the EEGLAB zip file' 10 ...
        'online, uncompress it on your computer and modify and save' 10 ...
        'the Matlab paths manually.' ], 'Install warning message', 'Continue', 'Abord', 'Abord');
end

% check for duplicate versions of EEGLAB (code copied from eeglab.m)
% --------------------------------------
eeglabpath = mywhich('eeglab.m');
eeglabpath = eeglabpath(1:end-length('eeglab.m'));
if nargin < 1
    eeglabpath2 = '';
    if strcmpi(eeglabpath, pwd) || strcmpi(eeglabpath(1:end-1), pwd)
        cd('functions');
        warning('off', 'MATLAB:rmpath:DirNotFound');
        rmpath(eeglabpath);
        warning('on', 'MATLAB:rmpath:DirNotFound');
        eeglabpath2 = mywhich('eeglab.m');
        cd('..');
    else
        try, rmpath(eeglabpath); catch, end
        eeglabpath2 = mywhich('eeglab.m');
    end
    if ~isempty(eeglabpath2)
        %evalin('base', 'clear classes updater;'); % this clears all the variables
        eeglabpath2 = eeglabpath2(1:end-length('eeglab.m'));
        tmpWarning = warning('backtrace');
        warning backtrace off;
        warndlg2( [ 'There are at least two versions of EEGLAB in your path' 10 ...
            sprintf('One is at %s', eeglabpath) 10 ...
            sprintf('The other one is at %s', eeglabpath2) 10 ...
            'You must at least removeone version from the Matlab path' 10 ...
            'before you can install a new version of EEGLAB. Abording instalation.' ] );
        return
    end
end

zipfilelink = 'https://sccn.ucsd.edu/eeglab/currentversion/eeglab_current.zip';
[~,zipfile,zipext] = fileparts(eeglabUpdater.downloadUrl);
zipfile = [ zipfile zipext ];

eeglabNewPath = fullfile( fileparts(fileparts(which('eeglab.m'))), [ 'eeglab' eeglabv ]);

if strcmpi(response, 'custom')
    cb_folder = 'tmpfolder = uigetdir; if ~isequal(tmpfolder, 0), set(findobj(gcbf, ''tag'', ''folder''), ''string'', tmpfolder); end; clear tmpfolder;';
    cb_path   = 'if get(gcbo, ''value'') == 0, warndlg2([ ''This means that you will need to update the paths manually'' ], ''Path warning''); end';
    cbradio   = 'set(findobj(gcbf, ''userdata'', ''radio''), ''value'', 0); set(gcbo, ''value'', 1);';
    
    uilist   = { { 'style' 'text' 'String' 'EEGLAB update interface' 'fontweight' 'bold' } ...
        { 'style' 'text' 'String' 'Install folder' } ...
        { 'style' 'edit' 'string' eeglabNewPath 'tag' 'folder' 'horizontalalignment' 'right' } ...
        { 'style' 'pushbutton' 'string' '...' 'callback' cb_folder } ...
        {} ...
        { 'style' 'text' 'String' 'What to do with current version main EEGLAB folder?' } ...
        {} { 'style' 'checkbox' 'String' 'Leave it as is (do nothing)' 'tag' 'oldnothing'  'userdata' 'radio' 'callback' cbradio } ...
        {} { 'style' 'checkbox' 'String' 'Rename with postfix _old (recommended)' 'value' 1 'tag' 'oldmove' 'userdata' 'radio' 'callback' cbradio  } ...
        {} { 'style' 'checkbox' 'String' 'Delete it and delete all files it contains' 'tag' 'olddelete' 'userdata' 'radio' 'callback' cbradio } ...
        {} ...
        { 'style' 'text' 'String' 'Other options' } ...
        {} { 'style' 'checkbox' 'string' 'Update and resave Matlab paths' 'tag' 'updatepath' 'value' 1 'callback' cb_path } ...
        {} { 'style' 'checkbox' 'string' 'Copy current EEGLAB version plugins' 'value' 1 'tag' 'copyplugins' } ...
        {} { 'style' 'checkbox' 'string' 'Copy EEGLAB preferences'             'value' 1 'tag' 'copyprefs'   } ...
        };
    
    geom     = { [1] [1 2 0.4] [1] [1] [0.1 1] [0.1 1] [0.1 1] [1] [1] [0.1 1] [0.1 1] [0.1 1] };
    geomvert = [ 1   1         0.5 1   1       1       1       0.5 1   1       1       1       ];
    [ res2, ~, ~, res] = inputgui( 'geometry', geom, 'geomvert', geomvert, 'uilist', uilist, ...
        'helpcom', 'pophelp(''eeglab_update'')', 'title', 'Update EEGLAB -- eeglab_update()');
    if isempty(res2), return; end
else
    res.folder      = eeglabNewPath;
    res.oldnothing  = false;
    res.oldmove     = true;
    res.olddelete   = false;
    res.updatepath  = true;
    res.copyplugins = true;
    res.copyprefs   = true;
end
[parentPath,eeglabFolder,ext] = fileparts(res.folder);
eeglabFolder = [ eeglabFolder ext ];

% check if target folder can be writen into
% -----------------------------------------
if ~exist(res.folder)
    try
        createDir = mkdir(res.folder);
    catch
        createDir = 0;
    end
    if ~createDir
        msg = [ 'Parent folder of EEGLAB folder is not writable, select another location. Operation aborded.' ];
        warndlg2(msg);
        return;
    end
else
    msg = [ 'EEGLAB target folder already exist ' res.folder 10 'Remove it manually or select another folder for installing EEGLAB.' ];
    warndlg2(msg);
    return;
end

% start downloading EEGLAB
% ------------------------
try
    disp('Downloading, hang in there...');
    plugin_urlwrite( zipfilelink, fullfile(parentPath, zipfile));
catch
    msg = [ 'Could not download EEGLAB. Server might be' 10 ...
        'unavailable, or your internet might be too slow or down.' 10 ...
        'Try again just in case.' ];
    warndlg2(msg);
    return;
end

% unzip EEGLAB
% ------------
disp([ 'Unzipping EEGLAB... ']);
unzip(fullfile(parentPath, zipfile), fullfile(parentPath, eeglabFolder));

disp('Cleaning EEGLAB zip file...');
delete(fullfile(parentPath, zipfile));

% seeing what is in the plugin and moving files if necessary
% ----------------------------------------------------------
eeglabContent = dir(fullfile(parentPath, eeglabFolder));
if length(eeglabContent) < 5
    for index = 1:length(eeglabContent)
        if ~strcmpi(eeglabContent(index).name, '.') && ~strcmpi(eeglabContent(index).name, '..')
            fullFolder = fullfile(parentPath, eeglabFolder, eeglabContent(index).name);
            if exist(fullFolder) == 7 % folder detected
                % move files from folder
                movefile(fullfile(fullFolder, '*'), fullfile(parentPath, eeglabFolder));
                rmdir(fullFolder, 's');
            end
        end
    end
    fprintf('EEGLAB is now installed\n');
end

% copying eeg_options but only if the one being used is in the EEGLAB current folder
% ----------------------------------------------------------------------------------
if res.copyprefs
    eeglab_options;
    if length(option_file) > length(eeglabpath) && strcmpi(eeglabpath, option_file(1:length(eeglabpath)))
        copyfile(option_file, fullfile(parentPath, eeglabFolder, 'functions', 'adminfunc', 'eeg_options.m'));
    end
end

% update the paths
% ----------------
if res.updatepath
    p = path;
    p = textscan(p, '%s', 'delimiter', ':');
    p = p{1};
    rmInd = [];
    for iPath = 1:length(p)
        if length(p{iPath}) >= length(eeglabpath) && strcmpi(eeglabpath, p{iPath}(1:length(eeglabpath)))
            rmInd = [ rmInd iPath ];
        end
    end
    p(rmInd) = [];
    p(:,2) = { ';' };
    p = p';
    p = strcat(p{:});
    path(p);
    addpath(fullfile(parentPath, eeglabFolder));
    savepath;
    disp('Paths updated and saved successfully...');
end

% copying plugins
% ---------------
if res.copyplugins
    pluginOri  = fullfile(eeglabpath, 'plugins');
    pluginDest = fullfile(parentPath, eeglabFolder, 'plugins');
    
    pluginOriList = dir(pluginOri);
    for iPlugin = 1:length(pluginOriList)
        if contains(lower(pluginOriList(iPlugin).name), 'dipfit') && ...
                contains(lower(pluginOriList(iPlugin).name), 'firfilt') && ...
                contains(lower(pluginOriList(iPlugin).name), 'iclabel') && ...
                contains(lower(pluginOriList(iPlugin).name), 'clean_rawdata')
            copyfile(fullfile( pluginOri, pluginOriList(iPlugin).name), pluginDest);
        end
    end
end

% renaming old EEGLAB folder
% --------------------------
cd(fullfile(parentPath, eeglabFolder));
if res.oldmove
    try
        moveRes = movefile(eeglabpath, [ eeglabpath '_old' ]);
    catch
        moveRes = 0;
    end
    if ~moveRes
        disp('Could not rename the old EEGLAB folder. Make sure you have permission to do so.');
        disp('This does not affect EEGLAB installation which is now complete.');
    else
        disp('Old EEGLAB folder renamed successfully...');
    end
end
if res.olddelete
    try
        delRes = rmdir(eeglabpath, 's');
    catch
        delRes = 0;
    end
    if ~delRes
        disp('Could not delete the old EEGLAB folder. Make sure you have permission to do so.');
        disp('This does not affect EEGLAB installation which is now complete.');
    else
        disp('Old EEGLAB folder deleted successfully...');
    end
end
evalin('base', 'eeglab');

function res = mywhich(varargin);
try
    res = which(varargin{:});
catch
    fprintf('Warning: permission error accesssing %s\n', varargin{1});
end