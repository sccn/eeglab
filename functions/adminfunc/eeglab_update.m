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
eeglab_options;
if nargin == 0
    [~,eeglabVersionUpdate] = plugin_getweb('update', []);
else
    eeglabVersionUpdate = varargin{1};
    if ~option_checkversion && isempty(eeglabVersionUpdate)
        return;
    end
end

% Cannot check version for some reason
if isempty(eeglabVersionUpdate)
    msg = [ 'Could not check for the latest EEGLAB version (internet may be disconnected).' 10 ...
        'To prevent long startup time, disable checking for new EEGLAB versions (File > Preferences).' ];
    stateWarning = warning('query', 'backtrace');
    warning('off', 'backtrace');
    warning(msg);
    warning(stateWarning.state, 'backtrace');
    return
end
if length(eeglabVersionUpdate) > 1
    eeglabVersionUpdate = eeglabVersionUpdate(1);
end

%% automatic update
eeglabVersionNumber = eeg_getversion;
if isempty(eeglabVersionNumber)
    eeglabVersionNumber = 'dev';
end

if ~isequal(eeglabVersionUpdate.version, eeglabVersionNumber)
    stateWarning = warning('query', 'backtrace');
    warning('off', 'backtrace');
    msg = sprintf(['\nA %s revision of EEGLAB (v%s) is available <a href="matlab:eeglab_update;">HERE</a>.\n%s\n' ...
        'See <a href="matlab: web(''%s'', ''-browser'')">Release notes</a> for more information\n' ...
        'You may disable this message in the File > Preferences menu.\n' ], ...
        fastif(eeglabVersionUpdate.critical, 'CRITICAL', 'newer'), eeglabVersionUpdate.version, ...
        eeglabVersionUpdate.releasenotes, eeglabVersionUpdate.webdoc);
    if nargin > 0, warning( msg ); end
    warning(stateWarning.state, 'backtrace');
else
    msg = 'You are using the latest version of EEGLAB.';
    if nargin > 0
        fprintf('%s\n', msg);
    else
        warndlg2(msg);
    end
    return
end

% return at once if users said to not show this interface again or no new version available
% -----------------------------------------------------------------------------------------
eeglab_options;
if nargin > 0
    return
end    

% try saving path
% ---------------
cb_notes = [ 'web(''' eeglabVersionUpdate.webdoc ''', ''-browser'');' ];
cb_install = 'set(gcbf, ''userdata'', ''install'');';
cb_custom  = 'set(gcbf, ''userdata'', ''custom'');';
cb_ignore  = 'set(gcbf, ''userdata'', ''ignore'');';
% cb_install = 'set(findobj(gcbf, ''tag'', ''eeglab''), ''userdata'', ''install'');';
% cb_custom  = 'set(findobj(gcbf, ''tag'', ''eeglab''), ''userdata'', ''custom'');';
% cb_ignore  = 'set(findobj(gcbf, ''tag'', ''eeglab''), ''userdata'', ''ignore'');';
uilist   = { ...
    { 'style' 'text' 'String' [ 'EEGLAB ' eeglabVersionUpdate.version ' now available' ] 'fontweight' 'bold' } ...
    { 'style' 'pushbutton' 'string' 'See release notes' 'callback' cb_notes 'tag' 'eeglab' 'userdata' []} ...
    {} ...
    {} { 'style' 'pushbutton' 'String' 'Automatic install' 'callback' cb_install } {} ...
    {} { 'style' 'pushbutton' 'String' 'Custom install'    'callback' cb_custom  } {} ...
    {} { 'style' 'pushbutton' 'String' 'Ignore for now'    'callback' cb_ignore  } {} ...
    {} ...
    { 'style' 'text' 'String' [ 'Note: requires 120MB of free space plus the space ' 10 'for copying plugins from the current version' ] } ...
    };
geom     = { [1.3 1] [1] [0.5 1 0.5] [0.5 1 0.5] [0.5 1 0.5] [1] [1] };
geomvert = [ 1     1   1           1           1           1     1.5   ];
res = supergui( 'geomhoriz', geom, 'geomvert', geomvert, 'uilist', uilist, 'title', 'Update EEGLAB -- eeglab_update()');
set(gcf, 'userdata', 'wait');
waitfor( gcf, 'userdata');
response = get(gcf, 'userdata');
if iscell(response) % closed figure
    return;
end

% process GUI output
% ------------------
close(gcf);
if strcmpi(response, 'ignore'), return; end

% try saving path
% ---------------
try
    res = savepath;
catch
    res = 1;
end
if res
    questdlg2( [ 'EEGLAB cannot modify and save the Matlab path file.' 10 ...
        'Although EEGLAB could still be updated, EEGLAB will not' 10 ...
        'be able to set paths in a way that is persistent after' 10 ...
        'you close Matlab. We therefore recommend that you abort and' 10 ...
        'update EEGLAB manually by downloading the EEGLAB zip file' 10 ...
        'online, uncompress it on your computer and modify and save' 10 ...
        'the Matlab paths manually.' ], 'Install warning message', 'Continue', 'Abord', 'Abord');
end
clear res;

% check for duplicate versions of EEGLAB (code copied from eeglab.m)
% --------------------------------------
eeglabpath = fileparts(mywhich('eeglab.m'));
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
        eeglabpath2 = eeglabpath2(1:end-length('eeglab.m'));
        tmpWarning = warning('backtrace');
        warning backtrace off;
        warndlg2( [ 'There are at least two versions of EEGLAB in your path' 10 ...
            sprintf('One is at %s', eeglabpath) 10 ...
            sprintf('The other one is at %s', eeglabpath2) 10 ...
            'You must at least remove one version from the Matlab path' 10 ...
            'before you can install a new version of EEGLAB. Abording installation.' ] );
        return
    end
end

% force HTTPS
eeglabVersionUpdate.zip = strrep(eeglabVersionUpdate.zip, 'http://', 'https://');
[~,zipfile,zipext] = fileparts(eeglabVersionUpdate.zip);
zipfile = [ zipfile zipext ];

eeglabNewPath = fullfile( fileparts(fileparts(which('eeglab.m'))), [ 'eeglab' eeglabVersionUpdate.version ]);

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
        {} { 'style' 'popupmenu' 'String' { 'Leave it as is (do nothing)' 'Rename with postfix _old (recommended)' 'Delete it and delete all files it contains' } 'value' 2 'tag' 'popmenuchoice' } ...
        {} ...
        { 'style' 'text' 'String' 'Other options' } ...
        {} { 'style' 'checkbox' 'string' 'Update and resave Matlab paths' 'tag' 'updatepath' 'value' 1 'callback' cb_path } ...
        {} { 'style' 'checkbox' 'string' 'Copy current EEGLAB version plugins' 'value' 1 'tag' 'copyplugins' } ...
        {} { 'style' 'checkbox' 'string' 'Copy EEGLAB preferences'             'value' 1 'tag' 'copyprefs'   } ...
        };
    
    geom     = { [1] [1 2 0.4] [1] [1] [0.1 1] [1] [1] [0.1 1] [0.1 1] [0.1 1] };
    geomvert = [ 1   1         0.5 1   1       0.5 1   1       1       1       ];
    [ res2, ~, ~, res] = inputgui( 'geometry', geom, 'geomvert', geomvert, 'uilist', uilist, ...
        'helpcom', 'pophelp(''eeglab_update'')', 'title', 'Update EEGLAB -- eeglab_update()');
    if isempty(res2), return; end
    res.oldnothing  = false;
    res.oldmove     = false;
    res.olddelete   = false;
    switch res.popmenuchoice
        case 1, res.oldnothing = true;
        case 2, res.oldmove    = true;
        case 3, res.olddelete  = true;
    end
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

% check if target folder can be written into
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
    restmp = questdlg2('This will download about 80Mb with no wait bar, be patient.', 'Warning', 'Cancel', 'OK', 'OK');
    if isempty(restmp) || strcmpi(restmp, 'Cancel'), return; end
    disp('Downloading about 80Mb, hang in there...');
    plugin_urlwrite( eeglabVersionUpdate.zip, fullfile(parentPath, zipfile));
catch
    msg = [ 'Could not download EEGLAB. Server might be unavailable, or your internet might' 10 ...
        'be down or too slow. Alternatively, your version of Matlab might not support HTTPS.' 10 ...
        'Try again just in case. Otherwise, download the new zip file from the Internet.' ];
    warndlg2(msg);
    return;
end

% unzip EEGLAB
% ------------
disp([ 'Unzipping EEGLAB... ']);
unzip(fullfile(parentPath, zipfile), fullfile(parentPath, eeglabFolder));

disp('Cleaning EEGLAB zip file...');
delete(fullfile(parentPath, zipfile));

% seeing what is in the plugin folder and moving files if necessary
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
    if ispc
        p = textscan(p, '%s', 'delimiter', ';');
    else
        p = textscan(p, '%s', 'delimiter', ':');
    end
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
disp('EEGLAB require some plugins to function properly.');
disp('Updated versions of plugins Dipfit, Firfilt, ICLabel, and clean_rawdata are automatically included');
if res.copyplugins
    pluginOri  = fullfile(eeglabpath, 'plugins');
    pluginDest = fullfile(parentPath, eeglabFolder, 'plugins');
    
    pluginOriList = dir(pluginOri);
    for iPlugin = 1:length(pluginOriList)
        if ~isequal(pluginOriList(iPlugin).name, '.') && ...
                ~isequal(pluginOriList(iPlugin).name, '..') && ...
                isempty(strfind(lower(pluginOriList(iPlugin).name), 'dipfit')) && ...
                isempty(strfind(lower(pluginOriList(iPlugin).name), 'firfilt')) && ...
                isempty(strfind(lower(pluginOriList(iPlugin).name), 'iclabel')) && ...
                isempty(strfind(lower(pluginOriList(iPlugin).name), 'clean_rawdata'))
            destPath = fullfile(pluginDest, pluginOriList(iPlugin).name);
            try
                mkdir(destPath);
                copyfile(fullfile( pluginOri, pluginOriList(iPlugin).name, '*'), destPath);
                fprintf('Plugin %s copied successfully\n', pluginOriList(iPlugin).name);
            catch
                fprintf('Issue with copying plugin %s - we suggest you reinstall it from the plugin manager\n', pluginOriList(iPlugin).name);
            end
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
    fprintf('Warning: permission error accessing %s\n', varargin{1});
end
