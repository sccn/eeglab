function eeglab_update(vers)

% return at once if users said to not show this interface again
% -------------------------------------------------------------
eeglab_options;
if ~option_updateeeglab
    return
end

% try saving path
% ---------------
cb_notes = '';
cb_install = 'set(gcbf, ''userdata'', ''install'');';
cb_custom  = 'set(gcbf, ''userdata'', ''custom'');';
cb_ignore  = 'set(gcbf, ''userdata'', ''ignore'');';
% cb_install = 'set(findobj(gcbf, ''tag'', ''eeglab''), ''userdata'', ''install'');';
% cb_custom  = 'set(findobj(gcbf, ''tag'', ''eeglab''), ''userdata'', ''custom'');';
% cb_ignore  = 'set(findobj(gcbf, ''tag'', ''eeglab''), ''userdata'', ''ignore'');';
uilist   = { ...
    { 'style' 'text' 'String' 'EEGLAB xxx now available' 'fontweight' 'bold' } ...
    { 'style' 'pushbutton' 'string' 'See release notes' 'callback' cb_notes 'tag' 'eeglab' 'userdata' []} ...
    {} ...
    {} { 'style' 'pushbutton' 'String' 'Automatic install' 'callback' cb_install } {} ...
    {} { 'style' 'pushbutton' 'String' 'Custom install'    'callback' cb_custom  } {} ...
    {} { 'style' 'pushbutton' 'String' 'Ignore for now'    'callback' cb_ignore  } {} ...
    {} ...
    { 'style' 'text' 'String' [ 'Note: requires 120MB of free space plus the space ' 10 'for copying plugins from the current version' ] 'tag' 'hidemsg' } ...
    { 'style' 'checkbox' 'String' 'Do not show this message again' 'tag' 'hidemsg' } ...
    };
geom     = { [1.3 1] [1] [0.5 1 0.5] [0.5 1 0.5] [0.5 1 0.5] [1] [1] [1] };
geomvert = [ 1     1   1           1           1           1     1.5 1   ];
res = supergui( 'geomhoriz', geom, 'geomvert', geomvert, 'uilist', uilist, 'title', 'Update EEGLAB -- eeglab_update()');
set(gcf, 'userdata', 'wait');
waitfor( gcf, 'userdata');

% hide interface next time
% ------------------------
doNotShow = get(findobj(gcf, 'tag', 'hidemsg'), 'value')
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

% check for duplicate versions of EEGLAB (code from main EEGLAB code)
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

zipfilelink = 'http://sccn.ucsd.edu/eeglab/currentversion/eeglab_current.zip';

eeglabPath = fullfile( fileparts(fileparts(which('eeglab.m'))), [ 'eeglab2019.1' ]);

if strcmpi(response, 'custom')
    cb_folder = 'tmpfolder = uigetdir; if ~isequal(tmpfolder, 0), set(findobj(gcbf, ''tag'', ''folder''), ''string'', tmpfolder); end; clear tmpfolder;';
    cb_path   = 'if get(gcbo, ''value'') == 0, warndlg2([ ''This means that you will need to update the paths manually'' ], ''Path warning''); end';
    cbradio   = 'set(findobj(gcbf, ''userdata'', ''radio''), ''value'', 0); set(gcbo, ''value'', 1);';

    uilist   = { { 'style' 'text' 'String' 'EEGLAB update interface' 'fontweight' 'bold' } ...
        { 'style' 'text' 'String' 'Install folder' } ...
        { 'style' 'edit' 'string' eeglabPath 'tag' 'folder' 'horizontalalignment' 'right' } ...
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
        {} { 'style' 'checkbox' 'string' 'Copy EEGLAB preferences'             'value' 1 'tag' 'copyrefs'   } ...
        };

    geom     = { [1] [1 2 0.4] [1] [1] [0.1 1] [0.1 1] [0.1 1] [1] [1] [0.1 1] [0.1 1] [0.1 1] };
    geomvert = [ 1   1         0.5 1   1       1       1       0.5 1   1       1       1       ];
    [ res2, ~, ~, res] = inputgui( 'geometry', geom, 'geomvert', geomvert, 'uilist', uilist, ...
        'helpcom', 'pophelp(''eeglab_update'')', 'title', 'Update EEGLAB -- eeglab_update()');
    if isempty(res2), return; end
else
    res.folder      = eeglabPath;
    res.oldnothing  = false;
    res.oldmove     = true;
    res.olddelete   = false;
    res.updatepath  = true;
    res.copyplugins = true;
    res.copyrefs    = true;
end
[parentPath,eeglabFolder] = fileparts(res.folder);

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
    plugin_urlwrite( zipfilelink, fullfile(parentPath, zipfilelink));
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
unzip(fullfile(parentPath, zipfile), eeglabFolder);

disp('Cleaning EEGLAB zip file...');
delete(fullfile(parentPath, zipfile));

% seeing what is in the plugin and moving files if necessary
% ----------------------------------------------------------
eeglabContent = dir(eeglabFolder);
if length(eeglabContent) < 3
    for index = 1:length(eeglabContent)
        if ~strcmpi(eeglabContent(index).name, '.') && ~strcmpi(eeglabContent(index).name, '..')
            fullFolder = fullfile(eeglabFolder, eeglabContent(index).name);
            if exist(fullFolder) == 7 % folder detected
                % move files from folder
                movefile(fullfile(fullFolder, '*'), eeglabFolder);
                rmdir(fullFolder, 's');
            end
        end
    end
    fprintf('EEGLAB is now installed\n');
end

% update the paths
% ----------------
if res.pdatepath
    p = path;
    p = textscan(p, '%s', 'delimiter', ':');
    rmInd = [];
    for iPath = 1:length(p)
        if length(p{iPath}) >= length(eeglabpath) && strcmpi(eeglabpath, p{iPath}(1:length(eeglabpath)))
            rmInd = [ rmInd iPath ];
        end
    end
    p(rmInd) = [];
    path(p);
    addpath(eeglabFolder);
    savepath;
    disp('Paths updated and saved successfully...');
end

% Copying plugins
% ---------------


% renaming old EEGLAB folder
% --------------------------
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
        delRes = rmdir(eeglabpath);
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

function res = mywhich(varargin);
try
    res = which(varargin{:});
catch
    fprintf('Warning: permission error accesssing %s\n', varargin{1});
end