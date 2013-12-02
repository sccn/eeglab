function result = plugin_install(zipfilelink, name, version, forceInstall);

    if nargin < 4, forceInstall = false; end;
    result = 1;

    % get plugin path
    % ---------------
    %version(find(version == '.'))  = '_';
    generalPluginPath = fullfile(fileparts(which('eeglab.m')), 'plugins');
    newPluginPath     = fullfile(generalPluginPath, [ name version ]);

    % download plugin
    % ---------------
    zipfile = 'tmp.zip';
%     [tmp zipfile ext] = fileparts(zipfilelink);
%     zipfile = [ zipfile ext ];
%     equalPos = find(zipfile == '=');
%     if ~isempty(equalPos) zipfile  = zipfile(equalPos(end)+1:end); end;
    depth = length(dbstack);
    try
        pluginSize = plugin_urlsize(zipfilelink);
        pluginSizeStr = num2str(round(pluginSize/100000)/10);
        if pluginSize > 500000 && depth > 1 && ~forceInstall
             res = questdlg2( [ 'Extension ' name ' size is ' pluginSizeStr 'MB. Are you sure' 10 ...
                                'you want to download this extension?' ], 'Warning', 'No', 'Yes', 'Yes');
             if strcmpi(res, 'no'), fprintf([ 'Skipping ' name ' extension instalation\n' ]); result = -1; return; end;               
        end;
    catch,
        warndlg2( [ 'Could not download extension. Host site might be' 10 'unavailable or you do not have permission' 10 'to write in the EEGLAB plugin folder' ]);
        result = -1;
        return;
    end;
    disp([ 'Downloading extension ' name '(' pluginSizeStr 'Mb)...' ]);
    
    try
        plugin_urlread(['http://sccn.ucsd.edu/eeglab/plugin_uploader/plugin_increment.php?plugin=' name '&version=' version ]);
        plugin_urlwrite( zipfilelink, fullfile(generalPluginPath, zipfile));
    catch,
        warndlg2( [ 'Could not download extension. Host site might be' 10 'unavailable or you do not have permission' 10 'to write in the EEGLAB plugin folder' ]);
        result = -1;
        return;
    end;

    % unzip plugin
    % ------------
    if ~exist(newPluginPath)
        mkdir(newPluginPath);
    else
        msg = [ 'Extension folder already exist ' newPluginPath 10 'Remove it manually before installing extension' ];
        if ~forceInstall
            warndlg2(msg);
        else
            disp(msg);
        end;
        result = -1;
        return;
    end;
    disp([ 'Unzipping extension file... ']);
    unzip(fullfile(generalPluginPath, zipfile), newPluginPath);
    
    disp('Cleaning up zip file...');
    delete(fullfile(generalPluginPath, zipfile));
    
    % seeing what is in the plugin and moving files if necessary
    % ----------------------------------------------------------
    pluginContent = dir(newPluginPath);
    if length(pluginContent) > 3
        return;
    end;
    for index = 1:length(pluginContent)
        if ~strcmpi(pluginContent(index).name, '.') && ~strcmpi(pluginContent(index).name, '..')
            fullFolder = fullfile(newPluginPath, pluginContent(index).name);
            if exist(fullFolder) == 7 % folder detected
                % move files from folder
                movefile(fullfile(fullFolder, '*'), newPluginPath);
                rmdir(fullFolder, 's');
            end;
        end;
    end;
    fprintf('Extension %s version %s now installed\n', name, version); 
    