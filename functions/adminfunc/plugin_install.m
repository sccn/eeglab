function result = plugin_install(zipfilelink, name, version);

    result = 1;

    % get plugin path
    % ---------------
    version(find(version == '.'))  = '_';
    generalPluginPath = fullfile(fileparts(which('eeglab.m')), 'plugins');
    newPluginPath     = fullfile(generalPluginPath, [ name version ]);

    % download plugin
    % ---------------
    disp([ 'Downloading ' zipfilelink ]);
    [tmp zipfile ext] = fileparts(zipfilelink);
    zipfile = [ zipfile ext ];
    try
        urlwrite( zipfilelink, fullfile(generalPluginPath, zipfile));
    catch,
        warndlg2( [ 'Could download ' zipfile ' in plugin folder.' 10 'Host site might be unavailable or you do not have' 10 'permission to write in the EEGLAB plugin folder' ]);
        result = -1;
        return;
    end;

    % unzip plugin
    % ------------
    if ~exist(newPluginPath)
        mkdir(newPluginPath);
    else
        warndlg2( [ 'Plugin folder already exist ' newPluginPath 10 'Remove it manually before installing plugin' ]);
        result = -1;
        return;
    end;
    disp([ 'Unzipping plugin file ' zipfile ]);
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
    fprintf('Plugin %s version %s now installed\n', name, version); 
    