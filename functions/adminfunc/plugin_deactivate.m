function plugin_deactivate(foldername);

    % get plugin path
    % ---------------
    fullpluginfolder = fullfile(fileparts(which('eeglab.m')), 'plugins', foldername);
    if ~exist(fullpluginfolder)
        error([ 'Could not find folder ' foldername 'in extension folder' ]);
    end;

    disp('Removing extension from path');
    allPaths = path;
    indSep   = [0 find(allPaths == ':') length(allPaths)+1];
    for index = 1:length(indSep)-1
        tmpPath = allPaths(indSep(index)+1:indSep(index+1)-1);
        if ~isempty(strfind(tmpPath, fullpluginfolder))
            rmpath(tmpPath);
            fprintf('Removing path %s\n', tmpPath);
        end;
    end;

    disp([ 'Moving extension ' foldername ' to deactivatedplugins folder' ]);
    fulldeactivatedpluginfolder = fullfile(fileparts(which('eeglab.m')), 'deactivatedplugins');
    if ~exist(fulldeactivatedpluginfolder), mkdir(fulldeactivatedpluginfolder); end;
    try
        movefile(fullpluginfolder, fulldeactivatedpluginfolder);
    catch
        eeglab_error;
    end;