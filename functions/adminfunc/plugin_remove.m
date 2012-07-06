function plugin_remove(foldername);

    % get plugin path
    % ---------------
    fullpluginfolder = fullfile(fileparts(which('eeglab.m')), 'plugins', foldername);
    if ~exist(fullpluginfolder)
        error([ 'Could not find folder ' foldername 'in plugin folder' ]);
    end;

    disp('Removing plugin from path');
    allPaths = path;
    indSep   = [0 find(allPaths == ':') length(allPaths)+1];
    for index = 1:length(indSep)-1
        tmpPath = allPaths(indSep(index)+1:indSep(index+1)-1);
        if ~isempty(strfind(tmpPath, fullpluginfolder))
            rmpath(tmpPath);
            fprintf('Removing path %s\n', tmpPath);
        end;
    end;
    
    disp([ 'Removing plugin folder ' foldername ]);
    rmdir(fullpluginfolder, 's');
    