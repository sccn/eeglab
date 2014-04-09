function plugin_remove(foldername);

    % get plugin path
    % ---------------
    fullpluginfolder = fullfile(fileparts(which('eeglab.m')), 'deactivatedplugins', foldername);
    if ~exist(fullpluginfolder)
        error([ 'Could not find folder ' foldername ' in deactivatedplugins folder' ]);
    end;

    disp([ 'Removing plugin folder ' foldername ]);
    try
        rmdir(fullpluginfolder, 's');
    catch
        eeglab_error;
    end;