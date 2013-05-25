function plugin_reactivate(foldername);

    % get plugin path
    % ---------------
    fullpluginfolder = fullfile(fileparts(which('eeglab.m')), 'deactivatedplugins', foldername);
    if ~exist(fullpluginfolder)
        error([ 'Could not find folder ' foldername 'in deactivatedplugins folder' ]);
    end;

    disp([ 'Moving plugin ' foldername ' to plugins folder' ]);
    fulldeactivatedpluginfolder = fullfile(fileparts(which('eeglab.m')), 'plugins');
    movefile(fullpluginfolder, fulldeactivatedpluginfolder);
    