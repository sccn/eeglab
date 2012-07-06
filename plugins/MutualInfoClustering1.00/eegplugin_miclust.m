function eegplugin_miclust( fig, try_strings, catch_strings);
plotmenu = findobj(fig, 'tag', 'plot');
uimenu(plotmenu, 'label', 'Cluster dataset ICs', 'callback', ...
    [try_strings.check_ica '[EEG LASTCOM]= pop_miclust(EEG);' catch_strings.store_and_hist '[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);' ]);
