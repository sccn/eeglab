% EEGLAB_NEW - script called when a new dataset is created and require 
%              input from user. Use workspace variables EEG, ALLEEG,
%              CURRENTSET, LASTCOM

functionsEEGLAB = { 'pop_chanedit' 'pop_editset' 'pop_comment' 'pop_headplot' 'pop_selectcomps' ...
              'pop_runica'   'pop_editeventfield' 'pop_editeventvals' 'pop_adjustevents' ...
              'pop_saveset'  'pop_importepoch'    'pop_importevent'   'pop_chanevent' ...
              'pop_importpres' 'pop_importerplab' 'topoplot' };

posEqual = find('=' == LASTCOM);
if ~isempty(LASTCOM) && ~isempty(EEG) && ~isempty(posEqual) && contains(LASTCOM(1:posEqual(1)), 'EEG') && ~contains(LASTCOM, functionsEEGLAB) % data assignment

    % Here the dataset was modified 
    EEG = eegh(LASTCOM, EEG); 
    [ALLEEG, EEG, CURRENTSET, LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET, 'study', ~isempty(STUDY)+0);
    eegh(LASTCOM);
    disp('Done');
    eeglab('redraw');

elseif ~isempty(LASTCOM)

    % Here the dataset was not modified (plotting only) 
    EEG = eegh(LASTCOM, EEG); 
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    eeglab('redraw');
    % even though EEG is modified by eegh call, do not keep eeg_store in history
    % so the history does not reflect modifications to EEG.history for
    % better redeability

end

clear functionsEEGLAB posEqual;