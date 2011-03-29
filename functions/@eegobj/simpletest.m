% simple dataset tests
p = fileparts(which('eeglab'));
EEG = pop_loadset('filename','eeglab_data_epochs_ica.set','filepath',fullfile(p, 'sample_data'));

