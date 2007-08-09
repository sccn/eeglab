% eeg_options() - eeglab option script 
%
% Note: DO NOT EDIT, instead use pop_editoptions() or the menu
%       /File/Maximize memory in EEGLAB gui

% Memory options 
option_storedisk = 0 ;   % If set, keep at most one dataset in memory. Allows processing many datasets at a time.
option_savematlab = 1 ;  % If set, write the data and dataset structure into the same file; if unset, write two files per dataset (.dat .set).
option_computeica = 1 ;  % If set, precompute ICA activations. This requires more RAM but allows faster plotting of component activations. 
% Folder options
option_rememberfolder = 1 ;  % If set, when reading dataset assume the folder/directory of previous dataset.
% ICA options
option_scaleicarms = 0 ;  % If set, scale ICA component activities to RMS (Root Mean Square) in microvolt
