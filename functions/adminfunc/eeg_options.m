% eeg_options() - eeglab option script 
%
% Note: DO NOT EDIT, instead use pop_editoptions() or the menu
%       /File/Maximize memory in EEGLAB gui

option_computeica = 1 ;  % [Set|Unset] -> [Do|Do not] Precompute ICA activations (+RAM) 
option_keepdataset = 1 ; % [Set|Unset] -> [Do|Do not] Retain parent datasets (+RAM)
option_savematlab = 1 ;  % [Set|Unset] -> Store EEG.data in [the .set|a .dat] file
option_rememberfolder = 1 ;  % [Set|Unset] -> When reading EEGLAB dataset, remember old folder
option_storedisk = 1 ;   % [Set|Unset] -> [Do|Do not] store non-current dataset on disk
