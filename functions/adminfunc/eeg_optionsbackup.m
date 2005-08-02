% eeg_options() - eeglab option script 
%
% Note: DO NOT EDIT, instead use pop_editoptions() or the menu
%       /File/Maximize memory in EEGLAB gui

% Memory options
option_storedisk = 0 ;   % [Set|Unset] -> [Do|Do not] keep at most one dataset in memory (this allow processing dozens of datasets at a time)
option_savematlab = 1 ;  % [Set|Unset] -> [Do|Do not] write data in a separate file from datasets (this allow some functions to read one channel at a time from disk)
option_computeica = 1 ;  % [Set|Unset] -> [Do|Do not] Precompute ICA activations (+RAM) 
% Folder options
option_rememberfolder = 1 ;  % [Set|Unset] -> [Do|Do not] remember old folder when reading dataset
