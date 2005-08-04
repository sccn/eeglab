% eeg_options() - eeglab option script 
%
% Note: DO NOT EDIT, instead use pop_editoptions() or the menu
%       /File/Maximize memory in EEGLAB gui

% Memory options 
option_storedisk = 0 ;   % If set, keep at most one dataset in memory (allow processing dozens of datasets at a time)
option_warningstore = 1 ;   % If set, warn user before overwriting a file on disk (only relevant when option above is set)
option_savematlab = 1 ;  % If set, write data in same file as dataset (unset - 2 files - allow reading one channel at a time from disk)
option_computeica = 1 ;  % If set, precompute ICA activations (requires RAM, but faster plotting of component activation) 
% Folder options
option_rememberfolder = 1 ;  % If set, remember old folder when reading dataset
