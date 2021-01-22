% eeg_options() - eeglab option script 
%
% Note: DO NOT EDIT, instead use pop_editoptions() or the menu
%       /File/Maximize memory in EEGLAB gui

% STUDY and file options (set the first checkbox if you intend to work with studies) 
option_storedisk     = 0 ;  % If set, keep at most one dataset in memory. This allows processing hundreds of datasets within studies.
option_savetwofiles  = 0 ;  % If set, save not one but two files for each dataset (header and data). No longer set by default as of 2021.
option_saveversion6  = 1 ;  % If set, write files in Matlab v6.5 (fastest and Octave compatible). If not, write files in Matlab v7.3 (for files > 2Gb).
option_saveasstruct  = 1 ;  % If set, save the fields of the EEG structure as individual variables in the file (new 2021 default).
% Memory options 
option_single        = 1 ;  % If set, use single precision number (32-bit instead of 64-bit) in memory.
option_memmapdata    = 0 ;  % If set, use memory mapped array under Matlab 7.x. This may slow down some computation (beta).
option_eegobject     = 0 ;  % If set, use the EEGLAB EEG object instead of the standard EEG structure (beta). 
% ICA options 
option_computeica    = 1 ;  % If set, precompute ICA activations. This requires more RAM but allows faster plotting of component activations. 
option_scaleicarms   = 1 ;  % If set, scale ICA component activities to RMS (Root Mean Square) in microvolt (recommended).
% Folder options
option_rememberfolder = 1 ;  % If set, when browsing to open a new dataset assume the folder/directory of the previous dataset.
% Toolbox options
option_donotusetoolboxes = 0 ;  % If set, do not use Matlab additional toolboxes functions even if they are present (need to restart EEGLAB).
% EEGLAB connectivity and support
option_showadvanced      = 0 ;  % If set, show advanced options (close and reopen this GUI to effect changes)
option_showpendingplugins = 0 ;  % If set, show plugins pending approval instead of approved plugins (for developers only) 
option_allmenus          = 0 ;  % If set, show all menu items from previous EEGLAB versions. You must restart EEGLAB for this to take effect.
option_checkversion      = 1 ;  % If set, check for new version of EEGLAB and EEGLAB extensions at startup.
option_cachesize         = 500 ;  % Size of cache in Mbytes for EEGLAB STUDY cache in RAM.
