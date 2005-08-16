% eeg_options() - eeglab option script 
%
% Note: DO NOT EDIT, instead use pop_editoptions() or the menu
%       Modifed by EEGLAB gui menu item: File/Maximize memory 

% Memory options 
option_storedisk = 0 ;   % If set, keep at most one dataset in memory. This allows processing many datasets at a time.
option_savematlab = 1 ;  % If set, write data in same file as dataset; if unset, use two files per dataset
option_computeica = 1 ;  % If set, precompute ICA activations. Requires more RAM, but allows faster plotting of component activations. 
% Folder options
option_rememberfolder = 0 ;  % If set, assume the folder/directory of the previous file read when reading a new dataset; if unset, look first in the present working directory (pwd)
