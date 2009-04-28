function [eegPath] = eeg_toolbox_path

% eeg_toolbox_path - locate the eeg_toolbox installation path
%
% [eegPath] = eeg_toolbox_path
%
% This function looks for the installation path of
% eeg_toolbox.m in the matlab path.  If not 
% found, update the matlab path (try addpath).
%


% $Revision: 1.1 $ $Date: 2009-04-28 22:13:53 $

% Licence:  GNU GPL, no implied or express warranties
% History:  12/2003, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eegPath = fileparts(which('eeg_toolbox'));
if isempty(eegPath),
    msg = sprintf('Cannot find eeg_toolbox on the matlab path.\nPlease use the addpath command.\n\n');
    error(msg);
else
    eegPath = strcat(eegPath,filesep);
end

return
