% eeg_retrieve() - Retrieve an EEG dataset from the variable
%                  containing all datasets (standard: ALLEEG).
%
% Usage: >> EEG = eeg_retrieve( ALLEEG, index );
%
% Inputs:
%   ALLEEG     - variable containing all datasets
%   index      - index of the dataset to retrieve
%
% Outputs:
%   EEG        - output dataset. The workspace variable EEG is also updated 
%   ALLEEG     - updated ALLEEG structure
%   CURRENTSET - workspace variable index of the current dataset
%
% Note: The function performs >> EEG = ALLEEG(index);
%       It also runs eeg_checkset() on it.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeg_store(), eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [EEG, ALLEEG, CURRENTSET] = eeg_retrieve( ALLEEG, CURRENTSET);

if nargin < 2
	help eeg_retrieve;
	return;
end;	

%try
    eeglab_options;

    try,   % check whether recent changes to this dataset have been saved or not
           %--------------------------------------------------------------------
           tmpsaved = { ALLEEG.saved };
           tmpsaved = tmpsaved(CURRENTSET);
    catch, tmpsaved = 'no';
    end;

    if length(CURRENTSET) > 1 & option_storedisk
        [ EEG tmpcom ] = eeg_checkset(ALLEEG(CURRENTSET)); % do not load data if several datasets
        if length(CURRENTSET) ~= length(ALLEEG)
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        else
            ALLEEG = EEG;
        end;
    else
        if CURRENTSET ~= 0
            [ EEG tmpcom ] = eeg_checkset(ALLEEG(CURRENTSET), 'loaddata');
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        else
            EEG = eeg_emptyset; % empty dataset
            return;
        end;
    end;
    
    % retain saved field
    % ------------------
    for index = 1:length(CURRENTSET)
        ALLEEG(CURRENTSET(index)).saved = tmpsaved{index};
        EEG(index).saved                = tmpsaved{index};
    end;
    
%catch
%	fprintf('Warning: cannot retrieve dataset with index %d\n', CURRENTSET); 
%	return;
%end;

return;

