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
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

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
    end

    if length(CURRENTSET) > 1 && option_storedisk
        [ EEG tmpcom ] = eeg_checkset(ALLEEG(CURRENTSET)); % do not load data if several datasets
        if length(CURRENTSET) ~= length(ALLEEG)
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        else
            ALLEEG = EEG;
        end
    else
        if CURRENTSET ~= 0
            [ EEG tmpcom ] = eeg_checkset(ALLEEG(CURRENTSET), 'loaddata');
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        else
            EEG = eeg_emptyset; % empty dataset
            return;
        end
    end
    
    % retain saved field
    % ------------------
    for index = 1:length(CURRENTSET)
        ALLEEG(CURRENTSET(index)).saved = tmpsaved{index};
        EEG(index).saved                = tmpsaved{index};
    end
    
%catch
%	fprintf('Warning: cannot retrieve dataset with index %d\n', CURRENTSET); 
%	return;
%end

return;

