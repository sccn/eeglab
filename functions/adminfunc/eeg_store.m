% eeg_store() - store specified EEG dataset(s) in the ALLEEG variable 
%               containing all current datasets, after first checking 
%               dataset consistency using eeg_checkset().
%
% Usage: >> [ALLEEG EEG index] = eeg_store(ALLEEG, EEG);
%        >> [ALLEEG EEG index] = eeg_store(ALLEEG, EEG, index);
%
% Inputs:
%   ALLEEG     - variable containing all current EEGLAB datasets
%   EEG        - dataset(s) to store - usually the current dataset. 
%                May also be an array of datasets; these will be 
%                checked and stored separately in ALLEEG.
%   index      - (optional), ALLEEG index (or indices) to use to store 
%                the new dataset(s). If no index is given, eeg_store() 
%                uses the lowest empty slot(s) in the ALLEEG array. 
% Outputs:
%   ALLEEG - array of all current datasets
%   EEG    - EEG dataset (after syntax checking)
%   index  - index of the new dataset
%
% Note: When 3 arguments are given, after checking the consistency of 
%       the input dataset structures this function simply performs 
%        >> ALLEEG(index) = EEG;
%
% Typical use:
%        >> [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
%  creates a new dataset in variable ALLEEG.
%        >> [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
%  overwrites the current dataset in variable ALLEEG.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab(), eeg_checkset(), eeg_retrieve()

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

% store set
% ------------------
function [ALLEEG, EEG, storeSetIndex] = eeg_store(ALLEEG, EEG, storeSetIndex, varargin)

% check parameter consistency
% ------------------------------
if nargin >= 3
    if length(EEG) ~= length(storeSetIndex) && storeSetIndex(1) ~= 0
        error('Length of input dataset structure must equal the length of the index array');
    end
end

% considering multiple datasets
% -----------------------------
if length(EEG) > 1
	TMPEEG = EEG;
    if nargin >= 3 && storeSetIndex(1) ~= 0
        for index=1:length(TMPEEG)
            EEG = TMPEEG(index);
            tmpsaved      = EEG.saved;
            if strcmpi(tmpsaved, 'justloaded'), tmpsaved = 'yes'; end
            [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, storeSetIndex(index), varargin{:});
            ALLEEG(storeSetIndex(index)).saved = tmpsaved;
            TMPEEG(index).saved                = tmpsaved;
        end
    else
        for index=1:length(TMPEEG)
            EEG = TMPEEG(index);
            tmpsaved      = EEG.saved;
            if strcmpi(tmpsaved, 'justloaded'), tmpsaved = 'yes'; end
            [ALLEEG, EEG, storeSetIndex(index)] = eeg_store(ALLEEG, EEG);
            ALLEEG(storeSetIndex(index)).saved = tmpsaved;
            TMPEEG(index).saved                = tmpsaved;
        end
    end
    EEG = TMPEEG;
	return
end

if nargin < 3
    % creating new dataset
    % -> erasing file information
    % ---------------------------
    EEG.filename = '';
    EEG.filepath = '';
    EEG.datfile  = '';
end

if isempty(varargin) % no text output and no check (study loading)
    [ EEG, com ]  = eeg_checkset(EEG);
else
    com = '';
end
if nargin > 2 
    if storeSetIndex == 0 || strcmpi(EEG.saved, 'justloaded')
        EEG.saved = 'yes'; % just loaded
    else 
        EEG.saved = 'no';
    end
elseif strcmpi(EEG.saved, 'justloaded')
    EEG.saved = 'yes';        
else
    EEG.saved = 'no';        
end
EEG = eeg_hist(EEG, com);

% find first free index
% ---------------------
findindex = 0;
if nargin < 3,             findindex = 1;
elseif storeSetIndex == 0, findindex = 1; 
end

if findindex
	i = 1;
	while (i<10000)
		try
			if isempty(ALLEEG(i).data)
				storeSetIndex = i; i = 10000;
			end
			i = i+1;	
		catch
			storeSetIndex = i; i = 10000;
		end
   end
   if isempty(varargin) % no text output and no check
       fprintf('Creating a new ALLEEG dataset %d\n', storeSetIndex);
   end
else
	if isempty(storeSetIndex) || storeSetIndex == 0
		storeSetIndex = 1;
	end
end

if ~isempty( ALLEEG )
	try
		ALLEEG(storeSetIndex) = EEG;
	catch
		allfields = fieldnames( EEG );
		for i=1:length( allfields )
			ALLEEG(storeSetIndex).(allfields{i}) = EEG.(allfields{i});
        end
        if ~isfield(EEG, 'datfile') && isfield(ALLEEG, 'datfile')
            ALLEEG(storeSetIndex).datfile = '';
        end
	end
else	
	ALLEEG = EEG;
	if storeSetIndex ~= 1
		ALLEEG(storeSetIndex+1) = EEG;
		ALLEEG(1) = ALLEEG(storeSetIndex); % empty
 		ALLEEG(storeSetIndex) = ALLEEG(storeSetIndex+1); 
 		ALLEEG = ALLEEG(1:storeSetIndex);
    end	
end	
