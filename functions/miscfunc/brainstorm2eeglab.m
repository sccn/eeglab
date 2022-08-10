% brainstorm2eeglab - convert Brainstorm structures to EEGLAB dataset
%
% EEG = brainstorm2eeglab(bst);
%
% Inputs:
%    bst - Brainstorm string for a folder containing data epoch files or Brainstorm 
%          data epoch structure
%
% Limitations:
%   - Only import data epochs
%   - Cannot import linked file
%
% Output:
%    EEG     - EEGLAB structure
%
% Author: Arnaud Delorme, UCSD

% Copyright (C) Arnaud Delorme, UCSD 2022
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

function EEG = brainstorm2eeglab(bst)

if nargin < 1
    help brainstorm2eeglab;
    return;
end

if ischar(bst)
    oriBst = bst;
    bst = dir(fullfile(bst, '*.mat'));
    if isempty(bst)
        error('Files not found')
    end
else
    oriBst = '';
    % Need brainstrorm
    try
        prot = bst_get('ProtocolInfo');
    catch
        brainstorm nogui
        try
            prot = bst_get('ProtocolInfo');
        catch
            error('Brainstorm not found. Brainstorm must be installed and running.')
        end
    end
    brainstorm_path = bst_get('BrainstormDbDir');
    protocol = prot.Comment;

    for iEpoch = 1:length(bst)
        bst(iEpoch).name = bst(iEpoch).FileName; 
        bst(iEpoch).folder = fullfile(brainstorm_path, protocol, 'data');
    end
end

disp('Reading data epochs...');
allEpochs = {};
for iEpoch = 1:length(bst)
    if ~isequal(bst(iEpoch).name, 'brainstormstudy.mat')
        epochStruct = load('-mat', fullfile(bst(iEpoch).folder, bst(iEpoch).name));
        % epochStruct.Comment(1) contain the type of event
        allEpochs{end+1} = epochStruct.F;
    end
end
data = reshape([allEpochs{:}], size(allEpochs{1},1), size(allEpochs{1},2), length(allEpochs));

% create EEG structure
EEG = eeg_emptyset;
EEG.srate = mean(1./diff(epochStruct.Time));
EEG.data  = data;

% channel locations
if isempty(oriBst)
    chanFile = fullfile(bst(1).folder, bst(1).ChannelFile);
else
    chanFile = fullfile(oriBst, '..', '@default_study', 'channel.mat');    
end
if exist(chanFile, 'file')
    chans = load('-mat', chanFile);
    for iChan = 1:length(chans.Channel)
        EEG.chanlocs(iChan).labels = chans.Channel(iChan).Name;
        if ~isempty(chans.Channel(iChan).Loc)
            EEG.chanlocs(iChan).X      = chans.Channel(iChan).Loc(1);
            EEG.chanlocs(iChan).Y      = chans.Channel(iChan).Loc(2);
            EEG.chanlocs(iChan).Z      = chans.Channel(iChan).Loc(3);
        end
        EEG.chanlocs(iChan).type   = chans.Channel(iChan).Type;
    end
    EEG.chanlocs = convertlocs(EEG.chanlocs, 'cart2all');
end
EEG = eeg_checkset(EEG);
