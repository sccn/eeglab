% eeg_epoch2continuous() - convert epoched dataset to continuous dataset
%                          with data epochs separated by boundary events.
% Usage:
%           >> EEGOUT = eeg_epoch2continuous(EEGIN);
%
% Inputs:
%   EEGIN  - a loaded epoched EEG dataset structure.
%
% Inputs:
%   EEGOUT - a continuous EEG dataset structure.
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, January, 2012

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, October 11, 2012, arno@sccn.ucsd.edu
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

function EEG = eeg_epoch2continuous(EEG)

if nargin < 1
    help eeg_epoch2continuous;
    return;
end

EEG.data = reshape(EEG.data, size(EEG.data,1), size(EEG.data,2)*size(EEG.data,3));

eeglab_options;
for index = 1:EEG.trials-1
    EEG.event(end+1).type = eeg_boundarytype(EEG);
    EEG.event(end  ).latency  = index*EEG.pnts-0.5;
    EEG.event(end  ).duration = NaN;
end

EEG.pnts   = size(EEG.data,2);
EEG.trials = 1;
if ~isempty(EEG.event) && isfield(EEG.event, 'epoch')
    EEG.event = rmfield(EEG.event, 'epoch');
end

EEG = eeg_checkset(EEG);
