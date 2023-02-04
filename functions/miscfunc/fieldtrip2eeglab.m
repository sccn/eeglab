% fieldtrip2eeglab - convert Fieldtrip structures to EEGLAB dataset
%
% EEG = fieldtrip2eeglab(header, rawdata, evt);
% EEG = fieldtrip2eeglab(data);
%
% Inputs:
%    header   - Fieldtrip data header 
%    rawdata  - Fieldtrip raw data
%    evt      - Fieldtrip event structure (optional)
%    data     - Fieldtrip data out of ft_preprocessing. Note that this uses
%               a legacy conversion method. It is better to use 
%               fieldtrip2eeglab(data.hdr, data.trial) to use the default
%               FileIO API.
%
% Output:
%    EEG     - EEGLAB structure
%
% Author: Arnaud Delorme, UCSD

% Copyright (C) Arnaud Delorme, UCSD 2018
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

function EEG = fieldtrip2eeglab(header,data,evt)

if nargin < 3
    evt = [];
end

if nargin >= 2
    EEG = pop_fileio(header, data, evt);
else
    fprintf(2, 'fieldtrip2eeglab: Use 2-input argument, header and data to use the fileio API to convert Fieldtrip data\n')

    if isfield(header, 'hdr')
        % use the hdr field in the data
        hdr = header.hdr;
        if isfield(header, 'label')
            hdr.label = header.label;
        end
    else
        % create a minimal header
        hdr.nChans  = numel(header.label);
        hdr.nSamplesPre = 0; % FIXME perhaps better to make it nan?; And what about nSamplesPst?
        if iscell(header.trial)
            hdr.nTrials  = numel(header.trial);
            hdr.nSamples = numel(header.time{1}); % variable length trials will be caught below 
            if isfield(header, 'fsample')
                hdr.Fs = header.fsample;
            else
                hdr.Fs = 1./mean(diff(header.time{1}));
            end
        else
            hdr.nTrials  = size(header.trial,1);
            hdr.nSamples = numel(header.time);
            if isfield(header, 'fsample')
                hdr.Fs = header.fsample;
            else
                hdr.Fs = 1./mean(diff(header.time));
            end
        end
        
    end
    
    EEG = pop_fileio(hdr, header);
    if iscell(EEG.data) && length(EEG.data) == 1
        EEG.data  = EEG.data{1};
        EEG.times = header.time{1};
    elseif iscell(EEG.data)
        len = cellfun(@(x)size(x,2), EEG.data);
        unLen = unique(len);
        if length(unLen) > 1
            if length(unLen) > 2 || diff(unLen) > 1
                error('Epochs of different length, conversion is not possible');
            else
                fprintf(2, 'Epochs of different length, but small difference in sample\n');
                fprintf(2, 'Assuming outlier due to fractional sample limits and selecting the most common epoch length\n');
                e1 = sum(len == unLen(1));
                e2 = sum(len == unLen(2));
                if min(e1, e2) > 5
                    error('Too many epochs to remove, check data');
                end
                if e1 > e2
                    fprintf(2, 'Truncating %d epochs\n', e2);
                    epochsToTruncate = find(len == unLen(2));
                    truncateLen = unLen(1);
                else
                    fprintf(2, 'Truncating %d epochs\n', e1);
                    epochsToTruncate = find(length == unLen(1));
                    truncateLen = unLen(2);
                end
                for iEpoch = 1:length(epochsToTruncate)
                    EEG.data{epochsToTruncate(iEpoch)} = EEG.data{epochsToTruncate(iEpoch)}(:,1:end-1);
                end
            end
        end
        % EEG.data = [ EEG.data{:} ];
        EEG.data = cat(3, EEG.data{:});
        EEG.pnts = len(1);
    else
        % error('Unknown fieldtrip data format');
        % the trial field was a 3D matrix
        EEG.data  = permute(EEG.data, [2 3 1]);
        
        % overrule the previously created metadata
        [EEG.nbchan, EEG.pnts, EEG.trials] = size(EEG.data);
    end
    if iscell(header.time) % assuming uniformity of trials
        EEG.xmin = header.time{1}(1);
        EEG.xmax = header.time{1}(end);
    else
        EEG.xmin = header.time(1);
        EEG.xmax = header.time(end);
    end
    EEG.trials = size(EEG.data,3);
    EEG = eeg_checkset(EEG);
end   
