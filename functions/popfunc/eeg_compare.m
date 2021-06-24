% eeg_compare - compare EEG structures. Differences are shown on the command
%               line.
% Usage:
%   eeg_compare(EEG1, EEG2);
%
% Input:
%  EEG1 - first EEGLAB structure
%  EEG2 - second EEGLAB structure
%
% Author: Arnaud Delorme, 2020
%
% See also: eeglab(), eegplot(), pop_rejepoch()

% Copyright (C) 2020 Arnaud Delorme
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

function eeg_compare(EEG, EEG2)

%% Assess difference between datasets
fields = fieldnames(EEG);
for iField = 1:length(fields)
    if ~isfield(EEG2, fields{iField})
        fprintf('Field %s missing in second dataset\n', fields{iField});
    else
        if ~isequal(EEG.(fields{iField}), EEG2.(fields{iField}))
            fprintf('Field %s differs\n', fields{iField});
        end
    end
end
if ~isequal(EEG.xmin, EEG2.xmin), fprintf('Difference between xmin is %1.6f sec\n', EEG.xmin-EEG2.xmin); end
if ~isequal(EEG.xmax, EEG2.xmax), fprintf('Difference between xmax is %1.6f sec\n', EEG.xmax-EEG2.xmax); end

% check chanlocs
[~,~,chanlocs1] = eeg_checkchanlocs( EEG.chanlocs, EEG.chaninfo);
[~,~,chanlocs2] = eeg_checkchanlocs( EEG2.chanlocs, EEG.chaninfo);
if length(chanlocs1) == length(chanlocs2)
    differ = 0;
    differLabel = 0;
    for iChan = 1:length(chanlocs1)
        if  sum(abs([ chanlocs1(iChan).X chanlocs1(iChan).Y chanlocs1(iChan).Z] - ...
                [ chanlocs2(iChan).X chanlocs2(iChan).Y chanlocs2(iChan).Z])) > 1e-12
            differ = differ+1;
        end
        if ~isequal(chanlocs1(iChan).labels, chanlocs2(iChan).labels)
            differLabel = differLabel+1;
        end
    end
    if differ
        fprintf('%d channel coordinates differ\n', differ);
    else
        disp('All channel coordinates are OK');
    end
    if differLabel
        fprintf('%d channel label(s) differ\n', differLabel);
    else
        disp('All channel labels are OK');
    end
else
    disp('Different numbers of channels');
end

% check events
if length(EEG.event) ~= length(EEG2.event)
    disp('Different numbers of events');
elseif isempty(EEG.event)
    disp('All events OK (empty)');
else
    fields1 = fieldnames(EEG.event);
    fields2 = fieldnames(EEG2.event);
    allFieldsOK = true;
    
    if ~isequal(sort(fields1), sort(fields2))
        disp('Not the same number of event fields');
        allFieldsOK = false;
    end
    
    for iField = 1:length(fields1)
        if isfield(EEG.event, fields1{iField}) && isfield(EEG2.event, fields1{iField})
            diffVal = zeros(1,length(EEG.event));
            if strcmpi(fields1{iField}, 'latency')
                for iEvent = 1:length(EEG.event)
                    diffVal(iEvent) = EEG.event(iEvent).(fields1{iField}) - EEG2.event(iEvent).(fields1{iField});
                end
            else
                for iEvent = 1:length(EEG.event)
                    diffVal(iEvent) = ~isequal(EEG.event(iEvent).(fields1{iField}), EEG2.event(iEvent).(fields1{iField}));
                end
            end
            if any(diffVal ~= 0)
                if strcmpi(fields1{iField}, 'latency')
                    fprintf('Event latency (%2.1f %%) are not OK (abs diff of these is %1.4f samples)\n', length(find(diffVal))/length(diffVal)*100, mean( abs(diffVal(diffVal ~=0 ))) );
                    fprintf(' ******** (see plot)\n');
                    figure; plot(diffVal);
                else
                    fprintf('Event fields "%s" are NOT OK (%2.1f %% of them)\n', fields1{iField}, length(find(diffVal))/length(diffVal)*100);
                end
                allFieldsOK = false;
            end
        end
    end
    disp('All other events OK (unless warning above)');
end

% check epochs
if length(EEG.epoch) == length(EEG2.epoch)
    if ~isempty(EEG.epoch)
        fields1 = fieldnames(EEG.epoch);
        fields2 = fieldnames(EEG2.epoch);
        allFieldsOK = true;
        if ~isequal(sort(fields1), sort(fields2))
            disp('Not the same number of event fields');
            allFieldsOK = false;
        else
            diffVal = [];
            for iField = 1:length(fields1)
                for iEpoch = 1:length(EEG.epoch)
                    diffVal(iEpoch) = ~isequal(EEG.epoch(iEpoch).(fields1{iField}), EEG2.epoch(iEpoch).(fields1{iField}));
                end
                if any(diffVal ~= 0)
                    fprintf('Epoch fields "%s" are NOT OK (%2.1f %% of them)\n', fields1{iField}, length(find(diffVal))/length(diffVal)*100);
                    allFieldsOK = false;
                end
            end
        end
        if allFieldsOK
            disp('All epoch and all epoch fields are OK');
        end        
    end
elseif ~isempty(EEG.epoch)
    disp('Different numbers of epochs');
end
