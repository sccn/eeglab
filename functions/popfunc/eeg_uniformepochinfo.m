% eeg_uniformepochinfo() - make uniform event field information within 
%                          epochs by filling in empty field values with 
%                          values of other events in the same epoch. This
%                          is particularly useful when time-locking event
%                          do not contain information of other events in
%                          the epoch.
% Usage: 
%   >> EEG = eeg_uniformepochinfo(EEG);  
%
% Inputs:
%   EEG      - EEGLAB dataset
%
% Optional input:
%   'checkepochs' - ['off'|'warning'|'error'] issue warning or errors when
%                   some epochs have conflicting information or are missing
%                   information. Default is 'off'.
%
% Inputs:
%   EEG      - EEGLAB dataset with updated event structure
%
% Authors: Arnaud Delorme, SCCN/INC/UCSD, November 2019

% Copyright (C) Arnaud Delorme arno@ucsd.edu
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

% copy field content within epochs so the time-locking event has more information
% -------------------------------------------------------------------------------
function EEG = eeg_uniformepochinfo(EEG, varargin)

tmpevent  = EEG.event;
if ~isfield(tmpevent, 'epoch')
    error('Data must contain epochs to use this function');
end
allepochs = [ tmpevent.epoch ];

opt = finputcheck(varargin, { 'checkepochs' 'string' { 'warning' 'error' 'off' } 'off' }, 'eeg_uniformepochinfo');
if ischar(opt), error(opt); end

% uniformize fields content for the different epochs
% --------------------------------------------------
difffield = fieldnames(EEG.event);
difffield = difffield(~(strcmp(difffield,'latency')|strcmp(difffield,'epoch')|strcmp(difffield,'type')|strcmp(difffield,'mffkeys')|strcmp(difffield,'mffkeysbackup')|strcmp(difffield,'begintime')));
for index = 1:length(difffield)
    tmpevent  = EEG.event;
    allvalues = { tmpevent.(difffield{index}) };
    try
        valempt = cellfun('isempty', allvalues);
    catch
        valempt = mycellfun('isempty', allvalues);
    end
    arraytmpinfo = cell(1,EEG.trials);
    if ischar(allvalues{1}), arraytmpinfo(:) = {''}; end % preserves uniformity of char or num
    
    % get the field content
    % ---------------------
    indexevent = find(~valempt); % non empty
    if length(unique(allepochs(indexevent))) > length(allepochs(indexevent))
        msg = 'Duplicate event information detected when uniformizing epoch information';
        if strcmpi(opt.checkepochs, 'warning')
            warning(msg);
        elseif strcmpi(opt.checkepochs, 'error')
            error(msg);
        end
    end
    arraytmpinfo(allepochs(indexevent)) = allvalues(indexevent);

    % uniformize content for all epochs
    % ---------------------------------
    indexevent = find(valempt); % empty
    tmpevent   = EEG.event;
    [tmpevent(indexevent).(difffield{index})] = arraytmpinfo{allepochs(indexevent)};
    if ~isequal(EEG.event, tmpevent) % might be slow?
        EEG.event  = tmpevent;
        fprintf(['eeg_uniformepochinfo: found empty values for field ''' difffield{index} ''' for %d epochs\n'], sum(valempt));
        fprintf(['                      filling with values of other events in the same epochs\n']);
        if mod(sum(valempt), EEG.trials) ~= 0
            msg = sprintf('                      note: epochs not uniform (should be multiple of %d epochs)', EEG.trials);
            disp(msg);
        end        
    end    
end