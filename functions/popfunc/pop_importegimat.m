% pop_importegimat() - import EGI Matlab segmented file
%
% Usage:    
%   >> EEG = pop_importegimat(filename, srate, latpoint0); 
%
% Inputs:
%   filename    - Matlab file name
%   srate       - sampling rate
%   latpoint0   - latency in sample ms of stimulus presentation.
%                 When data files are exported using Netstation, the user specify
%                 a time range (-100 ms to 500 ms for instance). In this 
%                 case, the latency of the stimulus is 100 (ms). Default is 0 (ms)
%
% Output:
%   EEG        - EEGLAB dataset structure
%
% Authors: Arnaud Delorme, CERCO/UCSD, Jan 2010

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function [EEG com] = pop_importegimat(filename, srate, latpoint0, dataField);

    EEG = [];
    com = '';
    if nargin < 3, latpoint0 = 0; end
    if nargin < 4, dataField = 'Session'; end
    if nargin < 1 
        % ask user
        [filename, filepath] = uigetfile('*.mat', 'Choose a Matlab file from Netstation -- pop_importegimat()'); 
        if filename == 0 return; end
        filename = fullfile(filepath, filename);
        tmpdata = load('-mat', filename);
        fieldValues = fieldnames(tmpdata);
        sessionPos = strmatch('Session', fieldValues);
        posFieldData = 1; 
        if ~isempty(sessionPos), posFieldData = sessionPos; end
        
        if ~isfield(tmpdata, 'samplingRate'), srate = 250; else srate = tmpdata.samplingRate; end
        % epoch data files only
        promptstr    = { { 'style' 'text'       'string' 'Sampling rate (Hz)' } ...
                         { 'style' 'edit'       'string'  int2str(srate) } ...
                         { 'style' 'text'       'string' 'Sample latency for stimulus (ms)' } ...
                         { 'style' 'edit'       'string' '0' } ...
                         { 'style' 'text'       'string' 'Field containing data' } ...
                         { 'style' 'popupmenu'  'string' fieldValues 'value' posFieldData } ...
                         };
        geometry = { [1 1] [1 1] [1 1] };
        result       = inputgui( 'geometry', geometry, 'uilist', promptstr, ...
                                 'helpcom', 'pophelp(''pop_importegimat'')', ...
                                 'title', 'Import a Matlab file from Netstation -- pop_importegimat()');

        if length(result) == 0 return; end;        
        srate = str2num(result{1});
        latpoint0 = str2num(result{2});
        dataField = fieldValues{result{3}};
        if isempty(latpoint0), latpoint0 = 0; end
    end
   
    EEG = eeg_emptyset;
    fprintf('Reading EGI Matlab file %s\n', filename);
    tmpdata = load('-mat', filename);
    if isfield(tmpdata, 'samplingRate') % continuous file
        srate = tmpdata.samplingRate;
    end

    fieldValues = fieldnames(tmpdata);
    if all(cellfun(@(x)isempty(findstr(x, 'Segment')), fieldValues))
        EEG.srate = srate;
        indData = strmatch(dataField, fieldValues);
        EEG.data = tmpdata.(fieldValues{indData(1)});
        EEG = eeg_checkset(EEG);
        EEG = readegilocs(EEG);
        
        com = sprintf('EEG = pop_importegimat(''%s'');', filename);
    else
        % get data types
        % --------------
        allfields = fieldnames(tmpdata);
        for index = 1:length(allfields)
            allfields{index} = allfields{index}(1:findstr(allfields{index}, 'Segment')-2);
        end
        datatypes = unique_bc(allfields);
        datatypes(cellfun(@isempty, datatypes)) = [];

        % read all data
        % -------------
        counttrial = 1;
        EEG.srate  = srate;
        latency = (latpoint0/1000)*EEG.srate+1;
        for index = 1:length(datatypes)
            tindex = 1;
            for tindex = 1:length(allfields)
                if isfield(tmpdata, sprintf('%s_Segment%d', datatypes{index}, tindex))
                    datatrial = getfield(tmpdata, sprintf('%s_Segment%d', datatypes{index}, tindex));
                    if counttrial == 1
                        EEG.pnts = size(datatrial,2);
                        EEG.data = repmat(single(0), [size(datatrial,1), size(datatrial,2), 1000]);
                    end
                    EEG.data(:,:,counttrial) = datatrial;

                    EEG.event(counttrial).type    = datatypes{index};
                    EEG.event(counttrial).latency = latency;
                    EEG.event(counttrial).epoch   = counttrial;

                    counttrial = counttrial+1;
                    latency = latency + EEG.pnts;
                end
            end
        end
        fprintf('%d trials read\n', counttrial-1);
        EEG.data(:,:,counttrial:end) = [];

        EEG.setname  = filename(1:end-4);
        EEG.nbchan   = size(EEG.data,1);
        EEG.trials   = counttrial-1;
        if latpoint0 ~= 1
            EEG.xmin     = -latpoint0/1000;
        end
        EEG = eeg_checkset(EEG);

        % channel location
        % ----------------
        if all(EEG.data(end,1:10) == 0)
            disp('Deleting empty data reference channel (reference channel location is retained)');
            EEG.data(end,:)   = [];
            EEG.nbchan        = size(EEG.data,1);
            EEG = eeg_checkset(EEG);
        end
        EEG = readegilocs(EEG);

        com = sprintf('EEG = pop_importegimat(''%s'', %3.2f, %3.2f, %d);', filename, srate, latpoint0, dataField);
    end
    
