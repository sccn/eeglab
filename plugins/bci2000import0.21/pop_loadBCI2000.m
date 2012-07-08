function EEG = pop_loadBCI2000(fileName, events)
% pop_loadBCI2000() - loads BCI2000 .dat files into EEGLAB
%
% Usage:
%   >> EEG = pop_loadBCI2000();
%   >> EEG = pop_loadBCI2000(fileName);
%   >> EEG = pop_loadBCI2000(fileName, events);
%
% Inputs:
%   fileName - File name(s), either a string (possibly containing wildcards) or
%              a cell array of strings containing multiple file names
%   events   - List of events to import, can either be numeric (a vector) or a
%              cell array of strings containing the names of the events to
%              import (default: all events are imported)
%
% Outputs:
%   EEG      - EEGLAB data structure
%
% Examples:
%   Open GUI dialog:
%     >> EEG = pop_loadBCI2000();
%
%   Load all .dat files in the current folder with all events:
%     >> EEG = pop_loadBCI2000('*.dat');
%
%   Load all .dat files with event numbers 2 and 5:
%     >> EEG = pop_loadBCI2000('*.dat', [2, 5]);
%
%   Load all .dat files with events 'StimulusCode' and 'Feedback':
%     >> EEG = pop_loadBCI2000('*.dat', {'StimulusCode', 'Feedback'});
%
%   Load two specific files:
%     >> EEG = pop_loadBCI2000({'set001.dat', 'set002.dat'});

% Copyright by Clemens Brunner <clbrunner@ucsd.edu>
% Revision: 0.30
% Date: 09/14/2011
% Parts based on BCI2000import.m from the BCI2000 distribution (www.bci2000.org)

% Revision history:
%   0.30: Create boundary events when loading multiple files
%   0.26: Import channel labels
%   0.25: Added EEG.urevent structure
%   0.20: Adapted to EEGLAB conventions
%   0.10: Initial version

% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 59 Temple
% Place - Suite 330, Boston, MA  02111-1307, USA.

if nargin < 1  % No input arguments specified, show GUI
    [fileName, filePath] = uigetfile('*.dat', 'Choose BCI2000 file(s) -- pop_loadBCI2000', 'multiselect', 'on');
    
    if ~iscell(fileName)  % If only one file was selected
        [signal, states, parms] = load_bcidat(fullfile(filePath, fileName));
    else  % Multiple files were selected
        files = struct('name', fileName);
        totalSamples = zeros(1, length(files));
        for k = 1:length(files)
            files(k).name = fullfile(filePath, files(k).name);  % Add full path to each file name
            [~, ~, ~, totalSamples(k)] = load_bcidat(files(k).name, [0, 0]);  % Get length of each individual file
        end
        
        % Determine start time of individual files
        fileBorders = zeros(1, length(files));
        fileBorders(1) = 1;
        for k = 2:length(files)
            fileBorders(k) = fileBorders(k - 1) + totalSamples(k - 1);  end
        
        [signal, states, parms] = load_bcidat(files.name);  % Load all files
    end
    
    % Build a string consisting of all events, separated by a "|" (necessary for
    % inputgui() below)
    eventList = fieldnames(states);
    eventString = [];
    for k = 1:length(eventList)
        if k > 1
            eventString = [eventString, '|', eventList{k}];
        else
            eventString = eventList{k};
        end
    end
    
    temp = inputgui('geometry', {1, 1}, 'geomvert', [1, 10], ...
        'uilist', {{'style', 'text', 'string', 'Select events'}, ...
        {'style', 'listbox', 'string', eventString, 'min', 0, 'max', 2}});
    
    if isempty(temp)  % Abort if user clicked "Cancel"
        error('Loading aborted.'); end
    
    events = temp{1};
    clear('temp');
    
else  % Input arguments specified
    if iscell(fileName)  % Multiple file names in a cell array
        files = cell2struct(fileName, 'name', 1);
    else  % Just a string (possibly containing wildcards)
        files = dir(fileName);
        totalSamples = zeros(1, length(files));
        for k = 1:length(files)
            files(k).name = fullfile(fileparts(fileName), files(k).name);  % Add full path to each file name
            [~, ~, ~, totalSamples(k)] = load_bcidat(files(k).name, [0, 0]);  % Get length of each individual file
        end
    end
    
    % Determine start time of individual files
    if length(files) > 1
        fileBorders = zeros(1, length(files));
        fileBorders(1) = 1;
        for k = 2:length(files)
            fileBorders(k) = fileBorders(k - 1) + totalSamples(k - 1);  end
    end
    
    [signal, states, parms] = load_bcidat(files.name);  % Load all files
    
    if ~exist('events', 'var')  % If no events were specified, select all events contained in the file
        events = 1:length(fieldnames(states)); end
end

% Select events
stateNames = fieldnames(states);

% Convert event names (strings) into numbers
if iscell(events)
    temp = [];  % This vector will hold the positions of all events contained in the signal
    for k = 1:length(events)
        [isState, stateNumber] = ismember(events{k}, stateNames);
        if isState
            temp = [temp, stateNumber];
        else
            warning('''%s'' is not an event in the signal. Therefore, this event is being ignored.', events{k});
        end
    end
    events = sort(temp);
    clear('temp');
end
selectedEvents = stateNames(events);

% Remove events that were not selected
for k = 1:length(stateNames)
    if ismember(stateNames{k}, selectedEvents)
        continue; end
    states = rmfield(states, stateNames{k});
end

states = compressData(states);
stateNames = fieldnames(states);

EEG = eeg_emptyset();
EEG.setname = 'Imported BCI2000 data set';
EEG.srate = parms.SamplingRate.NumericValue;
EEG.nbchan = size(signal, 2);
EEG.pnts = size(signal, 1);
EEG.trials = 1;
EEG.data = signal';
if ~isempty(parms.ChannelNames.Value)
    for k = 1:EEG.nbchan  % Assign channel labels
        EEG.chanlocs(k).labels = parms.ChannelNames.Value{k}; end
    EEG = eeg_checkset(EEG);
else
    for k = 1:EEG.nbchan  % Use channel numbers as labels
        EEG.chanlocs(k).labels = num2str(k);
    end
    EEG = eeg_checkset(EEG);
end

evCount = 1;
for k = 1:length(stateNames)
    st = getfield(states, stateNames{k});
    for ev = 1:length(st.latency)
        EEG.event(evCount).latency = st.latency(ev);
        EEG.event(evCount).position = st.value(ev);
        EEG.event(evCount).type = stateNames{k};
        evCount = evCount + 1;
    end
end

if exist('fileBorders', 'var')  % If multiple files were loaded and concatenated, create boundary events at the start of new files
    EEG.event = eeg_insertbound(EEG.event, EEG.pnts * EEG.trials, fileBorders');  end

EEG = eeg_checkset(EEG, 'eventconsistency');
EEG.urevent = EEG.event;  % Create urevent structure
for k = 1:length(EEG.event)  % Create urevent pointers
    EEG.event(k).urevent = k; end
EEG = eeg_checkset(EEG);



function data = compressData(data)
% Compresses event data time series by noting only changes in event values.
%
% TODO: Detailed description

% List of known event types to compress
compressStates = {'TargetCode','ResultCode','IntertrialInterval','Feedback','Dwelling','StimulusCode'};

if isstruct(data)
    stateNames = fieldnames(data);
    for k = 1:length(stateNames)
        if ismember(stateNames{k}, compressStates)  % If the event is in the list of known events, compress it
            value = getfield(data, stateNames{k});
            if isstruct(value)
                continue; end
            data = setfield(data, stateNames{k}, compressData(value));
            clear('value');
        else  % If the event is unknown, try to figure out if to compress it
            value = getfield(data, stateNames{k});
            if isstruct(value)  % I don't know when this field can be a struct; obviously, it is ignored then
                continue; end
            % Specifies the percentage of number of different event values relative to the length of the EEG signals (e.g. 0.05 is 5%)
            relativeEventLength = 0.05;
            U = unique(value);
            
            % Compress only if there are at least 2 different event values AND
            % if the number of different event values does not exceed the above set threshold
            if length(U) > 1 && length(U) < length(value)*relativeEventLength
                data = setfield(data, stateNames{k}, compressData(value));
            else
                warning('''%s'' is not a valid event for EEGLAB. This event is being discarded.', stateNames{k});
                data = rmfield(data, stateNames{k});
            end
        end
    end
    
else
    
    lat = find(diff(data) > 0) + 1;
    pos = data(lat);
    
    % TODO: Check calculation of duration, right now it just checks if the value
    % goes back to 0. However, we should probably take the time where it is
    % constant as its duration.
    dur = [];
    for k = 1:length(lat)
        b = find(data(lat(k):end) == 0);
        if ~isempty(b)
            dur(k) = b(1) - 1;
        else
            dur(k) = length(data) - lat(k);
        end
    end
    clear data;
    data.latency = lat;
    data.value = pos;
    data.duration = dur(:);
end
