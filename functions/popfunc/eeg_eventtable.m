% eeg_eventtable() - returns all events contained in the EEG structure (and
%                    optionally exports them to a CSV file)
%
% Usage:
%   >> events = eeg_eventtable(EEG);
%
% Inputs:
%   EEG          - EEG structure
%
% Optional inputs:
%   'unit'       - Unit of events containing time samples, can be either
%                  'samples' (default) or 'seconds'
%   'dispTable'  - Display an overview of all events if set to true (default);
%                  do not display if set to false
%   'exportFile' - Exports events as a CSV file (using tabs as delimiters); set
%                  this parameter to the file name you wish to export the events
%                  to (string); default: no file is exported
%   'indexCol'   - Adds a column with the event indices as the first row; note 
%                  that this is not an event field; default: true
%
% Outputs:
%   events       - event structure (cell array)
%
% Examples:
%   Get all events contained in the EEG structure and display an overview:
%     >> events = eeg_eventtable(EEG);
%
%   In addition to displaying an overview, export into a CSV file:
%     >> events = eeg_eventtable(EEG, 'exportFile', 'test.csv');

% Copyright by Clemens Brunner <clbrunner@ucsd.edu>, 2011
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

% Revision: 0.15
% Date: 10/13/2011

% Revision history:
%   0.15: Added option to add/remove index column
%   0.11: Changed function name and updated documentation
%   0.10: Initial version

function events = eeg_eventtable(EEG, varargin)

% Default parameters, can be overwritten by varargin
unit = 'samples';  % Unit of latencies is 'samples'
dispTable = true;  % Display the event table
exportFile = [];  % Export events to text file
indexCol = true;  % Add an index column as first column

if ~isempty(varargin)  % Are there optional parameters available?
    k = 1;
    while k <= length(varargin)
        if strcmp(varargin{k}, 'dispTable')
            dispTable = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'unit')
            unit = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'exportFile')
            exportFile = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'indexCol')
            indexCol = varargin{k + 1};
            k = k + 2;
        else  % Ignore unknown parameters
            k = k + 2;
        end
    end
end

fields = fieldnames(EEG.event);  % Get all event fields
nFields = length(fields);  % Number of event fields
nEvents = length(EEG.event);  % Number of events

events = cell(nFields, 1);
for k = 1:nFields
    events{k}.name = fields{k};  % Event name
    events{k}.values = cell(nEvents, 1);  % Will contain values for this event
    
    for l = 1:nEvents  % Write values
        if strcmp(events{k}.name, 'latency') && strcmp(unit, 'seconds')
            events{k}.values{l} = getfield(EEG.event(l), fields{k})/EEG.srate;
        elseif strcmp(events{k}.name, 'duration') && strcmp(unit, 'seconds')
            events{k}.values{l} = getfield(EEG.event(l), fields{k})/EEG.srate;
        else
            events{k}.values{l} = getfield(EEG.event(l), fields{k});
        end
    end

end

if dispTable || ~isempty(exportFile)  % Display overview or export to CSV file?
    
    nr = cell(nEvents, 1);  % Consecutive event numbers
    for l = 1:nEvents
        nr{l} = l;  end
    
    eventsTable = [];
    if indexCol  % Add index column
        titleTable{1} = 'number';
        titleOffset = 1;
    else  % Do not add index column
        titleOffset = 0;
    end
    for k = 1:nFields
        eventsTable = [eventsTable, events{k}.values];
        titleTable{k + titleOffset} = events{k}.name;
    end
    if indexCol
        eventsTable = [nr, eventsTable];  % Add event numbers in first column
    end
    finalTable = [titleTable; eventsTable];  % Add title in first row
    
    if dispTable
        format('shortG');
        disp(finalTable);
        format;
    end
    
    if ~isempty(exportFile)
        fid = fopen(exportFile, 'w');
        for rows = 1:size(finalTable, 1)
            for cols = 1:size(finalTable, 2)
                dataValue = finalTable(rows, cols);
                if cols < size(finalTable, 2)  % If we're not in the last column, print the delimiter
                    if iscellstr(dataValue)
                        fprintf(fid, '%s\t', dataValue{:});
                    elseif iscellnumeric(dataValue)
                        fprintf(fid, '%f\t', dataValue{:});
                    end
                else  % If we're in the last column, do not print the delimiter
                    if iscellstr(dataValue)
                        fprintf(fid, '%s', dataValue{:});
                    elseif iscellnumeric(dataValue)
                        fprintf(fid, '%f', dataValue{:});
                    end
                end
            end
            fprintf(fid, '\n');
        end
        fclose(fid);
    end
end

function b = iscellnumeric(C)
% Return 1 if all elements of cell array are numeric

b = all(cellfun(@(x) isnumeric(x),C));


