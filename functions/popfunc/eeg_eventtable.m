function events = eeg_eventtable(EEG, varargin)
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

% Copyright by Clemens Brunner <clbrunner@ucsd.edu>
% Revision: 0.15
% Date: 10/13/2011

% Revision history:
%   0.15: Added option to add/remove index column
%   0.11: Changed function name and updated documentation
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
        end;
    end;
end;

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

