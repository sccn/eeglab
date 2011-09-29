function com = pop_expevents(EEG, filename, unit)
% pop_expevents() - export events to CSV file
%
% Usage:
%   >> pop_expevents(EEG);
%   >> pop_expevents(EEG, filename);
%   >> pop_expevents(EEG, filename, unit);
%
% Inputs:
%   EEG      - EEGLAB dataset
%   filename - text file name
%   unit     - display latency and duration in 'samples' (default) or 'seconds'
%
% Outputs:
%   com      - Command to execute this function from the command line
%
% Examples:
%   Export all events and show a dialog window:
%     >> pop_expevents(EEG);
%
%   Export all events to a file, latency and duration in samples:
%     >> pop_expevents(EEG, filename);
%
%   Export all events to a file, latency and duration in seconds:
%     >> pop_expevents(EEG, filename, 'seconds');

% Copyright by Clemens Brunner <clbrunner@ucsd.edu>
% Revision: 0.10
% Date: 09/16/2011

% Revision history:
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

com = '';
if nargin < 1 
    help pop_expevents;
    return
end

if nargin < 2
	[filename, filepath] = uiputfile('*.csv', 'File name -- pop_expevents()'); 
    drawnow;
	if filename == 0
        return;  end
	filename = [filepath filename];
end;

if nargin < 3
    unit = 'samples';  end

if ~(strcmp(unit, 'samples') || strcmp(unit, 'seconds'))
    error('Unit must be either ''samples'' or ''seconds''.');  end

eeg_eventtable(EEG, 'dispTable', false, 'exportFile', filename, 'unit', unit);

com = sprintf('pop_expevents(%s, ''%s'', ''%s'');', inputname(1), filename, unit); 
