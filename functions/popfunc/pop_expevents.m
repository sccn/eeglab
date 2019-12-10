% pop_expevents() - export events to tab separated text file
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

% Revision: 0.10
% Date: 09/16/2011

% Revision history:
%   0.10: Initial version

function com = pop_expevents(EEG, filename, unit)

com = '';
if nargin < 1 
    help pop_expevents;
    return
end

if nargin < 2
	[filename, filepath] = uiputfile('*.txt', 'File name -- pop_expevents()'); 
    drawnow;
	if filename == 0
        return;  end
	filename = [filepath filename];
end

if nargin < 3
    unit = 'samples';  end

if ~(strcmp(unit, 'samples') || strcmp(unit, 'seconds'))
    error('Unit must be either ''samples'' or ''seconds''.');  end

eeg_eventtable(EEG, 'dispTable', false, 'exportFile', filename, 'unit', unit);

com = sprintf('pop_expevents(EEG, ''%s'', ''%s'');', filename, unit); 
