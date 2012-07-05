% pop_read_erpss() - interactively import an uncompressed ERPSS-format 
%                    data file (.RAW or .RDF) using read_erpss()
% Usage:
%   >> OUTEEG = pop_read_erpss( filename, srate );
%
% Inputs:
%   filename       - file name (with extension)
%
% Outputs:
%   OUTEEG         - EEGLAB data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 23 January 2003
%
% See also: eeglab(), read_erpss()

% Copyright (C) 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [EEG, command] = pop_read_erpss(filename, srate); 
EEG = [];
command = '';

if nargin < 1
	% ask user
	[filenametmp, filepath] = uigetfile('*.rdf;*.RDF;*.raw;*.RAW', 'Choose a ERPSS file -- pop_read_erpss'); 
    drawnow;
	if filenametmp == 0 return; end;
	filename = [filepath '/' filenametmp];
end;

% read ERPSS format
EEG = eeg_emptyset;
fprintf('pop_read_erpss: importing ERPSS file...\n');
[EEG.data,events,header] = read_erpss(filename);
EEG.nbchan = size(EEG.data,1);
if nargin < 1 & round(header.srate) == 0
    promptstr    = { 'Sampling rate' };
    inistr       = { '256' };
    result       = inputdlg2( promptstr, 'Import BCI2000 data -- pop_loadbci()', 1,  inistr, 'pop_loadbci');
    if length(result) == 0 return; end;
    srate   = eval( result{1} );
elseif round(header.srate) ~= 0 
    srate = header.srate;
elseif nargin ~= 2
    disp('WARNING: Unknown sampling rate.Use menu "Edit > Dataset info" to enter it.');
    srate = NaN;
end;
EEG.srate = srate;
EEG.setname = 'ERPSS data';
EEG.comments        = [ 'Original file: ' filename ];
if ~isempty(events)
    EEG.event = struct( 'type', { events.event_code }, 'latency', {events.sample_offset});
end;
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'makeur');

command = sprintf('EEG = pop_read_erpss(''%s'', %f);',filename, srate); 
return;
