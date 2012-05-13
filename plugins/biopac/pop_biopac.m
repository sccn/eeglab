% pop_biopac() - Read Biosemi 24-bit BDF file 
%
% Usage:
%   >> EEG = pop_biopac;             % an interactive window pops up
%   >> EEG = pop_biopac( filename ); % no pop-up window 
%
% Inputs:
%   filename       - BIOPAC MATLAB file
% Outputs:
%   EEG            - EEGLAB data structure
%
% Author: Arnaud Delorme, SCCN, UCSD, 2011-
%
% See also: openbdf(), biopac()

% Copyright (C) Arnaud Delorme, SCCN, UCSD, 2011-
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

function [EEG, command] = pop_biopac(filename); 
EEG = [];
command = '';

if nargin < 1 
	% ask user
	[filename, filepath] = uigetfile('*.MAT;*.mat', 'Choose a BIOPAC .MAT file -- pop_biopac()'); 
    drawnow;
    
	if filename == 0 return; end;
	filename = [filepath filename];
end;

% load data
% ---------
try
    tmp = load('-mat', filename);
catch, 
    error('Could not read file');
end;

% convert data
% ------------
EEG = eeg_emptyset;
EEG.data = tmp.data';
for index = 1:size(tmp.labels,1)
    EEG.chanlocs(index).labels = tmp.labels(index,:);
end;
if strcmpi(tmp.isi_units, 'ms')
    EEG.srate = 1000/tmp.isi;
elseif strcmpi(tmp.isi_units, 's')
    EEG.srate = 1/tmp.isi;
end;
command = sprintf('EEG = pop_biopac(''%s'');', filename);
