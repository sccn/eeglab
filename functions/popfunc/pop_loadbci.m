% pop_loadbci() - import BCI2000 ascii/Matlab file into eeglab
%
% Usage:
%   >> OUTEEG = pop_loadbci( filename, srate );
%
% Inputs:
%   filename       - file name
%   srate          - sampling rate
%
% Outputs:
%   OUTEEG         - EEGLAB data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 9 July 2002
%
% See also: eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $

function [EEG, command] = pop_loadbci(filename, srate); 
command = '';

if nargin < 1
	% ask user
	[filename, filepath] = uigetfile('*.DAT', 'Choose a BCI file -- pop_loadbci'); 
	if filename == 0 return; end;
	filename = [filepath '/' filename];
	promptstr    = { 'Sampling rate' };
	inistr       = { '256' };
	result       = inputdlg( promptstr, 'Import BCI2000 data -- pop_loadbci()', 1,  inistr);
	if length(result) == 0 return; end;
	srate   = eval( result{1} );
end;

% try to read as matlab
EEG = eeg_emptyset;
fprintf('Pop_loadbci: importing BCI file...\n');
try
	load( filename, '-mat');
catch
	
	fid = fopen(filename, 'r');
	allcollumns = fgetl(fid);
	colnames = {};
	while ~isempty(deblank(allcollumns))
		[colnames{end+1} allcollumns] = strtok(allcollumns);
	end;
	tmpdata = fscanf(fid, '%f', Inf);
	tmpdata = reshape(tmpdata, length(colnames), length(tmpdata)/length(colnames));
	%tmpdata = tmpdata';
	
	% find electrode indices
	% ----------------------
	indices = [];
	for index = 1:length(colnames)
		if strcmp( colnames{index}(1:2), 'ch')
			indices = [indices index ];
		end;
	end;
	eventindices = setdiff(1:length(colnames), indices);
	eventindices(3 + [ 1 2 3 4 7 9 10 11 12]) = [];
	eventindices(1:3) = []; % suppress these event 
	
	% process events
	% --------------
	fprintf('Pop_loadbci: importing events...\n');
	counte = 1; % event counter
	events(10000).latency = 0;
	for index = eventindices
		tmpevent = find( diff(tmpdata(index, :)) ~= 0);
		tmpevent = tmpevent+1;
		for tmpi = tmpevent
			events(counte).type = colnames{index};
			events(counte).latency = tmpi;
			events(counte).value   = tmpdata(index, tmpi);
			counte = counte+1;
			%if mod(counte, 100) == 0, fprintf('%d ', counte); end;
		end;
	end;
	EEG.event = events(1:counte-1);	
	EEG.data = tmpdata(indices,:);
	EEG.nbchan = size(EEG.data, 1);
	EEG.srate  = srate;
	EEG = eeg_checkset(EEG);
	EEG = pop_editeventvals( EEG, 'sort', { 'latency', [0] } );
end;

command = sprintf('EEG = pop_loadbci(''%s'', %f);',filename, srate); 
return;
