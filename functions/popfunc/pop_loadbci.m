% pop_loadbci() - import a BCI2000 ascii file into EEGLAB
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
% Revision 1.9  2003/04/10 18:05:22  arno
% default argument
%
% Revision 1.8  2002/11/13 17:09:12  scott
% help msg
%
% Revision 1.7  2002/10/15 17:00:49  arno
% drawnow
%
% Revision 1.6  2002/08/12 16:33:37  arno
% inputdlg2
%
% Revision 1.5  2002/07/11 17:18:15  arno
% updating header
%
% Revision 1.4  2002/07/10 23:44:17  arno
% debugging events .
%
% Revision 1.3  2002/07/09 23:30:50  arno
% selecting interesting events
%
% Revision 1.2  2002/07/08 19:28:15  arno
% changing function filter
%
% Revision 1.1  2002/07/08 19:20:53  arno
% Initial revision
%

function [EEG, command] = pop_loadbci(filename, srate); 
EEG = [];
command = '';

if nargin < 1
	% ask user
	[filename, filepath] = uigetfile('*.*', 'Choose a BCI file -- pop_loadbci'); 
    drawnow;
	if filename == 0 return; end;
	filename = [filepath '/' filename];
	promptstr    = { 'Sampling rate' };
	inistr       = { '256' };
	result       = inputdlg2( promptstr, 'Import BCI2000 data -- pop_loadbci()', 1,  inistr, 'pop_loadbci');
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
	ISIind = eventindices(3 + 9);
	eventindices(3 + [ 1 2 3 4 7 8 9 10 11 12]) = [];
	eventindices(1:3) = []; % suppress these event 
	
	% add the trial number
	% --------------------
	tmptrial = find( diff(tmpdata(ISIind, :)) ~= 0);
	tmptrial = tmptrial+1;

	% process events
	% --------------
	fprintf('Pop_loadbci: importing events...\n');
	counte = 1; % event counter
	events(10000).latency = 0;
	for index = eventindices
		counttrial = 1;
		tmpevent = find( diff(tmpdata(index, :)) ~= 0);
		tmpevent = tmpevent+1;
		for tmpi = tmpevent
			if tmpdata(index, tmpi)
				events(counte).type    = [ colnames{index} int2str(tmpdata(index, tmpi)) ];
				events(counte).latency = tmpi;
				%events(counte).value   = tmpdata(index, tmpi);
				%while tmpi > tmptrial(counttrial) & counttrial < length(tmptrial)
				%	counttrial = counttrial+1;
				%end;
				%events(counte).trial = counttrial;				
				counte = counte+1;
				%if mod(counte, 100) == 0, fprintf('%d ', counte); end;
			end;
		end;
	end;
	
	% add up or down events
	% ---------------------
	EEG.data = tmpdata(indices,:);
	EEG.nbchan = size(EEG.data, 1);
	EEG.srate  = srate;
	EEG = eeg_checkset(EEG);
	EEG.event = events(1:counte-1);	
	EEG = pop_editeventvals( EEG, 'sort', { 'latency', [0] } );
	for index=1:length(EEG.event)
		if strcmp(EEG.event(index).type(1:6), 'Target')
			targetcode = str2num(EEG.event(index).type(end));
			if targetcode == 1
				EEG.event(index).type = 'toptarget';
			else
				EEG.event(index).type = 'bottomtarget';
			end;
		else 
			if strcmp(EEG.event(index).type(1:6), 'Result')
				resultcode = str2num(EEG.event(index).type(end));
				if resultcode == 1
					EEG.event(index).type = 'topresp';
				else
					EEG.event(index).type = 'bottomresp';
				end;
				EEG.event(end+1).latency = EEG.event(index).latency;
				if (resultcode == targetcode) 
					EEG.event(end).type = 'correct';
				else
					EEG.event(end).type = 'miss';
				end;
			end;
		end;
	end;
	EEG = pop_editeventvals( EEG, 'sort', { 'latency', [0] } );
	%EEG.data = tmpdata([72 73 75],:);
end;

command = sprintf('EEG = pop_loadbci(''%s'', %f);',filename, srate); 
return;
