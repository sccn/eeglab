% pop_loaddat() - merge a neuroscan DAT file with input dataset
%                (pop out window if no arguments).
%
% Usage:
%   >> OUTEEG = pop_loaddat( INEEG ); % pop-up window mode
%   >> OUTEEG = pop_loaddat( INEEG, filename, no_rt);
%
% Graphic interfance:
%   "Code signifying no event ..." - [edit box] reaction time 
%                    no event code. See 'no_rt' command line equivalent
%                    help.
% Inputs:
%   filename       - file name
%   INEEG          - input EEGLAB data structure
%   no_rt          - no reaction time integer code (ex: 1000). Since 
%                    a number has to be defined for each reaction 
%                    time, epochs with no reaction time usually have
%                    a stereotyped reaction time value (such as 1000).
%                    Default none.
% Outputs:
%   OUTEEG         - EEGLAB data structure
%
% Author: Arnaud Delorme, CNL/Salk Institute, 2001
%
% See also: loaddat(), eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% 01-25-02 reformated help & license -ad 
% 13/02/02 removed the no latency option -ad

function [EEG, command] = pop_loaddat(EEG, filename, no_rt); 
command = '';

if nargin < 1
	help  pop_loaddat;
	return;
end;	

if nargin < 2 
	% ask user
	[filename, filepath] = uigetfile('*.DAT', 'Choose a DAT file -- pop_loaddat'); 
    drawnow;
	if filename == 0 return; end;
	result       = inputdlg2( { strvcat('Code signifying no event in a trial ([]=none)', ...
									 '(none=all latencies are imported)')}, ...
									 'Load Neuroscan DATA file -- pop_loaddat()', 1,  {'1000'}, 'pop_loaddat');
	if length(result) == 0 return; end;
	no_rt = eval( result{1} );
end;
if exist('no_rt') ~= 1 | isempty(no_rt)
	no_rt = NaN;
end;

% load datas
% ----------
if exist('filepath')
	fullFileName = sprintf('%s%s', filepath, filename);
else
	fullFileName = filename;
end;	
disp('Loading dat file...');
[typeeeg, rt, response, n] = loaddat( fullFileName );

if n ~= EEG.trials
	error('pop_loaddat, number of trials in input dataset and DAT file different, aborting');
end;	  

for index = 1:length(EEG.event)
	EEG.event(index).eegtype  = typeeeg (EEG.event(index).epoch);
	EEG.event(index).response = response(EEG.event(index).epoch);
end;

for index = 1:n
	if rt(index) ~= no_rt
		EEG.event(end+1).type   = 'rt';
		EEG.event(end).latency  = eeg_lat2point(rt(index)/1000, index, EEG.srate, [EEG.xmin EEG.xmax]);
		EEG.event(end).epoch    = index;
		EEG.event(end).eegtype  = typeeeg(index);
		EEG.event(end).response = response(index);
	end
end;
tmpevent = EEG.event;
tmp = [ tmpevent.latency ];
[tmp indexsort] = sort(tmp);
EEG.event = EEG.event(indexsort);
EEG = eeg_checkset(EEG, 'eventconsistency');

command = sprintf('%s = pop_loaddat(%s, %s, %d);', inputname(1), inputname(1), fullFileName, no_rt); 

return;
