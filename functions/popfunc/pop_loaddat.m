% pop_loaddat() - merge a neuroscan DAT file with input dataset
%                (pop out window if no arguments).
%
% Usage:
%   >> [ OUTEEG ] = pop_loaddat( INEEG, filename);
%
% Inputs:
%   filename       - file name
%   INEEG          - input EEGLAB data structure
%
% Outputs:
%   OUTEEG         - EEGLAB data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: loaddat(), eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $

% 01-25-02 reformated help & license -ad 
% 13/02/02 removed the no latency option -ad

function [EEG, command] = pop_loaddat(EEG, filename); 
command = '';

if nargin < 1
	help  pop_loaddat;
	return;
end;	

if nargin < 2 
	% ask user
	[filename, filepath] = uigetfile('*.DAT', 'Choose a DAT file -- pop_loaddat'); 
	if filename == 0 return; end;
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
	error('pop_loaddat, number of trials in input dataset and DAT file different, abording');
end;	   

EEG.epoch = eeg_epochformat([rt(:) typeeeg(:) response(:) ], 'struct', {'eventlatency', 'type', 'response'});
EEG = eeg_checkset(EEG);

command = sprintf('%s = pop_loaddat(''%s'', %s);', inputname(1), inputname(1), fullFileName); 

return;
