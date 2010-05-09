% pop_copyset() - Copy the current EEG dataset into another dataset.
%
% Usage:
%   >> ALLEEG = pop_copyset(ALLEEG, index1); % pop-up
%   >> [ ALLEEG EEG CURRENTSET ] = pop_copyset(ALLEEG, index1, index2 );
%
% Inputs:
%   ALLEEG     - array of dataset structure
%   index1     - input dataset number
%   index2     - index of dataset to copy into
%
% Inputs:
%   ALLEEG     - array of dataset structures
%   EEG        - new copied structure
%   CURRENTSET - index of the new dataset
%
% Note: this function performs ALLEEG(index2) = ALLEEG(index1);
%       with dataset checks
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeg_store(), pop_delset(), eeglab() 

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

function [ALLEEG, EEG, CURRENTSET, com] = pop_copyset(ALLEEG, set_in, set_out);

com = '';
if nargin < 2
	help pop_copyset;
	return;
end;
if isempty(ALLEEG)
	error(['Pop_copyset error: cannot copy' 10 'single dataset mode']);
end;	
if set_in == 0
    error('Pop_copyset error: cannot copy dataset'); return;
end;
if isempty(ALLEEG(set_in).data)
    error('Pop_copyset error: cannot copy empty dataset'); return;
end;

if nargin < 3
	% which set to save
	% -----------------
	promptstr    = { 'Index of the new dataset:'};
	inistr       = { int2str(set_in+1) };
	result       = inputdlg2( promptstr, 'Copy dataset -- pop_copyset()', 1,  inistr, 'pop_copyset');
	if size( result ) == 0, EEG = []; CURRENTSET = 0; return; end;
	set_out   	 = eval( result{1} );
end;
ALLEEG = eeg_store(ALLEEG, eeg_retrieve(ALLEEG, set_in), set_out);
EEG    = eeg_retrieve(ALLEEG, set_out);
CURRENTSET = set_out;

com = sprintf('[ALLEEG EEG CURRENTSET] = pop_copyset( %s, %d, %d);', inputname(1), set_in, set_out);
return;
