% pop_copyset() - Copy the current dataset into another dataset.
%
% Usage:
%   >> ALLEEG = pop_copyset(ALLEEG, index1); % pop-up
%   >> ALLEEG = pop_copyset(ALLEEG, index1, index2 );
%
% Inputs:
%   ALLEEG - array of dataset structure
%   index1  - input dataset number
%   index2  - index of dataset to copy into
%
% Inputs:
%   ALLEEG - array of dataset structure
%
% Note: this function performs ALLEEG(index2) = ALLEEG(index1);
%       with dataset checks
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeg_store(), pop_delset(), eeglab() 

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
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 

function [ALLEEG, com] = pop_copyset(ALLEEG, set_in, set_out);

com = '';
if nargin < 2
	help pop_copyset;
	return;
end;
if isempty(EEG(set_in).data)
    disp('Pop_copyset error: cannot copy empty dataset'); return;
end;    
if nargin < 1
	% which set to save
	% -----------------
	promptstr    = { 'Enter the destination dataset:' };
	inistr       = { int2str(set_in+1) };
	result       = inputdlg( promptstr, 'Copy dataset -- pop_copyset()', 1,  inistr);
	size_result  = size( result );
	if size_result(1) == 0 return; end;
	set_out   	 = eval( result{1} );
end;
ALLEEG = eeg_store(ALLEEG, eeg_retrieve(ALLEEG, set_in), set_out);

com = sprintf('%s = pop_copyset( %s, %d, %d);', inputname(1), inputname(1), set_in, set_out);
return;
