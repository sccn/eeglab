% pop_copyset() - Copy the current EEG dataset into another dataset.
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
%   ALLEEG - array of dataset structures
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
% Revision 1.6  2002/08/13 23:57:45  arno
% update error message
%
% Revision 1.5  2002/08/13 17:22:39  arno
% text
%
% Revision 1.4  2002/08/12 16:20:06  arno
% debug
%
% Revision 1.3  2002/08/12 16:12:32  arno
% updating messages
%
% Revision 1.2  2002/04/23 22:05:43  arno
% making the function standalone
% ,
%
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
if isempty(ALLEEG)
	error(['Pop_copyset error: cannot copy' 10 'dataset in single dataset mode']);
end;	
if isempty(ALLEEG(set_in).data)
    error('Pop_copyset error: cannot copy empty dataset'); return;
end;
if set_in == 0
    error('Pop_copyset error: cannot copy dataset'); return;
end;

if nargin < 3
	% which set to save
	% -----------------
	promptstr    = { 'Index of the new dataset:'};
	inistr       = { int2str(set_in+1) };
	result       = inputdlg2( promptstr, 'Copy dataset -- pop_copyset()', 1,  inistr, 'pop_copyset');
	size_result  = size( result );
	if size_result(1) == 0 return; end;
	set_out   	 = eval( result{1} );
end;
ALLEEG = eeg_store(ALLEEG, eeg_retrieve(ALLEEG, set_in), set_out);

com = sprintf('%s = pop_copyset( %s, %d, %d);', inputname(1), inputname(1), set_in, set_out);
return;
