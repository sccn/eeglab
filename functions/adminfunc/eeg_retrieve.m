% eeg_retrieve() - Retrieve an EEG dataset from the global variable
%                 containing all datasets
%
% Usage: >> EEG = eeg_retrieve( index );
%
% Inputs:
%   index      - index of the dataset to retrieve
%
% Outputs:
%   EEG        - output dataset. The global variable EEG is also updated 
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeg_store(), eeglab()

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
% Revision 1.1  2002/04/05 17:46:04  jorn
% Initial revision
%
% 01-25-02 reformated help & license -ad 

function EEG = eeg_retrieve( retrieveSetIndex);
eeg_global;
eeg_consts;

if nargin < 1
	error('Not enought arguments');
	return;
end;	

try
	EEG = ALLEEG(retrieveSetIndex);
catch
	fprintf('Warning: cannot retrieve dataset with index %d\n', retrieveSetIndex); 
	return;
end;

CURRENTSET = retrieveSetIndex;

return;

