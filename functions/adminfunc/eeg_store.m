% eeg_store() - store a dataset into the variable containing
%               all datasets
%
% Usage: >> [ALLEEG EEG index] = eeg_store(ALLEEG, EEG);
%        >> [ALLEEG EEG index] = eeg_store(ALLEEG, EEG, index);
%
% Inputs:
%   ALLEEG     - variable containing all datasets
%   EEG        - current dataset to store. EEG can also contain an 
%                array of dataset that will be appended to the current
%                array of dataset.
%   index      - (optional), index of where to store the new
%                dataset. If no index is given, the first 
%                empty slot in the ALLEEG array is selected.
%
% Outputs:
%   ALLEEG - variable containing all datasets
%   EEG    - EEG dataset after syntax checking 
%   index  - index of the new dataset
%
% Note: given 3 arguments, after checking the input dataset 
%       structure syntax consistency, this function simply perfroms 
%       >> ALLEEG(index) = EEG;
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeg_retrieve(), eeglab()

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

% uses the global variable EEG ALLEEG CURRENTSET 

% $Log: not supported by cvs2svn $
% Revision 1.10  2002/04/25 18:53:12  arno
% tabs
%
% Revision 1.9  2002/04/25 17:10:07  scott
% editting msg
%
% Revision 1.8  2002/04/25 00:21:54  scott
% improving wording -sm
%
% Revision 1.7  2002/04/23 23:48:34  arno
% standalone version
%
% Revision 1.6  2002/04/23 21:28:34  arno
% correcting typo
%
% Revision 1.5  2002/04/23 17:57:52  arno
% allowing store of multiple datasets
%
% Revision 1.4  2002/04/20 17:18:33  arno
% removing "Done."
%
% Revision 1.3  2002/04/18 20:01:34  arno
% retrIeve
%
% Revision 1.2  2002/04/18 14:52:31  scott
% added "Done." -Scott
%
% Revision 1.1  2002/04/05 17:46:04  jorn
% Initial revision
%
% 01-25-02 reformated help & license -ad 
% 03-07-02 add the eeglab options -ad

% store set
% ------------------
function [ALLEEG, EEG, storeSetIndex] = eeg_store(ALLEEG, EEG, storeSetIndex);

% considering multiple datasets
% -----------------------------
if length(EEG) > 1
	TMPEEG = EEG;
	for index=1:length(TMPEEG)
		EEG = TMPEEG(index);
		[ALLEEG, EEG, storeSetIndex] = eeg_store(ALLEEG, EEG);
	end;
	return;
end;
EEG = eeg_checkset(EEG);

% test options
% ------------
eeg_options; 
if option_keepdataset == 0
    disp('Current dataset changed.');
	storeSetIndex = 0;
    return;
end;    

% find first free index
% ---------------------
if nargin < 3
	i = 1;
	while (i<200)
		try
			if isempty(ALLEEG(i).data);
				storeSetIndex = i; i = 200;
			end;
			i = i+1;	
		catch
			storeSetIndex = i; i = 200;
		end;
   end;
   fprintf('Creating a new dataset with index %d\n', storeSetIndex);
else
	if isempty(storeSetIndex) | storeSetIndex == 0
		storeSetIndex = 1;
	end;
end;

if ~isempty( ALLEEG )
	try
		ALLEEG(storeSetIndex) = EEG;
	catch
		allfields = fieldnames( EEG );
		for i=1:length( allfields )
			eval( ['ALLEEG(' int2str(storeSetIndex) ').' allfields{i} ' = EEG.' allfields{i} ';' ]);
		end;	
	end;
else	
	ALLEEG = EEG;
	if storeSetIndex ~= 1
		ALLEEG(storeSetIndex+1) = EEG;
		ALLEEG(1) = ALLEEG(storeSetIndex); % empty
 		ALLEEG(storeSetIndex) = ALLEEG(storeSetIndex+1); 
 		ALLEEG = ALLEEG(1:storeSetIndex);
 	end;	
end;	
return;
