% eeg_store() - Store the current dataset (global variable)
%              into an global variable array
%
% Usage: >> index = eeg_store(index);
%
% Inputs:
%   index  - index for storing the dataset. If no index is given,
%            the first free index is used.
%
% Outputs:
%   EEG    - the global name of the current EEG dataset is updated 
%   index  - if automatically generated, the function returns this
%            index.
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
function storeSetIndex = eeg_store( storeSetIndex);
eeg_global;

% considering multiple datasets
% -----------------------------
if length(EEG) > 1
	TMPEEG = EEG;
	for index=1:length(TMPEEG)
		EEG = TMPEEG(index);
		eeg_store;
	end;
	return;
end;
EEG = eeg_checkset(EEG);

% test options
% ------------
eeg_options; % changed from eeglaboptions 3/30/02 -sm
if option_keepdataset == 0
    disp('Old dataset overwrote (see eeg_options.m file to change feature)');
    return;
end;    

% find first free index
% ---------------------
if nargin < 1
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
CURRENTSET = storeSetIndex;
return;
