% eeg_store() - store specified EEG dataset(s) in the ALLEG variable 
%               containing all current datasets, after first checking 
%               dataset consistency using eeg_checkset().
%
% Usage: >> [ALLEEG EEG index] = eeg_store(ALLEEG, EEG);
%        >> [ALLEEG EEG index] = eeg_store(ALLEEG, EEG, index);
%
% Inputs:
%   ALLEEG     - variable containing all current EEGLAB datasets
%   EEG        - dataset(s) to store - usually the current dataset. 
%                May also be an array of datasets; these will be 
%                checked and stored separately in ALLEEG.
%   index      - (optional), ALLEEG index (or indices) to use to store 
%                the new dataset(s). If no index is given, eeg_store() 
%                uses the lowest empty slot(s) in the ALLEEG array. 
% Outputs:
%   ALLEEG - array of all current datasets
%   EEG    - EEG dataset (after syntax checking)
%   index  - index of the new dataset
%
% Note: When 3 arguments are given, after checking the consistency of 
%       the input dataset structures this function simply performs 
%        >> ALLEEG(index) = EEG;
%
% Typical use:
%        >> [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
%  creates a new dataset in variable ALLEEG.
%        >> [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
%  overwrites the current dataset in variable ALLEEG.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab(), eeg_checkset(), eeg_retrieve()

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
% Revision 1.45  2006/05/07 18:06:27  arno
% cosmetic changes
%
% Revision 1.44  2006/04/19 14:18:59  arno
% multiple dataset file
%
% Revision 1.43  2006/03/12 16:58:55  arno
% store new dataset
%
% Revision 1.42  2006/03/10 21:39:40  arno
% debug last change
%
% Revision 1.41  2006/03/10 21:36:58  arno
% storing datasets
%
% Revision 1.40  2006/02/03 00:27:02  arno
% remove dat file reference in dataset when creating new dataset
%
% Revision 1.39  2006/02/01 06:54:37  arno
% remove all the savedata stuff
%
% Revision 1.38  2006/01/31 20:50:41  arno
% eeglab options
%
% Revision 1.37  2006/01/26 18:53:15  scott
% looking at help & msgs  -sm
%
% Revision 1.36  2006/01/26 00:29:57  arno
% message
%
% Revision 1.35  2006/01/26 00:26:25  arno
% better message
%
% Revision 1.34  2005/11/02 18:35:04  arno
% h -> eegh
%
% Revision 1.33  2005/09/13 16:34:56  arno
% add a clear function
%
% Revision 1.32  2005/09/08 21:53:29  arno
% same
%
% Revision 1.31  2005/09/08 21:37:17  arno
% fix storing multiple datasets
% /
%
% Revision 1.30  2005/09/08 16:54:41  arno
% fixed saved field
%
% Revision 1.29  2005/08/16 17:35:49  scott
% EEG.changes_not_saved 'no'   ->   EEG.saved 'yes'   -sm
%
% Revision 1.28  2005/08/16 17:21:06  scott
% edited help message and all disp fprint and detailed text messages  -sm
%
% Revision 1.27  2005/08/08 18:10:14  arno
% fix various problems
%
% Revision 1.26  2005/08/08 17:42:16  arno
% comments
%
% Revision 1.25  2005/08/08 17:40:34  arno
% erasing file information for newly created datasets
%
% Revision 1.24  2005/08/04 23:37:58  arno
% new pop-up
%
% Revision 1.23  2005/08/04 17:21:44  arno
% only set changes not saved if updating dataset
%
% Revision 1.22  2005/08/04 17:09:50  arno
% remving ALLEEG stuff
%
% Revision 1.21  2005/08/04 16:28:28  arno
% reprogrammed function
%
% Revision 1.20  2005/08/04 15:35:03  arno
% remove option to keep only one dataset
%
% Revision 1.19  2005/08/04 15:28:40  arno
% change command to automatically save dataset
%
% Revision 1.18  2005/08/02 16:47:40  arno
% savedata
%
% Revision 1.17  2005/08/01 22:43:03  arno
% eeg_options call
%
% Revision 1.16  2005/07/28 18:11:20  arno
% allowing storing multiple datasets
%
% Revision 1.15  2003/12/05 20:07:48  arno
% eeg_hist
%
% Revision 1.14  2003/12/05 00:42:18  arno
% debug history
%
% Revision 1.13  2003/12/05 00:26:20  arno
% [6~[6~dataset history
%
% Revision 1.12  2002/07/24 00:51:47  arno
% editing header
%
% Revision 1.11  2002/04/25 19:05:12  arno
% default output for single dataset option
%
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
function [ALLEEG, EEG, storeSetIndex] = eeg_store(ALLEEG, EEG, storeSetIndex, varargin);

% check parameter consistency
% ------------------------------
if nargin == 3
    if length(EEG) ~= length(storeSetIndex) & storeSetIndex(1) ~= 0
        error('Length of input dataset structure must equal the length of the index array');
    end;
end;

% considering multiple datasets
% -----------------------------
if length(EEG) > 1
	TMPEEG = EEG;
    if nargin >= 3 & storeSetIndex(1) ~= 0
        for index=1:length(TMPEEG)
            EEG = TMPEEG(index);
            tmpsaved      = EEG.saved;
            [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, storeSetIndex(index), varargin{:});
            ALLEEG(storeSetIndex(index)).saved = tmpsaved;
            TMPEEG(index).saved                = tmpsaved;
        end;        
    else
        for index=1:length(TMPEEG)
            EEG = TMPEEG(index);
            tmpsaved      = EEG.saved;
            if strcmpi(tmpsaved, 'justloaded'), tmpsaved = 'yes'; end;
            [ALLEEG, EEG, storeSetIndex(index)] = eeg_store(ALLEEG, EEG);
            ALLEEG(storeSetIndex(index)).saved = tmpsaved;
            TMPEEG(index).saved                = tmpsaved;
        end;
    end;
    EEG = TMPEEG;
	return;
end;

if nargin < 3
    % creating new dataset
    % -> erasing file information
    % ---------------------------
    EEG.filename = '';
    EEG.filepath = '';
    if isfield(EEG, 'datfile')
        EEG = rmfield(EEG, 'datfile');
    end;
end;

[ EEG com ]  = eeg_checkset(EEG);
if nargin > 2, 
    if storeSetIndex == 0
        EEG.saved = 'yes'; % just loaded
    else 
        EEG.saved = 'no';
    end;
else
    EEG.saved = 'no';        
end;
EEG = eeg_hist(EEG, com);

% find first free index
% ---------------------
findindex = 0;
if nargin < 3,             findindex = 1;
elseif storeSetIndex == 0, findindex = 1; 
end;

if findindex
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
   if isempty(varargin)
       fprintf('Creating a new ALLEEG dataset %d\n', storeSetIndex);
   end;
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
        if ~isfield(EEG, 'datfile') & isfield(ALLEEG, 'datfile')
            ALLEEG(storeSetIndex).datfile = '';
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
