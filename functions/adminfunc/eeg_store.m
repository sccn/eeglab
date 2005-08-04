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
% Typical use:
% [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
% creates a new dataset in variable ALLEEG
% [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
% erase current dataset in variable ALLEEG
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

% check of parameter consistency
% ------------------------------
if nargin == 3
    if length(EEG) ~= length(storeSetIndex),
        error('Length of input dataset structure must be equal to length of index array');
    end;
end;

% considering multiple datasets
% -----------------------------
if length(EEG) > 1
	TMPEEG = EEG;
    if nargin >= 3
        for index=1:length(TMPEEG)
            EEG = TMPEEG(index);
            [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, storeSetIndex(index), varargin{:});
            TMPEEG(index) = EEG;
        end;        
    else
        for index=1:length(TMPEEG)
            EEG = TMPEEG(index);
            [ALLEEG, EEG, storeSetIndex(index)] = eeg_store(ALLEEG, EEG);
            TMPEEG(index) = EEG;
        end;
    end;
    EEG = TMPEEG;
	return;
end;
if nargin < 4 
    [ EEG com ]  = eeg_checkset(EEG);
    EEG = ALLEEG(CURRENTSET);
    EEG.changes_not_saved = 'yes';
else % savedata
     % --------
    eeg_optionsbackup;
    eeg_options;
    if option_storedisk & strcmpi(EEG.changes_not_saved, 'yes')
        if option_warningstore & strcmpi(varargin{1}, 'savegui')
            if ~isempty(EEG.filename) & ~isempty(EEG.filepath)                
                res = questdlg2(strvcat( 'You set the option to keep at most one dataset in memory at a time, so the previous', ...
                                         'dataset are automatically saved on disk (and overwrite previous dataset file).', ...
                                         'Do you still wish to proceed?', ...
                                         '(Use menu item "File > Maximize memory" to disable this functionality or', ...
                                         'to prevent this warning from appearing.)'), 'Make exception for this dataset', ...
                                         'Save as new file', 'Yes resave previous datase');
                if strcmpi(res, 'Yes resave previous datase'), option_save = 'resave';
                elseif strcmpi(res, 'Save as new file'),       option_save = 'new';
                else                                           option_save = 'exception';
                end;
            else
                res = questdlg2(strvcat( 'You set the option to keep at most one dataset in memory at a time, so the previous', ...
                                         'dataset has to be saved on disk. Do you still wish to proceed?', ...
                                         'If yes, you will be prompted to enter a file name.', ...
                                         '(Use menu item "File > Maximize memory" to disable this functionality or', ...
                                         'to prevent this warning from appearing.)'), 'Make exception for this dataset', ...
                                         'Yes save previous dataset');                
                if strcmpi(res, 'Yes save previous datase'), option_save = 'new';
                else                                         option_save = 'exception';
                end;
            end;
            
            % save mode
            % ---------
            if strcmpi(option_save, 'new')
                EEG.changes_not_saved = 'no';
                [ EEG com] = pop_saveset(EEG);
                EEG = update_datafield(EEG);
                h(com);
                if isempty(com)
                    disp('No file name given. Dataset not saved (making an exception).');
                    option_save = 'exception';
                end;
            end;
            if strcmpi(option_save, 'exception')
                EEG.changes_not_saved = 'yes';
                [ EEG com ] = eeg_checkset(EEG);
                EEG = ALLEEG(CURRENTSET);
            else
                EEG.changes_not_saved = 'no';
                [ EEG com] = pop_saveset(EEG, 'savemode', 'resave');
                EEG = update_datafield(EEG);
            end;
        else
            EEG.changes_not_saved = 'no';
            [ EEG com] = pop_saveset(EEG, 'savemode', 'resave');
            EEG = update_datafield(EEG);
        end;
    else
        if strcmpi(EEG.changes_not_saved, 'no')
            disp('Dataset not modified since last save, no need for resaving it');
        end;
        [ EEG com ] = eeg_checkset(EEG);
        EEG = ALLEEG(CURRENTSET);
        EEG.changes_not_saved = 'yes';
    end;
end;
EEG = eeg_hist(EEG, com);

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

function EEG = update_datafield(EEG);
    if isfield(EEG, 'datfile')
        if ~isempty(EEG.datfile)
            EEG.data = EEG.datfile;
        else
            EEG = rmfield(EEG, 'datfile');
        end;
    else 
        EEG.data = 'in set file';
    end;
