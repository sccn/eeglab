% pop_saveset() - save dataset(s) structure
%
% Usage:
%   >> pop_saveset( EEG, filename, filepath);
%   >> pop_saveset( indices, filename, filepath);
%
% Inputs:
%   EEG        - data structure
%   indices    - indices of datasets to save from ALLEEG global
%                structure
%   filename   - name of the file to save (optional)
%   filepath   - path of the file to save (optional)
%
% Note: individual dataset should be save with '.set' file extension
%       and multiple datasets with '.sets' file extension
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: pop_loadset(), eeglab()
  
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
% Revision 1.3  2002/04/11 03:35:48  arno
% *** empty log message ***
%
% Revision 1.2  2002/04/11 02:04:53  arno
% adding further checks
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 

function com = pop_saveset( inarg, curfilename, curfilepath);
eeg_global;

if nargin > 1 & isnumeric(inarg)
	indices = inarg;
else
	indices = [];
	for index = 1:length(ALLEEG)
		if ~isempty(ALLEEG(index).data)
			indices = [indices index];
		end;
	end;
end;

if nargin < 2
	% ask user
	if nargin < 1 | isnumeric(inarg)
		result = inputdlg( {'Indices of the datasets to save'}, 'Save several datasets -- pop_saveset()', 1, {int2str(indices)});
		if length(result) == 0 return; end;
		indices = eval( [ '[' result{1} ']' ] );
	end;

	if length(indices) == 1 | (nargin >= 1 & isstruct(inarg))
		[curfilename, curfilepath] = uiputfile('*.set', 'Save dataset with .set extension -- pop_saveset()'); 	
	else 
		[curfilename, curfilepath] = uiputfile('*.sets', 'Save dataset with .sets extension -- pop_saveset()'); 	
	end;
	if curfilename == 0 return; end;
end;

if (nargin >= 1 & isstruct(inarg)) | length(indices) == 1
	if nargin >= 1 & isstruct(inarg)
		EEG = inarg;
	else
		if indices < length(ALLEEG.data) & ~isempty(ALLEEG(indices).data)
			EEG = eeg_retrieve(indices)	
		else 
			error('Pop_saveset: index out-of-bound');			
		end;
	end;
	EEG.filename    = curfilename;
	EEG.filepath    = curfilepath;
	EEG.icaact      = [];
	disp('Saving dataset...');
	command = sprintf('save -MAT %s EEG', [ curfilepath curfilename ]);
	try, eval(command);
	catch, error('Pop_saveset: saving error, check permission on file or directory');
	end;
	
	EEG = eeg_checkset(EEG);
	com = sprintf('pop_saveset( %s, ''%s'', ''%s'');', inputname(1), curfilename, curfilepath);
else
	if max(indices) > length(ALLEEG)
		error('Pop_saveset: index out-of-bound');
	end;
	
	% checking datasets
	% -----------------
	disp('Pop_saveset: checking datasets syntax...');
	for index = 1:length(indices)
		if ~isempty(ALLEEG(indices(index)).data)
			try, ALLEEG(indices(index)) = eeg_checkset(ALLEEG(indices(index)), 'eventconsistency');
			catch
				if nargin < 2
					if ~popask( [ 'Warning: dataset ' int2str(indices(index)) ' has non-consistent events' 10 ...
								  'Do you want to continue ?' ])
						error( ['dataset ' int2str(indices(index)) ' has non-consistent events. You must fix the problem.']);
					end;	
				else
					disp( ['Warning: dataset ' int2str(indices(index)) ' has non-consistent events. You must fix the problem.']);
				end;
			end;
			ALLEEG(indices(index)).icaact = [];
		else
			disp(['Pop_saveset warning: dataset ' int2str(indices(index)) ' is empty']);
		end;
	end;

	% saving
	% ------
	disp('Pop_saveset: saving datasets...');
	TMPSAVE = ALLEEG;
	ALLEEG = ALLEEG(indices);
	try, save([ curfilepath curfilename ], '-mat', 'ALLEEG');
	catch, error('Pop_saveset: saving error, check permission on file or directory');
	end;
	ALLEEG = TMPSAVE;

	% restoring
	% ---------
	for index = 1:length(indices)
		if ~isempty(ALLEEG(indices(index)).data)
			ALLEEG(indices(index)) = eeg_checkset(ALLEEG(indices(index)));
		end;
	end;

	com = sprintf('pop_saveset( %s, ''%s'', ''%s'');', vararg2str(indices), curfilename, curfilepath);
end;
return;

function num = popask( text )
	 ButtonName=questdlg( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end;


