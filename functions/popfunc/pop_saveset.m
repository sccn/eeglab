% pop_saveset() - save dataset(s) structure
%
% Usage:
%   >> pop_saveset( EEG ); % pop_up
%   >> pop_saveset( ALLEEG ); % pop_up
%   >> EEG = pop_saveset( EEG, filename, filepath);
%   >> ALLEEG = pop_saveset( ALLEEG, indices, filename, filepath);
%
% Inputs:
%   EEG        - dataset structure
%   ALLEEG     - array of dataset structure
%   indices    - indices of datasets to save from ALLEEG global
%                structure
%   filename   - name of the file to save (optional)
%   filepath   - path of the file to save (optional)
%
% Outputs:
%   EEG        - dataset after extended syntax check
%   ALLEEG     - datasets after extended syntax check
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
% Revision 1.4  2002/08/12 18:34:13  arno
% questdlg2
%
% Revision 1.3  2002/08/12 02:27:13  arno
% inputdlg2
%
% Revision 1.2  2002/04/23 21:45:48  arno
% making the function standalone
% ,
%
% Revision 1.1  2002/04/18 19:57:34  arno
% Initial revision
%
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

function [EEG, com] = pop_saveset( EEG, inarg, curfilename, curfilepath);

com = '';
if nargin < 1
	help pop_saveset;
	return;
end;

if length(EEG) > 1
	mode = 1; % multiple datasets
else
	mode = 0; % single datasets
end;

if nargin < 2
	% ask user
	if mode == 1 % multiple datasets
		indices = [];
		for index = 1:length(EEG)
			if ~isempty(EEG(index).data)
				indices = [indices index];
			end;
		end;
		result = inputdlg2( {'Indices of the datasets to save'}, 'Save several datasets -- pop_saveset()', 1, {int2str(indices)}, 'pop_saveset');
		if length(result) == 0 return; end;
		indices = eval( [ '[' result{1} ']' ] );
		[curfilename, curfilepath] = uiputfile('*.sets', 'Save dataset with .sets extension -- pop_saveset()'); 	
	else
		[curfilename, curfilepath] = uiputfile('*.set', 'Save dataset with .set extension -- pop_saveset()'); 
	end;
	if curfilename == 0 return; end;	
else 
	if mode == 1
		indices = inarg;
	else
		curfilepath = curfilename;
		curfilename = inarg;
	end;
end;

if mode == 0  % single datasets
	disp('Pop_saveset: extended datasets syntax check...');
	EEG = eeg_checkset(EEG, 'eventconsistency');
	EEG.filename    = curfilename;
	EEG.filepath    = curfilepath;
	tmpica = EEG.icaact;
	EEG.icaact      = [];
	disp('Saving dataset...');
	command = sprintf('save -MAT %s EEG', [ curfilepath curfilename ]);
	try, eval(command);
	catch, error('Pop_saveset: saving error, check permission on file or directory');
	end;
	EEG.icaact = tmpica;
	
	com = sprintf('EEG = pop_saveset( %s, ''%s'', ''%s'');', inputname(1), curfilename, curfilepath);
else
	ALLEEG = EEG; clear EEG;
	
	if max(indices) > length(ALLEEG)
		error('Pop_saveset: index out-of-bound');
	end;

	% checking datasets
	% -----------------
	disp('Pop_saveset: extended datasets syntax check...');
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
		else
			disp(['Pop_saveset warning: dataset ' int2str(indices(index)) ' is empty']);
		end;
	end;
	TMPALLEEG = ALLEEG;
	for index = 1:length(indices)
		ALLEEG(indices(index)).icaact = [];
	end;
	
	% saving
	% ------
	disp('Pop_saveset: saving datasets...');
	ALLEEG = ALLEEG(indices);
	try, save([ curfilepath curfilename ], '-mat', 'ALLEEG');
	catch, error('Pop_saveset: saving error, check permission on file or directory');
	end;
	ALLEEG = TMPALLEEG;
	EEG = ALLEEG;
	
	com = sprintf('ALLEEG = pop_saveset( %s, %s, ''%s'', ''%s'');', inputname(1), vararg2str(indices), curfilename, curfilepath);
end;
return;

function num = popask( text )
	 ButtonName=questdlg2( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end;


