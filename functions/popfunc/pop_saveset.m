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
%   indices    - indices of datasets to save from ALLEEG
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
% Revision 1.15  2002/11/14 23:03:04  arno
% debugging save from the command line
%
% Revision 1.14  2002/10/15 23:42:31  arno
% error if no delimiter for last character of directory
%
% Revision 1.13  2002/10/15 16:59:57  arno
% drawnow for windows
%
% Revision 1.12  2002/10/14 00:46:31  arno
% debugging .sets filename check
%
% Revision 1.11  2002/10/10 16:43:26  arno
% debugging fdt files
%
% Revision 1.10  2002/10/03 16:26:27  arno
% debugging extension
%
% Revision 1.9  2002/09/23 16:43:48  arno
% [Aimplementing saving as float
%
% Revision 1.8  2002/08/19 22:09:28  arno
% debugging save for MAC
%
% Revision 1.7  2002/08/14 00:12:52  arno
% new error message
%
% Revision 1.6  2002/08/14 00:11:28  arno
% update header
%
% Revision 1.5  2002/08/14 00:06:16  arno
% empty command as default
%
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
if isempty(EEG)
	error('Cannot save multiple datasets');
end;

if length(EEG) > 1
	mode = 1; % multiple datasets
else
	mode = 0; % single datasets
end;

if (nargin < 2 & mode == 0) | (nargin < 3 & mode == 1)
	% ask user
	if mode == 1 % multiple datasets
		indices = [];
		for index = 1:length(EEG)
			if ~isempty(EEG(index).data)
				indices = [indices index];
			end;
		end;
		result = inputdlg2( {'Indices of the datasets to save'}, 'Save several datasets -- pop_saveset()', 1, {int2str(indices)}, 'pop_saveset');
        drawnow;
		if length(result) == 0 return; end;
		indices = eval( [ '[' result{1} ']' ] );
		[curfilename, curfilepath] = uiputfile('*.sets', 'Save dataset with .sets extension -- pop_saveset()'); 	
	else
		[curfilename, curfilepath] = uiputfile('*.set', 'Save dataset with .set extension -- pop_saveset()'); 
	end;
    drawnow;
	if curfilename == 0 return; end;	
else 
	if mode == 1
		indices = inarg;
        if nargin < 4
            curfilepath = '';
        end;
	else
        if nargin < 3
            curfilepath ='';
        else 
            curfilepath = curfilename;
        end;
		curfilename = inarg;
	end;
end;

if ~isempty(curfilepath)
    if curfilepath(end) ~= ':' & curfilepath(end) ~= '/' & curfilepath(end) ~= '\'
        error('Last character of filepath must be a directory delimiter');
    end;
end;

% currentfilename without the .set
% --------------------------------
if mode == 0 & ( length(curfilename)<=3 | ~strcmp(lower(curfilename(end-3:end)), '.set'))
	disp('Adding ''.set'' extension to the file');
	curfilename = [ curfilename '.set' ];
end;
if mode == 1 & ( length(curfilename)<=4 | ~strcmp(lower(curfilename(end-4:end)), '.sets'))
	if length(curfilename)>3 & strcmp(lower(curfilename(end-3:end)), '.set')
		disp('Changing file extension to ''.sets''');
		curfilename = [ curfilename 's' ];
	else
		disp('Adding ''.sets'' extension to the file');
		curfilename = [ curfilename '.sets' ];
	end;
end;

if length(curfilename)>4
	if strcmp(lower(curfilename(end-3:end)), '.set') 
		noextcurfilename = curfilename(1:end-4);
	elseif strcmp(lower(curfilename(end-4:end)), '.sets')
		noextcurfilename = curfilename(1:end-5);
	else
		noextcurfilename = curfilename;
	end;
else 
	noextcurfilename = curfilename;
end;

if mode == 0  % single datasets
	disp('Pop_saveset: extended datasets syntax check...');
	EEG = eeg_checkset(EEG, 'eventconsistency');
	EEG.filename    = curfilename;
	EEG.filepath    = curfilepath;
	tmpica = EEG.icaact;
	EEG.icaact      = [];
	disp('Saving dataset...');
	command = sprintf('save(''%s'',''-MAT'',''EEG'');', [ curfilepath curfilename ]);
	
	% saving data as float or as Matlab
	eeg_options;
	if exist('option_savematlab') == 1 & option_savematlab == 0
		tmpdata = EEG.data;
		EEG.data = [ noextcurfilename '.fdt' ];
		try, 
			eval(command);
			floatwrite( tmpdata, [EEG.filepath EEG.data], 'ieee-le');
		catch, 
			error('Pop_saveset: saving error, check permission on file or directory');
		end;
		EEG.data = tmpdata; 
		EEG.icaact = tmpica;
	else % saving data as a single Matlab file 
		try, 
			eval(command);
		catch, 
			error('Pop_saveset: saving error, check permission on file or directory');
		end;
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
	
	% saving data as float or as Matlab
	eeg_options;
	if exist('option_savematlab') == 1 & option_savematlab == 0
		for index = 1:length(ALLEEG)
			tmpdata = ALLEEG(index).data;
			ALLEEG(index).data = [ noextcurfilename '.fdt' int2str(index) ];
			try, 
				floatwrite( tmpdata, [ curfilepath ALLEEG(index).data], 'ieee-le');
			catch, 
				error('Pop_saveset: saving error, check permission on file or directory');
			end;
		end;
	end;
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


