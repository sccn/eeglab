% pop_editoptions() - Edit eeglab() options stored in the eeg_options()
%                    Matlab file. With no argument, a window pops-up to
%                    query the user about option values.
%
% Usage: >> pop_editoptions;
%        >> pop_editoptions( 'key1', value1, 'key2', value2, ...);
%
% Optional inputs:
%   'option_computeica' - [0|1] If 1, compute the ICA component activity and
%                   store it into a new variable. If 0, compute ICA activations
%                   only when needed (only partially, if possible) and do not 
%                   store the result). 0 may be used to process large datasets.
%   'option_keepdataset' - [0|1]. If 1, keep datasets so that the user can undo 
%                   any EEGLAB operations by returning to previous datasets.
%                   The user may also work on several datasets at a time.
%                   If 0, only one dataset is stored in memory, replacing the
%                   input (EEG) dataset. 
% Outputs:
%   In the output workspace, variables 'option_computeica', 
%   'option_keepdataset' and 'option_usedisk' are updated.
%
% Note:
%   Keep a copy of 'editoptions.m' in your working directory to overwrite system
%   defaults (assuming that '.' is the first non-Matlab directory in your 
%   MATLABPATH).
%
% Author: Arnaud Delorme, CNL / Salk Institute, 09 March 2002
%
% See also: eeg_options()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 09 March 2002, arno@salk.edu
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
% Revision 1.6  2002/08/12 18:40:17  arno
% questdlg2
%
% Revision 1.5  2002/07/25 17:22:31  arno
% adding a clear functions statement
%
% Revision 1.4  2002/04/25 16:56:13  arno
% copying file if read only
%
% Revision 1.3  2002/04/21 01:07:59  scott
% edited help msg -sm
%
% Revision 1.2  2002/04/18 18:09:57  arno
% updating error message
%
% Revision 1.1  2002/04/05 17:46:04  jorn
% Initial revision
%
%02/19/2002 debuging function -ad

function com = pop_editoptions(varargin);

com = '';
% parse the eeg_options file
% ----------------------------
filename = which('eeg_options.m');
fid = fopen( filename, 'r+');
storelocal = 0;
if	fid == -1
	if exist(filename) == 2 
		if ~popask(['Can not modify read-only file ''' filename '''' 10 'do you want EEGLAB to store a copy in the current directory ?']);
			return;
		else 
			fid = fopen( filename, 'r');
			if	fid == -1
				error('Can not open file');
			end;
			storelocal = 1;
		end;
	else
		error('File not found');
	end;
end;

% store header
% ------------
header = '';
str = fgets( fid );
while (str(1) == '%')
    header = [ header str];
    str = fgets( fid );
end;

% read variables values and description
% --------------------------------------
str = fgetl( fid ); % jump a line
index = 1;
while (str(1) ~= -1)
    [varname{index} tmp1 tmp2 count1] = sscanf(str, '%s',1);
    [equal          tmp1 tmp2 count2] = sscanf(str(count1:end), '%s',1);
    [value{index}   tmp1 tmp2 count3] = sscanf(str((count1+count2):end), '%d',1);
    description{index} = str((count1+count2+count3+3):end-1);

    str = fgets( fid ); % jump a line
    index = index+1;
end;
fclose(fid);

if nargin == 0
    geometry = { [4 1] };
    uilist = { ...
         { 'Style', 'text', 'string', 'Desciption', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'Set/Unset', 'fontweight', 'bold'   } };

    % add all fields to graphic interface
    % -----------------------------------
    for index = 1:length(varname)
        try, % format the description to fit a help box THIS DOES NOT WORK (CONFLICT BETWEEN WAITFOR & QUESTLDG IN INPUTDLG)
             % ----------------------------------------
            tmptext = description{ index };
            if length(tmptext) > 25,    stringtext = [ tmptext(1:25) '...' ]; 
            else                        stringtext = tmptext; 
            end;
            descrip = { 'string', stringtext, 'callback', ['questdlg2([''' choptext( tmptext ) '''],''Description of field ' varname{index} ''', ''OK'', ''OK'');' ] }; 
        catch, descrip = { 'string', 'no-description' }; end;

        descrip = { 'string', description{ index } };
           
        % create the gui for this variable
        % --------------------------------
        geometry = { geometry{:} [4 0.7 0.35 0.5] };
        uilist   = { uilist{:}, ...
         { 'Style', 'text', descrip{:}, 'horizontalalignment', 'left' }, ...
         { }, { 'Style', 'checkbox', 'string', '    ', 'value', value{index} } { } }; 
    end;

    results = inputgui( geometry, uilist, 'pophelp(''editeegoptions'');' );
    if length(results) == 0, return; end;
   
    % decode inputs
    % -------------
    args = {};
    for index = 1:length(varname)
        args = {  args{:}, varname{index}, results{index} }; 
    end;
else % no interactive inputs
    args = varargin;
    for index = 1:2:length(varargin)
        if isempty(strmatch(varargin{index}, varname, 'exact'))
            error(['Variable name ''' varargin{index} ''' is invalid']);
        end;
    end;        
end;

% write to eeg_options file
% -------------------------
if storelocal
	delimloc = findstr(filename, '/');
	filename = filename(delimloc(end)+1:end);
end;
fid = fopen( filename, 'w');
if fid == -1
	error('File not found');
end;
fprintf(fid, '%s\n', header);
for index = 1:2:length(args)
    fprintf( fid, '%s = %d ;%% %s\n', args{index}, args{index+1}, description{(index-1)/2+1});
end;
fclose(fid);    

% generate the output text command
% --------------------------------
com = 'editeegoptions(';
for index = 1:2:length(args)
    com = sprintf( '%s ''%s'', %d,', com, args{index}, args{index+1});
end;
com = [com(1:end-1) ');'];   
clear functions

% ---------------------------
function  chopedtext = choptext( tmptext )
    chopedtext = '';
    while length(tmptext) > 30
          blanks = findstr( tmptext, ' ');
          [tmp I] = min( abs(blanks - 30) );
          chopedtext = [ chopedtext ''' 10 ''' tmptext(1:blanks(I)) ];
          tmptext  = tmptext(blanks(I)+1:end);
    end;    
    chopedtext = [ chopedtext ''' 10 ''' tmptext];
    chopedtext = chopedtext(7:end);
return;

function num = popask( text )
	 ButtonName=questdlg2( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end;
