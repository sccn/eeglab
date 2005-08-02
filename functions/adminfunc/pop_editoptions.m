% pop_editoptions() - Edit memory-saving eeglab() options. These are stored in 
%                     a file 'eeg_options.m'. With no argument, pop up a window 
%                     to allow the user to set/unset these options. Store
%                     user choices in a new 'eeg_options.m' file in the 
%                     working directory.
%
% Usage: >> pop_editoptions;
%        >> pop_editoptions( 'key1', value1, 'key2', value2, ...);
%
% Graphic interface inputs:
%   "Precompute ICA activations" - [checkbox] If set, all the ICA activation
%                   time courses are precomputed (this requires more RAM). 
%                   Command line equivalent: option_computeica.
%   "Retain parent dataset" - [checkbox] If set, EEGLAB will retain parents of 
%                   new datasets in memory (requiring more RAM). 
%                   Command line equivalent: option_keepdataset.
%   "Store data in the .set or .dat file" - [checkbox] Set -> dataset data (EEG.data) are 
%                   saved in the EEG structure in the standard Matlab dataset (.set) file. 
%                   Unset -> The EEG.data are saved as a transposed stream of 32-bit 
%                   floats in a separate binary file. As of Matlab 4.51, the order 
%                   of the data in the binary file is as in the transpose of EEG.data 
%                   (i.e., as in EEG.data', frames by channels). This allows quick 
%                   reading of single channels from the data, e.g. when comparing 
%                   channels across datasets. The stored files have the extension 
%                   .dat instead of the pre-4.51, non-transposed .fdt. Both file types 
%                   are read by the dataset load function. Command line equivalent: 
%                   option_savematlab.
% Commandline keywords:
%   'option_computeica' - [0|1] If 1, compute the ICA component activitations and
%                   store them in a new variable. If 0, compute ICA activations
%                   only when needed (& only partially, if possible) and do not
%                   store the results).
%   'option_keepdataset' - [0|1] If 1, keep datasets in memory so that the user 
%                   can undo any EEGLAB operation by returning to the parent
%                   dataset.  The user may work on several datasets at a time.
%                   If 0, only one dataset is stored in memory, any changes 
%                   overwriting the current (EEG) dataset. 
%   NOTE: Turn OFF the options above when working with very large datasets or on 
%                   computers with limited memory.
%   'option_savematlab' - [0|1] If 1, datasets are saved as single Matlab .set files. 
%                   If 0, dataset data are saved in separate 32-bit binary float 
%                   .dat files.  See the corresponding GUI option above for details. 
% Outputs:
%   In the output workspace, variables 'option_computeica', 'option_keepdataset',
%   and 'option_savematlab'  are updated, and a new 'eeg_options.m' file may be
%   written to the working directory. The copy of 'eeg_options.m' placed in your 
%   working directory overwrites system defaults whenever EEGLAB operates in this
%   directory (assuming your working directory is in your MATLABPATH - see path()).
%   To adjust these options system-wide, edit the master "eeg_options.m" file in the
%   EEGLAB directory heirarchy.
%
% Author: Arnaud Delorme, SCCN / INC / UCSD, March 2002
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
% Revision 1.24  2005/08/02 01:45:39  arno
% dealing with newer option files
%
% Revision 1.23  2005/08/01 14:36:03  arno
% debug wrong alignment etc...
%
% Revision 1.22  2005/07/31 23:01:56  arno
% debuging option writing
%
% Revision 1.21  2004/11/22 17:30:23  scott
% edited help message - now tells the whole story
%
% Revision 1.20  2004/11/21 02:45:57  scott
% help message - new transposed binary file option
%
% Revision 1.19  2003/07/31 22:28:41  arno
% *** empty log message ***
%
% Revision 1.18  2003/04/10 17:31:03  arno
% header edit
%
% Revision 1.17  2002/12/04 19:13:19  arno
% debuging for MacOS6 Matlab 6 directories
%
% Revision 1.16  2002/11/15 01:38:52  scott
% same
%
% Revision 1.15  2002/10/23 15:01:50  arno
% isppc -> computer
%
% Revision 1.14  2002/09/26 14:00:44  scott
% help msg -sm
%
% Revision 1.13  2002/08/20 22:36:47  arno
% debug for windows
%
% Revision 1.12  2002/08/19 21:55:40  arno
% add Mac statement
%
% Revision 1.11  2002/08/14 21:30:48  arno
% debug for windows
%
% Revision 1.10  2002/08/13 23:47:11  arno
% debugging message
%
% Revision 1.9  2002/08/13 18:25:57  scott
% help msg
%
% Revision 1.8  2002/08/13 16:10:43  arno
% debugging
%
% Revision 1.7  2002/08/13 00:31:31  scott
% text
%
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
		if ~popask(['Cannot modify read-only file ''' filename '''' 10 'Should EEGLAB use a writable copy in the current directory ?']);
			return;
		else 
			fid = fopen( filename, 'r');
			if	fid == -1
				error('Cannot open file');
			end;
			storelocal = 1;
		end;
	else
		error('File not found');
	end;
end;

% read variables values and description
% --------------------------------------
fid2 = fopen('eeg_optionsbackup.m', 'r');
try 
    [ header varname1 value1 description1 ] = readoptionfile( fid  );
catch, varname1 = {}; value1 = {}; end;
[ header varname  value  description  ] = readoptionfile( fid2 );

% fuse the two informations
% -------------------------
varname(1:length(varname1)) = varname1;
value  (1:length(varname1)) = value1;

if nargin == 0
    geometry = { [4 1] };
    uilist = { ...
         { 'Style', 'text', 'string', 'Description', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'Set/Unset', 'fontweight', 'bold'   } };

    % add all fields to graphic interface
    % -----------------------------------
    for index = 1:length(varname)
        try, % format the description to fit a help box THIS DOES NOT WORK (CONFLICT BETWEEN WAITFOR & QUESTLDG IN INPUTDLG)
             % ----------------------------------------
            tmptext = description{ index };
            if length(tmptext) > 40,    stringtext = [ tmptext(1:40) '...' ]; 
            else                        stringtext = tmptext; 
            end;
            descrip = { 'string', stringtext, 'callback', ['questdlg2([''' ...
                                choptext( tmptext ) '''],''Description of field ' varname{index} ''', ''OK'', ''OK'');' ] }; 
        catch, descrip = { 'string', 'no-description' }; end;

        descrip = { 'string', description{ index } };
           
        % create the gui for this variable
        % --------------------------------
        geometry = { geometry{:} [4 0.7 0.45 0.5] };
        uilist   = { uilist{:}, ...
         { 'Style', 'text', descrip{:}, 'horizontalalignment', 'left' }, ...
         { }, { 'Style', 'checkbox', 'string', '    ', 'value', value{index} } { } }; 
    end;

    results = inputgui( geometry, uilist, 'pophelp(''pop_editoptions'');', 'Memory options - pop_editoptions()' );
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
    warning off;
	if isunix | strcmp(computer,'MAC')
        delimloc = sort(union(findstr(filename, ':'), findstr(filename, '/')));
	else
        delimloc = findstr(filename, '\');
	end;
    warning on;
	filename = filename(delimloc(end)+1:end);
end;
fid = fopen( filename, 'w');
if fid == -1
	error('File error, check writing permission');
end;
fprintf(fid, '%s\n', header);
for index = 1:2:length(args)
    fprintf( fid, '%s = %d ; %% %s\n', args{index}, args{index+1}, description{(index-1)/2+1});
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

% read option file
% ----------------
function [ header, varname, value, description ] = readoptionfile( fid );
    
    % skip header
    % -----------
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
        [ varname{index} str ] = strtok(str); % variable name
        [ equal          str ] = strtok(str); % =
        [ value{index}   str ] = strtok(str); % value
        [ tmp            str ] = strtok(str); % ;
        [ tmp            dsc ] = strtok(str); % comment
        dsc = deblank( dsc(end:-1:1) );
        description{index} = deblank( dsc(end:-1:1) );
        value{index}       = str2num( value{index} );
        
        str = fgets( fid ); % jump a line
        index = index+1;
    end;
    fclose(fid);
