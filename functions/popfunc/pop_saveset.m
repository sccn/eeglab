% pop_saveset() - save one or more EEG dataset structures
%
% Usage:
%   >> pop_saveset( EEG ); % use an interactive pop-up window 
%   >> EEG = pop_saveset( EEG, 'key', 'val', ...); % no pop-up
%                                              
% Inputs:
%   EEG        - EEG dataset structure. May only contain one dataset.
%
% Optional inputs:
%   'filename' - [string] name of the file to save to
%   'filepath' - [string] path of the file to save to
%   'check'    - ['on'|'off'] perform extended syntax check. Default 'off'.
%   'savemode' - ['resave'|'onefile'|'twofiles'] 'resave' resave the 
%                current dataset using the filename and path stored
%                in the dataset; 'onefile' saves the full EEG 
%                structure in a Matlab '.set' file, 'twofiles' saves 
%                the structure without the data in a Matlab '.set' file
%                and the transposed data in a binary float '.dat' file.
%                By default the option from the eeg_options.m file is 
%                used.
%
% Outputs:
%   EEG        - saved dataset (after extensive syntax checks)
%   ALLEEG     - saved datasets (after extensive syntax checks)
%
% Note: An individual dataset should be saved with a '.set' file extension,
%       multiple datasets (at once) with a '.sets' file extension
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
% Revision 1.50  2005/07/29 22:55:30  arno
% saving several datasets
%
% Revision 1.49  2005/07/07 22:26:54  hilit
% cleared the 'savemode' help message
%
% Revision 1.48  2005/04/01 01:02:43  hilit
% fixing a bug
%
% Revision 1.47  2005/04/01 00:56:20  hilit
% disable data loading in savemode option 2
%
% Revision 1.46  2005/03/31 23:40:16  hilit
% added another option to savemode, saves the data structure without the data
%
% Revision 1.45  2005/02/16 19:42:51  hilit
% in save changing '-V6' to '-v6'
%
% Revision 1.44  2004/12/17 17:34:02  hilit
% correcting 'savemode' help message
%
% Revision 1.43  2004/12/16 22:33:02  arno
% removing EEGDATA
%
% Revision 1.42  2004/12/16 19:29:58  arno
% new option savemode
%
% Revision 1.41  2004/12/03 23:55:38  hilit
% uiputfile -> uiputfile2
%
% Revision 1.40  2004/11/15 22:49:29  arno
% saving in .dat files
%
% Revision 1.39  2004/11/05 19:28:15  arno
% remove uiputfile2
%
% Revision 1.38  2004/11/05 19:27:34  arno
% uiputfile -> uiputfile2
%
% Revision 1.37  2004/09/22 21:29:56  hilit
% save EEG.filepath as current path
%
% Revision 1.36  2004/09/14 17:04:32  arno
% saving 2 separate ariables
%
% Revision 1.35  2004/08/18 21:30:16  arno
% debug last
%
% Revision 1.34  2004/08/18 21:26:47  hilit
% removed asdf bug
%
% Revision 1.33  2004/08/18 21:21:17  arno
% message
%
% Revision 1.32  2004/07/30 17:03:23  arno
% same
%
% Revision 1.31  2004/07/30 17:02:17  arno
% debug last
%
% Revision 1.30  2004/07/30 16:57:26  arno
% allow to save under Matlab 7
%
% Revision 1.29  2003/12/05 18:50:27  arno
% remove trailing blanks
%
% Revision 1.28  2003/10/31 18:42:01  arno
% removing extra spaces
%
% Revision 1.27  2003/10/31 01:12:57  scott
% same
%
% Revision 1.26  2003/10/31 01:11:40  scott
% message text
%
% Revision 1.25  2003/10/31 00:52:30  scott
% commandline messages
%
% Revision 1.24  2003/07/24 17:54:03  arno
% changing message
%
% Revision 1.23  2003/07/24 16:15:55  arno
% deleting old .fdt files
%
% Revision 1.22  2003/07/22 15:43:10  arno
% *** empty log message ***
%
% Revision 1.21  2003/05/29 18:54:59  arno
% debug command line call
%
% Revision 1.20  2003/03/05 19:47:15  arno
% adding done message
%
% Revision 1.19  2003/03/03 21:35:59  arno
% correcting save location problem
%
% Revision 1.18  2003/02/26 02:20:07  scott
% header edits -sm
%
% Revision 1.17  2003/02/26 02:12:42  arno
% always forcing filepath to []
%
% Revision 1.16  2002/11/15 01:40:45  scott
% Can not -> cannot
%
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

function [EEG, com] = pop_saveset( EEG, varargin);

com = '';
if nargin < 1
	help pop_saveset;
	return;
end;
if isempty(EEG)  , error('Cannot save empty datasets'); end;
if length(EEG) >1, error('For reasons of consistency, this function  does not save multiple datasets any more'); end;

if nargin < 2
    % pop up window to ask for file type
    % ----------------------------------
    [filename, filepath] = uiputfile2('*.set', 'Save dataset with .set extension -- pop_saveset()'); 
    drawnow;
    options = { 'filename' filename 'filepath' filepath };
else
    % account for old calling format
    % ------------------------------
   if strcmpi(inputname, 'filename'),
        options = { 'filename' varargin{1} };
        if nargin > 1
            options = { options{:} 'filepath' varargin{2} };
        end;
    else
        options = varargin;
    end;
end;

% decode input parameters
% -----------------------
g = finputcheck(options,  { 'filename'   'string'   []     '';
                            'filepath'   'string'   []     '';
                            'check'      'string'   { 'on' 'off' }     'off';
                            'savemode'   'string'   { 'resave' 'onefile' 'twofiles' '' } '' });
if isstr(g), error(g); end;

% current filename without the .set
% ---------------------------------
[g.filepath filenamenoext ext] = fileparts( fullfile(g.filepath, g.filename) ); ext = '.set';
g.filename = [ filenamenoext ext ];

% performing extended syntax check
% --------------------------------
if strcmpi(g.check, 'on')
    fprintf('Pop_saveset: Performing extended dataset syntax check...\n');
    EEG = eeg_checkset(EEG, 'eventconsistency');
end

% default saving otion
% --------------------
save_as_dat_file = 0;
if strcmpi(g.savemode, 'resave')
    g.filename = EEG.filename;
    g.filepath = EEG.filepath;
    if isfield(EEG, 'datfile')
        if isempty(EEG.datfile)
            EEG = rmfield(EEG, 'datfile');
        else
            save_as_dat_file = 1;
        end;
    end;
else
    EEG.filename    = g.filename;
    EEG.filepath    = g.filepath;
    eeg_optionsbackup;
    eeg_options;
    if isempty(g.savemode)
        if option_savematlab, g.savemode = 'onefile';
        else                  g.savemode = 'twofiles';
        end;
    end;
    if strcmpi(g.savemode, 'twofiles')
        save_as_dat_file = 1;
        EEG.datfile = [ filenamenoext '.dat' ];
    end;
end;

% Saving data as float and Matlab
% -------------------------------
tmpica       = EEG.icaact;
EEG.icaact   = [];
tmpdata      = single(reshape(EEG.data, EEG.nbchan,  EEG.pnts*EEG.trials));
v = version;
%try, 
    fprintf('Saving dataset...\n');
    if save_as_dat_file
        EEG.data = EEG.datfile;
        floatwrite( tmpdata', fullfile(EEG.filepath, EEG.data), 'ieee-le');
    end;
    if v(1) > 6, save(fullfile(EEG.filepath, EEG.filename), '-v6', '-mat', 'EEG');
    else         save(fullfile(EEG.filepath, EEG.filename), '-mat', 'EEG');
    end;
    if save_as_dat_file
        EEG.data = tmpdata;
    end;
%catch,
%    error('Pop_saveset: save error, out of space or file permission problem');
%end;

% try to delete old .fdt or .dat files
% ------------------------------------
tmpfilename = fullfile(EEG.filepath, [ filenamenoext '.fdt' ]);
if exist(tmpfilename) == 2
    disp('Old .fdt file format detected on disk, now replaced by .dat file; trying to erase file...');
    try,
        delete(tmpfilename);
        disp('Delete sucessfull.');
    catch, disp('Error while attempting to remove file'); 
    end;
end;
if save_as_dat_file == 0
    tmpfilename = fullfile(EEG.filepath, [ filenamenoext '.dat' ]);
    if exist(tmpfilename) == 2
        disp('Old .dat file detected on disk for this dataset, deleting file to avoid confusion...');
        try,
            delete(tmpfilename);
            disp('Delete sucessfull.');
        catch, disp('Error while attempting to remove file'); 
        end;
    end;
end;

% recovering variables
% --------------------
EEG.icaact = tmpica;
if isnumeric(EEG.data) & v(1) < 7
    EEG.data   = double(reshape(tmpdata, EEG.nbchan,  EEG.pnts, EEG.trials));
end;

com = sprintf('%s = pop_saveset( %s, %s);', inputname(1), inputname(1), vararg2str(options));
return;

function num = popask( text )
	 ButtonName=questdlg2( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end;


