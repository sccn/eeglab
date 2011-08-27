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
%   'version' - ['6'|'7.3'] save Matlab file as version 6 (default) or
%               '7.3' (large files).
%
% Outputs:
%   EEG        - saved dataset (after extensive syntax checks)
%   ALLEEG     - saved datasets (after extensive syntax checks)
%
% Note: An individual dataset should be saved with a '.set' file extension
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: pop_loadset(), eeglab()
  
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

% 01-25-02 reformated help & license -ad 
 
function [EEG, com] = pop_saveset( EEG, varargin);

com = '';
if nargin < 1
	help pop_saveset;
	return;
end;
if isempty(EEG)  , error('Cannot save empty datasets'); end;

% empty filename (resave file)
emptyfilename = 0;
if nargin > 1
    if isempty(varargin{1}) & isempty(EEG.filename), emptyfilename = 1; end;
end;

if nargin < 2 | emptyfilename
    if length(EEG) >1, error('For reasons of consistency, this function  does not save multiple datasets any more'); end;
    % pop up window to ask for file type
    % ----------------------------------
    [filename, filepath] = uiputfile2('*.set', 'Save dataset with .set extension -- pop_saveset()'); 
    if ~isstr(filename), return; end;
    drawnow;
    options = { 'filename' filename 'filepath' filepath };
else
    % account for old calling format
    % ------------------------------
   if isempty(strmatch( lower(varargin{1}), { 'filename' 'filepath' 'savemode' 'check' }))
        options = { 'filename' varargin{1} };
        if nargin > 2
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
                            'version'    'string'   { '6','7.3' } '6';
                            'check'      'string'   { 'on','off' }     'off';
                            'savemode'   'string'   { 'resave','onefile','twofiles','' } '' });
if isstr(g), error(g); end;

% current filename without the .set
% ---------------------------------
if emptyfilename == 1, g.savemode = ''; end;
[g.filepath filenamenoext ext] = fileparts( fullfile(g.filepath, g.filename) ); ext = '.set';
g.filename = [ filenamenoext ext ];

% performing extended syntax check
% --------------------------------
if strcmpi(g.check, 'on')
    fprintf('Pop_saveset: Performing extended dataset syntax check...\n');
    EEG = eeg_checkset(EEG, 'eventconsistency');
    EEG = eeg_checkset(EEG);
else
    EEG = eeg_checkset(EEG);
end

% check for change in saving mode
% -------------------------------
eeglab_options;
if length(EEG) == 1
    if strcmpi(g.savemode, 'resave') & isfield(EEG, 'datfile') & option_savematlab
        disp('Note that your memory options for saving datasets does not correspond')
        disp('to the format of the datasets on disk (ignoring memory options)')
% $$$         but = questdlg2(strvcat('This dataset has an associated ''.dat'' file, but since you have', ...
% $$$                           'changed of saving mode, all the data will now be saved within the', ...
% $$$                           'Matlab file and the ''.dat'' file will be deleted.', ...
% $$$                           '(Note: Just press ''No'' if you do not know what you are doing)'), ...
% $$$                           'Warning: saving mode changed', 'Cancel', 'No, save as before', 'Yes, do it', 'Yes, do it');
% $$$         switch but
% $$$             case 'Cancel', return;
% $$$             case 'No, save as before', % nothing
% $$$             case 'Yes, do it', g.savemode = 'onefile';
% $$$         end;
% $$$         g.filename = EEG.filename;
% $$$         g.filepath = EEG.filepath;
    elseif strcmpi(g.savemode, 'resave') & ~isfield(EEG, 'datsave(file') & ~option_savematlab
        disp('Note that your memory options for saving datasets does not correspond')
        disp('to the format of the datasets on disk (ignoring memory options)')
% $$$         but = questdlg2(strvcat('This dataset does not have yet an associated ''.dat'' file, but since you have', ...
% $$$                           'changed of saving mode, all the data will now be saved within the ''.dat''', ...
% $$$                           'file and not in the Matlab file (as it is currently the case).', ...
% $$$                           '(Note: Just press ''No'' if you do not know what you are doing)'), ...
% $$$                           'Warning: saving mode changed', 'Cancel', 'No, save as before', 'Yes, do it', 'Yes, do it');
% $$$         switch but
% $$$             case 'Cancel', return;
% $$$             case 'No, save as before', % nothing
% $$$             case 'Yes, do it', g.savemode = 'twofiles';
% $$$         end;
% $$$         g.filename = EEG.filename;
% $$$         g.filepath = EEG.filepath;
    end;
end;

% default saving otion
% --------------------
save_as_dat_file = 0;
data_on_disk     = 0;
if strcmpi(g.savemode, 'resave')
    % process multiple datasets
    % -------------------------
    if length(EEG) > 1
        for index = 1:length(EEG)
            pop_saveset(EEG(index), 'savemode', 'resave');
            EEG(index).saved = 'yes';
        end;
        com = sprintf('%s = pop_saveset( %s, %s);', inputname(1), inputname(1), vararg2str(options));
        return;
    end;
    
    if strcmpi( EEG.saved, 'yes'), disp('Dataset has not been modified; No need to resave it.'); return; end;
    g.filename = EEG.filename;
    g.filepath = EEG.filepath;
    if isfield(EEG, 'datfile')
        if ~isempty(EEG.datfile)
            save_as_dat_file = 1;
        end;
    end;
    if isstr(EEG.data) & ~save_as_dat_file % data in .set file
        TMP = pop_loadset(EEG.filename, EEG.filepath);
        EEG.data = TMP.data;
        data_on_disk = 1;
    end;
else
    if length(EEG) >1, error('For reasons of consistency, this function  does not save multiple datasets any more'); end;
    if ~strcmpi(EEG.filename, g.filename) | ~strcmpi(EEG.filepath, g.filepath)
         EEG.datfile = '';
    end;
    EEG.filename    = g.filename;
    EEG.filepath    = g.filepath;
    if isempty(g.savemode)
        if option_savematlab, g.savemode = 'onefile';
        else                  g.savemode = 'twofiles';
        end;
    end;
    if strcmpi(g.savemode, 'twofiles')
        save_as_dat_file = 1;
        EEG.datfile = [ filenamenoext '.fdt' ];
    end;
end;

% Saving data as float and Matlab
% -------------------------------
tmpica       = EEG.icaact;
EEG.icaact   = [];
if ~isstr(EEG.data)
    if ~strcmpi(class(EEG.data), 'memmapdata') & ~strcmpi(class(EEG.data), 'single')
        tmpdata       = single(reshape(EEG.data, EEG.nbchan,  EEG.pnts*EEG.trials));
    else 
        tmpdata = EEG.data;
    end;
    no_resave_dat = 'no';
else 
    no_resave_dat = 'yes';
end;
v = version;
try, 
    fprintf('Saving dataset...\n');
    EEG.saved = 'yes';
    if save_as_dat_file
        if ~isstr(EEG.data)
            EEG.data = EEG.datfile;
            floatwrite( tmpdata, fullfile(EEG.filepath, EEG.data), 'ieee-le');
        end;
    else
        if isfield(EEG, 'datfile')
            if ~isempty(EEG.datfile)
                if exist(fullfile(EEG.filepath, EEG.datfile))
                    try, 
                        delete(fullfile(EEG.filepath, EEG.datfile));
                        disp('Deleting .dat/.fdt file on disk (all data is within the Matlab file)');
                    catch, end;
                end;
            end;
            EEG.datfile = [];
        end;
    end;
    
    if str2num(v(1)) > 6, 
        if strcmpi(g.version, '6') save(fullfile(EEG.filepath, EEG.filename), '-v6',   '-mat', 'EEG');
        else                       save(fullfile(EEG.filepath, EEG.filename), '-v7.3', '-mat', 'EEG');
        end;
    else                  save(fullfile(EEG.filepath, EEG.filename), '-mat', 'EEG');
    end;
    if save_as_dat_file & strcmpi( no_resave_dat, 'no' )
        EEG.data = tmpdata;
    end;
    
    % save ICA activities
    % -------------------
    icafile = fullfile(EEG.filepath, [EEG.filename(1:end-4) '.icafdt' ]);
    if isempty(EEG.icaweights) & exist(icafile)
        disp('ICA activation file found on disk, but no more ICA activities. Deleting file.');
        delete(icafile);
    end;
    if ~option_saveica & exist(icafile)
        disp('Options indicate not to save ICA activations. Deleting ICA activation file.');
        delete(icafile);
    end;
    if option_saveica & ~isempty(EEG.icaweights)
        if ~exist('tmpdata')
            TMP = eeg_checkset(EEG, 'loaddata');
            tmpdata = TMP.data;
        end;
        if isempty(tmpica)
             tmpica2 = (EEG.icaweights*EEG.icasphere)*tmpdata(EEG.icachansind,:);
        else tmpica2 = tmpica;
        end;
        tmpica2 = reshape(tmpica2, size(tmpica2,1), size(tmpica2,2)*size(tmpica2,3));
        floatwrite( tmpica2, icafile, 'ieee-le');
        clear tmpica2;
    end;
    
catch,
    rethrow(lasterror);
end;

% try to delete old .fdt or .dat files
% ------------------------------------
tmpfilename = fullfile(EEG.filepath, [ filenamenoext '.dat' ]);
if exist(tmpfilename) == 2
    disp('Deleting old .dat file format detected on disk (now replaced by .fdt file)');
    try,
        delete(tmpfilename);
        disp('Delete sucessfull.');
    catch, disp('Error while attempting to remove file'); 
    end;
end;
if save_as_dat_file == 0
    tmpfilename = fullfile(EEG.filepath, [ filenamenoext '.fdt' ]);
    if exist(tmpfilename) == 2
        disp('Old .fdt file detected on disk, deleting file the Matlab file contains all data...');
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
if data_on_disk
    EEG.data = 'in set file';
end;
if isnumeric(EEG.data) & v(1) < 7
    EEG.data   = double(reshape(tmpdata, EEG.nbchan,  EEG.pnts, EEG.trials));
end;
EEG.saved = 'justloaded';

com = sprintf('%s = pop_saveset( %s, %s);', inputname(1), inputname(1), vararg2str(options));
return;

function num = popask( text )
	 ButtonName=questdlg2( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end;


