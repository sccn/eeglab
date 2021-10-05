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
%                and the transposed data in a binary float '.fdt' file.
%                By default the option from the eeg_options.m file is 
%                used.
%   'version' - ['6'|'7'|'7.3'] save Matlab file as version 6, 7, or
%               '7.3' (default; as defined in eeg_option file).
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
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

% 01-25-02 reformated help & license -ad 
 
function [EEG, com] = pop_saveset( EEG, varargin);

com = '';
if nargin < 1
	help pop_saveset;
	return;
end
if isempty(EEG)  , error('Cannot save empty datasets'); end

% empty filename (resave file)
emptyfilename = 0;
if nargin > 1
    if isempty(varargin{1}) && isempty(EEG.filename), emptyfilename = 1; end
    if strcmpi(varargin{1},'savemode') 
        if length(EEG) == 1
            if isempty(EEG(1).filename), varargin{2} = ''; emptyfilename = 1; end
        else
            if any(cellfun(@isempty, { EEG.filename }))
                error('Cannot resave files who have not been saved previously');
            end
        end
    end
end

if nargin < 2 || emptyfilename
    if length(EEG) >1, error('For reasons of consistency, this function  does not save multiple datasets any more'); end
    % pop up window to ask for file
    [filename, filepath] = uiputfile2('*.set', 'Save dataset with .set extension -- pop_saveset()'); 
    if ~ischar(filename), return; end
    drawnow;
    options = { 'filename' filename 'filepath' filepath };
else
    % account for old calling format
    % ------------------------------
   if isempty(strmatch( lower(varargin{1}), { 'filename' 'filepath' 'savemode' 'check' }))
        options = { 'filename' varargin{1} };
        if nargin > 2
            options = { options{:} 'filepath' varargin{2} };
        end
    else
        options = varargin;
    end
end

% decode input parameters
% -----------------------
eeglab_options;
defaultSave = fastif(option_saveversion6, '6', '7.3');
g = finputcheck(options,  { 'filename'   'string'   []     '';
                            'filepath'   'string'   []     '';
                            'version'    'string'   { '6','7','7.3' } defaultSave;
                            'check'      'string'   { 'on','off' }     'off';
                            'savemode'   'string'   { 'resave','onefile','twofiles','' } '' });
if ischar(g), error(g); end

% current filename without the .set
% ---------------------------------
if emptyfilename == 1, g.savemode = ''; end
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
if length(EEG) == 1
    if strcmpi(g.savemode, 'resave') && isfield(EEG, 'datfile') && ~isempty(EEG.datfile) && ~option_savetwofiles
        disp('Note that your memory options for saving datasets does not correspond')
        disp('to the format of the datasets on disk (ignoring memory options)')
    elseif strcmpi(g.savemode, 'resave') && ~isfield(EEG, 'datfile') && option_savetwofiles
        disp('Note that your memory options for saving datasets does not correspond')
        disp('to the format of the datasets on disk (ignoring memory options)')
    end
end

% default saving option
% ---------------------
save_as_dat_file = 0;
data_on_disk     = 0;
if strcmpi(g.savemode, 'resave')
    % process multiple datasets
    % -------------------------
    
    % check if called from EEGLAB
    % if this is not the case, resave the file anyway
    calledFromEEGLABFlag = true;
    try
        db = dbstack;
        eeglabp = fileparts(which('eeglab.m'));
        if length(db) > 1
            if isempty(strfind(which(db(2).file), eeglabp))
                calledFromEEGLABFlag = false;
            end
        end
    catch
    end
    
    if length(EEG) > 1
        for index = 1:length(EEG)
            if strcmpi( EEG(index).saved, 'yes') && calledFromEEGLABFlag
                disp('Dataset has not been modified; No need to resave it.'); 
            else
                pop_saveset(EEG(index), 'savemode', 'resave');
                EEG(index).saved = 'yes';
            end
        end
        if nargout > 1
            com = sprintf('EEG = pop_saveset( EEG, %s);', vararg2str(options));
        end
        return;
    end
    
    if strcmpi( EEG.saved, 'yes') && calledFromEEGLABFlag
        disp('Dataset has not been modified; No need to resave it.'); 
        return; 
    end
    g.filename = EEG.filename;
    g.filepath = EEG.filepath;
    if isfield(EEG, 'datfile')
        if ~isempty(EEG.datfile)
            save_as_dat_file = 1;
        end
    end
    if ischar(EEG.data) && ~save_as_dat_file % data in .set file
        TMP = pop_loadset(EEG.filename, EEG.filepath);
        EEG.data = TMP.data;
        data_on_disk = 1;
    end
else
    if length(EEG) >1, error('For reasons of consistency, this function  does not save multiple datasets any more'); end
    if ~strcmpi(EEG.filename, g.filename) || ~strcmpi(EEG.filepath, g.filepath)
         EEG.datfile = '';
    end
    EEG.filename    = g.filename;
    EEG.filepath    = g.filepath;
    if isempty(g.savemode)
        if option_savematlab, g.savemode = 'onefile';
        else                  g.savemode = 'twofiles';
        end
    end
    if strcmpi(g.savemode, 'twofiles')
        save_as_dat_file = 1;
        EEG.datfile = [ filenamenoext '.fdt' ];
        option_savetwofiles = true; % BUG: Missing, thus it will save data in .set file as well as fdt
    end
end

% Saving data as float and Matlab
% -------------------------------
tmpica       = EEG.icaact;
EEG.icaact   = [];
if ~ischar(EEG.data)
    if ~strcmpi(class(EEG.data), 'memmapdata') && ~strcmpi(class(EEG.data), 'mmo') && ~strcmpi(class(EEG.data), 'single')
        tmpdata       = single(reshape(EEG.data, EEG.nbchan,  EEG.pnts*EEG.trials));
    else 
        tmpdata = EEG.data;
    end
    no_resave_dat = 'no';
else 
    no_resave_dat = 'yes';
end
v = version;
try, 
    fprintf('Saving dataset...\n');
    EEG.saved = 'yes';
    if save_as_dat_file
        if ~ischar(EEG.data)
            EEG.data = EEG.datfile;
            tmpdata = floatwrite( tmpdata, fullfile(EEG.filepath, EEG.data), 'ieee-le');
        end
    end
    
    % reload daata before saving if necessary
    tmpdata2 = [];
    if ~option_savetwofiles && ischar(EEG.data)
        tmpdata2 = EEG.data;
        EEG.data = eeg_getdatact(EEG);
        EEG.datfile = [];
    end
    
    try
        if option_saveasstruct
            if strcmpi(g.version, '6') save(fullfile(EEG.filepath, EEG.filename), '-v6',   '-mat', '-struct', 'EEG');
            elseif strcmpi(g.version, '7') save(fullfile(EEG.filepath, EEG.filename), '-v7', '-mat', '-struct', 'EEG');
            else                       save(fullfile(EEG.filepath, EEG.filename), '-v7.3', '-mat', '-struct', 'EEG');
            end
        else
            if strcmpi(g.version, '6') save(fullfile(EEG.filepath, EEG.filename), '-v6',   '-mat', 'EEG');
            elseif strcmpi(g.version, '7') save(fullfile(EEG.filepath, EEG.filename), '-v7', '-mat', 'EEG');
            else                       save(fullfile(EEG.filepath, EEG.filename), '-v7.3', '-mat', 'EEG');
            end
        end
    catch
        disp('Cannot save in v6 format (Unicode characters?), using default (possible x20 slow down)');
        save(fullfile(EEG.filepath, EEG.filename), '-mat', 'EEG');
    end
    
    % put back the EEG structure field as they were
    if ~isempty(tmpdata2)
        EEG.data = tmpdata2;
    end
    
    if save_as_dat_file && strcmpi( no_resave_dat, 'no' )
        EEG.data = tmpdata;
    end
        
catch
    rethrow(lasterror);
end

% try to delete old .fdt or .dat files
% ------------------------------------
tmpfilename = fullfile(EEG.filepath, [ filenamenoext '.dat' ]);
if exist(tmpfilename) == 2
    disp('Deleting old .dat file format detected on disk (now replaced by .fdt file)');
    try,
        delete(tmpfilename);
        disp('Delete successful.');
        EEG.datfile = [];
    catch, disp('Error while attempting to remove file'); 
    end
end
if save_as_dat_file == 0
    tmpfilename = fullfile(EEG.filepath, [ filenamenoext '.fdt' ]);
    if exist(tmpfilename) == 2
        disp('Old .fdt file detected on disk, deleting file since the Matlab file now contains all the data');
        try
            delete(tmpfilename);
            disp('Delete successful.');
            EEG.datfile = [];
        catch, disp('Error while attempting to remove file'); 
        end
    end
end

% recovering variables
% --------------------
EEG.icaact = tmpica;
if data_on_disk
    EEG.data = 'in set file';
end
if isnumeric(EEG.data) && v(1) < 7
    EEG.data   = double(reshape(tmpdata, EEG.nbchan,  EEG.pnts, EEG.trials));
end
EEG.saved = 'justloaded';
if nargout > 1
    com = sprintf('EEG = pop_saveset( EEG, %s);', vararg2str(options));
end
return;

function num = popask( text )
	 ButtonName=questdlg2( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end


