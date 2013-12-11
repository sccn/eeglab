% pop_newset() - Edit/save EEG dataset structure information.
%
% Usage:
%   >> [ALLEEG EEG CURRENTSET] = pop_newset( ALLEEG, EEG, CURRENTSET,...
%                                            'key', val,...);
% Inputs and outputs:
%   ALLEEG     - array of EEG dataset structures
%   EEG        - current dataset structure or structure array
%   CURRENTSET - index(s) of the current EEG dataset(s) in ALLEEG
%
% Optional inputs:
%   'setname'     - ['string'] name for the new dataset
%   'comments'    - ['string'] comments on the new dataset
%   'overwrite'   - ['on'|'off'] overwrite the old dataset
%   'saveold'     - ['filename'] filename in which to save the old dataset
%   'savenew'     - ['filename'] filename in which to save the new dataset
%   'retrieve'    - [index] retrieve the old dataset (ignore recent changes)
%
% Note: Calls eeg_store() which may modify the variable ALLEEG 
%       containing the current dataset(s).
%
% Author: Arnaud Delorme, CNL / Salk Institute, 23 Arpil 2002
%
% See also: eeg_store(), pop_editset(), eeglab()

% Copyright (C) 23 Arpil 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

%   'aboutparent' - ['on'|'off'] insert reference to parent dataset in the comments

% testing: both options must be tested with datasets that have 2
%          files (.set and .dat) and with dataset using a single file
%
% Option when only one dataset is kept in memory
% 
% *********** When STUDY present ***********
% Study selected -> select single dataset
% Study selected -> select multiple datasets
% Several dataset selected -> select single dataset
% Several dataset selected -> select other several datasets
% Several dataset selected -> select study
% Dataset (not modified) selected -> select study 
% Dataset (not modified) selected -> select multiple datasets 
% Dataset (not modified) selected -> select other single dataset 
% Dataset (not modified) selected -> create new dataset (eg. resample) 
% Dataset (modified) selected -> select study 
% Dataset (modified) selected -> select multiple datasets 
% Dataset (modified) selected -> select other single dataset 
% Dataset (modified) selected -> create new dataset (eg. resample) 
% *********** When study absent ************
% Several dataset selected -> select single dataset
% Several dataset selected -> select other several datasets
% Dataset (not modified) selected -> select multiple datasets 
% Dataset (not modified) selected -> select other single dataset 
% Dataset (not modified) selected -> create new dataset (eg. resample) 
% Dataset (modified) selected -> select multiple datasets 
% Dataset (modified) selected -> select other single dataset 
% Dataset (modified) selected -> create new dataset (eg. resample) 
%
% Option when several datasets can be kept in memory
% 
% *********** When STUDY present ***********
% Study selected -> select single dataset
% Study selected -> select multiple datasets
% Several dataset selected -> select single dataset
% Several dataset selected -> select other several datasets
% Several dataset selected -> select study
% Dataset (not modified) selected -> select study 
% Dataset (not modified) selected -> select multiple datasets 
% Dataset (not modified) selected -> select other single dataset 
% Dataset (not modified) selected -> create new dataset (eg. resample) 
% Dataset (modified) selected -> select study 
% Dataset (modified) selected -> select multiple datasets 
% Dataset (modified) selected -> select other single dataset 
% Dataset (modified) selected -> create new dataset (eg. resample) 
% *********** When study absent ************
% Several dataset selected -> select single dataset
% Several dataset selected -> select other several datasets
% Dataset (not modified) selected -> select multiple datasets 
% Dataset (not modified) selected -> select other single dataset 
% Dataset (not modified) selected -> create new dataset (eg. resample) 
% Dataset (modified) selected -> select multiple datasets 
% Dataset (modified) selected -> select other single dataset 
% Dataset (modified) selected -> create new dataset (eg. resample) 

function [ALLEEG, EEG, CURRENTSET, com] = pop_newset( ALLEEG, EEG, OLDSET, varargin);
% pop_newset( ALLEEG, EEG, 1, 'retrieve', [], 'study', [1] (retreiving a study)

verbose = 0;
if nargin < 3
   help pop_newset;
   return;
end;
CURRENTSET = OLDSET;
com = sprintf('[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, %s); ', vararg2str( { OLDSET varargin{:} } ));

[g varargin] = finputcheck(varargin, { ...
                    'gui'           'string'     { 'on';'off' }   'on'; % []=none; can be multiple numbers
                    'retrieve'      'integer'    []               []; % []=none; can be multiple numbers
                    'study'         'integer'    [0 1]            0;  % important because change behavior for modified datasets
                    }, 'pop_newset', 'ignore');
if isstr(g), error(g); end;
eeglab_options;
if length(EEG) > 1
    % ***************************************************
    % case 1 -> multiple datasets in memory (none has to be saved), retrieving single dataset
    % ***************************************************
    if ~isempty(g.retrieve)
        [EEG, ALLEEG, CURRENTSET] = eeg_retrieve( ALLEEG, g.retrieve);
    elseif length(OLDSET) == length(EEG)
    % ***************************************************
    % case 2 -> multiple datasets processed, storing them
    % ***************************************************
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, OLDSET);
    else
    % ***************************************************
    % case 3 -> multiple datasets processed, storing new copies (not used in EEGLAB)
    % ***************************************************
        if strcmpi(EEG(1).saved, 'justloaded')
            for ieeg = 1:length(EEG)
                [ALLEEG TMP OLDSET] = pop_newset(ALLEEG, EEG(ieeg), OLDSET);
            end;
            CURRENTSET = OLDSET;
            EEG = TMP;
        else
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        end;
    end;
    return;
elseif ~isempty(g.retrieve) % command line call
    % ***************************************************
    % case 4 -> single dataset, does not have to be saved, 
    %           retrieving another dataset
    % ***************************************************
    if verbose, disp('Case 4'); end;
    if ~(option_storedisk & strcmpi(EEG.saved, 'no'))
        if strcmpi(EEG.saved, 'yes') & option_storedisk
            fprintf('pop_newset(): Dataset %d has not been modified since last save, so did not resave it.\n', OLDSET);
            EEG = update_datafield(EEG);
        
            tmpsave = EEG.saved;
            EEG = eeg_hist(EEG, com);        
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, OLDSET);
            EEG.saved            = tmpsave; % eeg_store automatically set it to 'no'
            ALLEEG(OLDSET).saved = tmpsave;
        end;
        
        if ~isempty(g.retrieve)
            [EEG, ALLEEG, CURRENTSET] = eeg_retrieve( ALLEEG, g.retrieve);
        end;
        return;
    end;
end;

if isempty(EEG)
    args = { 'retrieve', OLDSET }; % cancel    
elseif length(varargin) == 0 & length(EEG) == 1 & strcmpi(g.gui, 'on') % if several arguments, assign values 
    % popup window parameters	
    % -----------------------
    text_new   = 'What do you want to do with the new dataset?';
    comcomment = ['tmpuserdat = get(gcbf, ''userdata'');' ...
				  'tmpuserdat = pop_comments(tmpuserdat, ''Edit dataset comments'');' ...
				  'set(gcbf, ''userdata'', tmpuserdat); clear tmpuserdat;'];
	comsavenew   = ['[tmpfile tmppath] = uiputfile(''*.set'', ''Enter filename''); drawnow;' ...
				  'if tmpfile ~= 0,' ...
				  '    set(findobj(''parent'', gcbf, ''tag'', ''filenamenew''), ''string'', fullfile(tmppath, tmpfile));' ...
				  'end;' ...
				  'clear tmpuserdat tmpfile tmppath;'];
    cb_savenew   = [ 'set(findobj(gcbf, ''userdata'', ''filenamenew''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));' ];
    cb_saveold   = [ 'set(findobj(gcbf, ''userdata'', ''filenameold''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));' ];
	comsaveold   = ['[tmpfile tmppath] = uiputfile(''*.set'', ''Enter filename''); drawnow;' ...
				  'if tmpfile ~= 0,' ...
				  '    set(findobj(''parent'', gcbf, ''tag'', ''filenameold''), ''string'', fullfile(tmppath, tmpfile));' ...
				  'end;' ...
				  'clear tmpuserdat tmpfile tmppath;'];
    enable_saveold = 'off';
    enable_savenew = 'off';
    value_saveold  = 0;
    value_savenew  = 0;
    value_owrt   = 0;
    cb_owrt      = '';
	userdat = EEG.comments;
    
    % status of parent dataset etc...
    saved = 1;
    filenameold = '';
    filenamenew = '';
    if ~isempty(ALLEEG) & any(OLDSET ~= 0) & length(OLDSET) == 1
        if strcmpi(ALLEEG(OLDSET).saved, 'no')
            saved = 0;
            if ~isempty(ALLEEG(OLDSET).filename)
                filenameold = fullfile(ALLEEG(OLDSET).filepath, ALLEEG(OLDSET).filename);
            end;
        end;
    end;
    overwrite_or_save = 0;
    have_to_save_new  = 0;

    % ***************************************************
    % case 5 -> single dataset, has to be saved, one dataset to retreive and study present or several dataset to retrieve
    % ***************************************************
    if length(g.retrieve) > 1 | ( g.study & ~isempty(g.retrieve)) % selecting several datasets or a study is present
        if verbose, disp('Case 5'); end;
        text_new = 'Current dataset has not been saved. Saved it or reload it from disk.';
        text_old = '';
        have_to_save_new = 1;
        value_savenew    = 1;
        enable_savenew = 'on';
        filenamenew    = fullfile(EEG.filepath, EEG.filename);
        cb_savenew     = [ 'if ~get(gcbo, ''value'') & ~get(findobj(gcbf, ''tag'', ''cb_loadold''), ''value''),' ...
                       '   set(gcbo, ''value'', 1);' ...
                       '   warndlg2(strvcat(''You must enter a filename for the dataset or use the copy on disk!'','' '',' ...
                       '   ''It must be saved or you must use the old dataset on disk because, by your current memory'',' ...
                       '   ''option, only one full dataset or study can be kept in memory. Note that if you choose to'',' ...
                       '   ''save the dataset, this will be taken into account in the study.''));' ...
                       'else, ' cb_savenew ...
                       'end;' ];
    % ***************************************************
    % case 6 -> single dataset modified or not, study is present (the old
    %           dataset has to be replaced to preserve study consistency)
    % ***************************************************
    elseif g.study == 1 & isempty(g.retrieve)
        if verbose, disp('Case 6'); end;
        if saved
            text_old       = 'The old dataset has not been modified since last saved. What do you want to do with it?';
            cb_saveold     = '';
        else
            text_old = 'Some changes have not been saved. What do you want to do with the old dataset?';
        end;
        cb_owrt      = [ ...
                       '   set(gcbo, ''value'', 1);' ...
                       '    warndlg2(strvcat(''Cannot unset the overwrite checkbox!'','' '',' ...
                           '''The old dataset must be overwriten since all datasets'',' ...
                           '''must be in the STUDY.''), ''warning'');' ];
        value_owrt   = 1;
        filenamenew    = fullfile(EEG.filepath, EEG.filename);
    % ***************************************************
    % case 7 -> single dataset modified, study is absent, old copy has to
    %           be flush to disk or overwritten 
    % ***************************************************
    elseif ~saved & option_storedisk
        if verbose, disp('Case 7'); end;
        text_old = 'Some changes have not been saved. What do you want to do with the old dataset?';
        cb_saveold     = [ 'if ~get(findobj(gcbf, ''tag'', ''cb_owrt''), ''value''),' ...
                       '   set(gcbo, ''value'', 1);' ...
                       '    warndlg2(strvcat(''Cannot unset the save checkbox!'','' '',' ...
                           '''By your selected memory option, only one full dataset'',' ...
                           '''can be kept in memory. Thus, you must either'',' ...
                           '''save or delete/overwrite the old dataset.''));' ... 
                       'else, ' cb_saveold ...
                       'end;' ];
        cb_owrt      = [ 'if ~get(findobj(gcbf, ''tag'', ''cb_saveold''), ''value''),' ...
                       '   set(gcbo, ''value'', 1);' ...
                       '    warndlg2(strvcat(''Cannot unset the overwrite checkbox!'','' '',' ...
                           '''By your memory option, only one full dataset'',' ...
                           '''can be kept in memory. Thus, you must either'',' ...
                           '''save or delete/overwrite the old dataset.''));' ... 
                       'end;' ];
        enable_saveold = 'on';
        value_saveold  = 1;
        overwrite_or_save = 1;
    elseif ~saved
    % ***************************************************
    % case 8 -> single dataset modified, study is absent, no constraint on saving
    % ***************************************************
        if verbose, disp('Case 8'); end;
        text_old = 'Some changes have not been saved. What do you want to do with the old dataset?';
    else
    % ***************************************************
    % case 9 -> single dataset not modified, study is absent, no constraint on saving
    % ***************************************************
        if verbose, disp('Case 9'); end;
        text_old       = 'What do you want to do with the old dataset (not modified since last saved)?';
        cb_saveold     = '';
        cb_overwrite   = 'Overwrite current dataset|New dataset';
    end;
        
    geometry    = { [1] [0.15 0.5 1 0.5] [0.15 0.5 1 0.5] [1] [1] [0.15 1.8 0.1 0.1] [0.15 0.5 1 0.5] };
    geomvert    = [ ];
    uilist = { ...
         { 'style', 'text',       'string', text_new, 'fontweight', 'bold' } ...
         {} ...
         { 'Style', 'text',       'string', 'Name it:' } ...
		 { 'Style', 'edit',       'string', EEG.setname 'tag' 'namenew' } ...
		 { 'Style', 'pushbutton', 'string', 'Edit description', 'callback', comcomment } ...
         { 'Style', 'checkbox'  , 'string', '', 'callback', cb_savenew 'value' value_savenew 'tag' 'cb_savenew' } ...
         { 'Style', 'text',       'string', 'Save it as file:' } ...
         { 'Style', 'edit',       'string', filenamenew, 'tag', 'filenamenew' 'userdata' 'filenamenew'  'enable' enable_savenew } ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', comsavenew 'userdata' 'filenamenew' 'enable' enable_savenew } ...
         { } ...
         { 'style', 'text',       'string', text_old 'fontweight' 'bold' } ...
         { 'Style', 'checkbox'  , 'string', '' 'tag' 'cb_owrt' 'callback' cb_owrt 'value' value_owrt } ...
         { 'Style', 'text'      , 'string', 'Overwrite it in memory (set=yes; unset=create a new dataset)' } {} ...
         { } ...
         { 'Style', 'checkbox'  , 'string', '', 'callback', cb_saveold 'value' value_saveold 'tag' 'cb_saveold' } ...
         { 'Style', 'text'      , 'string', 'Save it as file:' } ...
         { 'Style', 'edit'      , 'string', filenameold, 'tag', 'filenameold'   'userdata' 'filenameold' 'enable' enable_saveold } ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', comsaveold 'userdata' 'filenameold' 'enable' enable_saveold } };
    
    % remove old dataset if not present
    % ---------------------------------
    if OLDSET == 0
        uilist = uilist(1:9);
        geometry = geometry(1:3);
    end;
    if isempty(cb_saveold)
        uilist(end-3:end) = [];
        geometry(end)     = [];
    end;        
    % update GUI for selecting multiple datasets
    % ------------------------------------------
    if length(g.retrieve) > 1 | ( g.study & ~isempty(g.retrieve)) % selecting several datasets or a study is present
        uilist = uilist(1:9);
        geometry = geometry(1:3);
        if ~isempty(EEG.filename)
            cb_loadold = [ 'if ~get(gcbo, ''value'') & ~get(findobj(gcbf, ''tag'', ''cb_savenew''), ''value''),' ...
                           '   set(gcbo, ''value'', 1);' ...
                           '   warndlg2(strvcat(''You must enter a filename for the dataset or use the copy on disk!'','' '',' ...
                           '   ''It must be saved or you must use the old dataset on disk because, by your current memory option, only one full'',' ...
                           '   ''dataset or study can be kept in memory. This will also affect the'',' ...
                           '   ''dataset at this position in the study.''));' ...
                           'end;' ];
            uilist{end+1} = { 'Style', 'checkbox'  , 'string', '', 'callback', cb_loadold 'tag' 'cb_loadold' };
            uilist{end+1} = { 'Style', 'text',       'string', 'Reload copy from disk (will be done after optional saving above)' };
            uilist{end+1} = {};
            uilist{end+1} = {};
             geometry = { geometry{:} [0.12 1.6 0.2 0.2] };
        end;
    end;
            
    % remove new dataset if already saved
    % -----------------------------------
    if ~isfield(EEG, 'saved')
        EEG = eeg_checkset(EEG);
    end;
    if strcmpi(EEG.saved, 'justloaded') | ~isempty(g.retrieve)
        if overwrite_or_save % only pop-up a window if some action has to be taken
            uilist = uilist(11:end);
            geometry = geometry(5:end);
            uilist(3) = {{ 'Style', 'text'      , 'string', 'Delete it from memory (set=yes)' }};
        elseif isempty(g.retrieve) % just loaded from disk
            
            % remove data from file for old dataset
            % -------------------------------------
            if option_storedisk & ~isempty(ALLEEG) & OLDSET ~= 0
                if ~isfield(ALLEEG(OLDSET), 'datfile'), ALLEEG(OLDSET).datfile = ''; end;
                ALLEEG(OLDSET) = update_datafield(ALLEEG(OLDSET));
            end;
            
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0); % 0 means that it is saved on disk
            com = '[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );';
            return;
        end;
    end;
    
    % show GUI (do not return if old dataset has to be saved or overwritten)
    % ----------------------------------------------------------------------
    cont = 1;
    while cont
        [result userdat tmp tags] = inputgui( 'geometry', geometry, 'uilist', uilist, 'helpcom', 'pophelp(''pop_newset'');', ...
                                     'title', 'Dataset info -- pop_newset()', ...
                                     'userdata', userdat, 'geomvert', geomvert);
        try, tags.cb_owrt;      catch, tags.cb_owrt = 0;    end;
        try, tags.cb_savenew;   catch, tags.cb_savenew = 0; end;
        try, tags.cb_saveold;   catch, tags.cb_saveold = 0; end;
        try, tags.cb_loadold;   catch, tags.cb_loadold = 0; end;
        try, tags.namenew;      catch, tags.namenew = EEG.setname; end;
        try, tags.filenamenew;  catch, tags.filenamenew = ''; end;
        cont = 0;
        if ~isempty(result)
            if overwrite_or_save & tags.cb_saveold % save but not overwrite
                if isempty(tags.filenameold)
                    warndlg2(strvcat('Error: You must enter a filename for the old dataset!',' ', ...
                           'The old dataset must be saved because, by your', ...
                           'current memory option, only one full dataset', ...
                           'can be kept in memory. Thus, you must either', ...
                           'save or overwrite the old dataset.'));
                    cont = 1;
                end;
            end;
        end;
        if have_to_save_new
            if isempty(result) % cancel
                com = '';
                drawnow;
                return;
            else
                if isempty(tags.filenamenew) & tags.cb_savenew
                    warndlg2(strvcat('Error: You must enter a filename for the dataset!',' ', ...
                           'It must be saved because, by your current memory option, only one full', ...
                           'dataset or study can be kept in memory. Note that all changes will be', ...
                           'taken into account when processing the STUDY.'));
                    cont = 1;
                end;
            end;
        end;
    end;
    drawnow;

    % decode parameters
    % -----------------
    args = {};
    if length(result) == 0,
        if isempty(g.retrieve)
            if isempty(OLDSET), error('Cancel operation'); end;
            args = { 'retrieve', OLDSET }; % cancel
        else
            com = '';
            return;
        end;
	else 
        % new dataset
        % -----------
        if ~strcmp(EEG.setname, tags.namenew )       
            args = { 'setname', tags.namenew };
        end;
        if tags.cb_savenew 
            if ~isempty(tags.filenamenew) 
                args = { args{:} 'savenew', tags.filenamenew };
            else
                disp('Warning: no filename given for new dataset, so it will not be saved to disk.');
            end;
        end;
        if ~strcmp(EEG.comments, userdat)
            args = { args{:} 'comments', userdat };
        end;
        if tags.cb_loadold    
            args = { args{:} 'reload' 'on' };
        end;
        
        % old dataset
        % -----------
        if tags.cb_owrt
            args = { args{:} 'overwrite' 'on' };
        end;
        if tags.cb_saveold 
            if ~isempty(tags.filenameold) 
                args = { args{:} 'saveold', tags.filenameold };
            else
                disp('Warning: no file name given for the old dataset, so it will not be saved to disk.');
            end;
        end;
    end;
elseif length(EEG) > 1
    % processing multiple datasets
    % ----------------------------
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, OLDSET ); % it is possible to undo the operation here
    com = '[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET );';
    return;
else
    % no interactive inputs
    args = varargin;
end;

% assigning values
% ----------------
overWflag    = 0;
if isempty(g.retrieve) & ~isempty(EEG)
    if strcmpi(EEG.saved, 'justloaded')
        EEG.saved = 'yes';
    else
        EEG.saved = 'no';
    end;
end;

for ind = 1:2:length(args)
    switch lower(args{ind})
	 case 'setname'   , EEG.setname = args{ind+1}; EEG = eeg_hist(EEG, [ 'EEG.setname=''' EEG.setname ''';' ]); 
	 case 'comments'  , EEG.comments = args{ind+1};
	 case 'reload'    , EEG = pop_loadset('filename', EEG.filename, 'filepath', EEG.filepath, 'loadmode', 'info');
                        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, OLDSET);
                        ALLEEG(OLDSET).saved = 'yes';
	 case 'retrieve'  , if ~isempty(ALLEEG) & args{ind+1} ~= 0
                            EEG = eeg_retrieve(ALLEEG, args{ind+1}); 
                        else
                            EEG = eeg_emptyset;
                        end;
                        com = ''; return;
	 case { 'save' 'savenew' }, [filepath filename ext] = fileparts( args{ind+1} );
                        EEG = pop_saveset(EEG, [ filename ext ], filepath);
	 case 'saveold',    [filepath filename ext] = fileparts( args{ind+1} );
                        TMPEEG = pop_saveset(ALLEEG(OLDSET), [ filename ext ], filepath);
                        [ALLEEG] = eeg_store(ALLEEG, TMPEEG, OLDSET);
                        ALLEEG(OLDSET).saved = 'yes';
	 case 'overwrite' , if strcmpi(args{ind+1}, 'on') | strcmpi(args{ind+1}, 'yes')
                            overWflag = 1; % so it can be done at the end
                        end;
	 otherwise, error(['pop_newset error: unrecognized key ''' args{ind} '''']); 
    end;
end;

% remove data from file if necessary
% ----------------------------------
if option_storedisk & ~isempty(ALLEEG) & OLDSET ~= 0
    if ~isfield(ALLEEG, 'datfile'), ALLEEG(OLDSET).datfile = ''; end;
    ALLEEG(OLDSET) = update_datafield(ALLEEG(OLDSET));
end;

% moving/erasing/creating datasets
% --------------------------------
if ~isempty(g.retrieve)
    
    % in case the old dataset was modified
    % ------------------------------------
    if strcmpi(EEG.saved, 'yes')
        [ALLEEG, EEG] = eeg_store( ALLEEG, EEG, OLDSET);
        ALLEEG(OLDSET).saved = 'yes';
        EEG.saved            = 'yes';
    else
        [ALLEEG, EEG] = eeg_store( ALLEEG, EEG, OLDSET);
    end;

    % dataset retrieval
    % -----------------
    if overWflag % delete old dataset
        ALLEEG = pop_delset( ALLEEG, OLDSET);
    end;
    [EEG, ALLEEG, CURRENTSET] = eeg_retrieve( ALLEEG, g.retrieve);
else
    % new dataset
    % -----------
    if overWflag
        if strcmpi(EEG.saved, 'yes')
            [ALLEEG, EEG] = eeg_store( ALLEEG, EEG, OLDSET);
            ALLEEG(OLDSET).saved = 'yes';
            EEG.saved            = 'yes';
        else
            [ALLEEG, EEG] = eeg_store( ALLEEG, EEG, OLDSET);
        end;
    else
        if strcmpi(EEG.saved, 'yes') || strcmpi(EEG.saved, 'justloaded')
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0); % 0 means that it is saved on disk
        else
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG);
        end;
    end;
end;        

com = sprintf('[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, %s); ', vararg2str( { OLDSET args{:} 'gui' 'off' } ));
return;

function num = popask( text )
	 ButtonName=questdlg2( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end;

function EEG = update_datafield(EEG);
    if ~isfield(EEG, 'datfile'), EEG.datfile = ''; end;
    if ~isempty(EEG.datfile)
        EEG.data = EEG.datfile;
    else 
        EEG.data = 'in set file';
    end;
    EEG.icaact = [];
