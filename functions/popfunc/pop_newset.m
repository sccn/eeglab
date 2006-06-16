% pop_newset() - Edit/save EEG dataset structure information.
%
% Usage:
%   >> [ALLEEG EEG CURRENTSET] = pop_newset( ALLEEG, EEG, OLDSET, NEWSET,...
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

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.73  2006/05/10 14:06:09  arno
% better handling of condition when study is present
%
% Revision 1.72  2006/04/19 14:18:29  arno
% multiple dataset file
%
% Revision 1.71  2006/04/12 05:42:03  arno
% cancel button
%
% Revision 1.70  2006/04/10 21:04:35  arno
% remove dbug mssage
%
% Revision 1.69  2006/03/23 16:44:00  scott
% text and help
%
% Revision 1.68  2006/03/18 14:45:01  arno
% y
% cancel button
%
% Revision 1.67  2006/03/10 22:58:33  arno
% testing length of result
%
% Revision 1.66  2006/03/10 21:51:28  arno
% fixing retreive
%
% Revision 1.65  2006/03/10 21:40:35  arno
% nothing
%
% Revision 1.64  2006/03/10 21:38:13  arno
% use eeg_store
%
% Revision 1.63  2006/03/10 21:21:18  arno
% same
%
% Revision 1.62  2006/03/10 21:20:01  arno
% saving new datasets
%
% Revision 1.61  2006/03/10 20:58:24  arno
% fix shift
%
% Revision 1.60  2006/03/10 20:54:27  arno
% nothing
%
% Revision 1.59  2006/03/10 19:18:54  arno
% change text for edit box
%
% Revision 1.58  2006/03/02 23:34:31  arno
% nothing
%
% Revision 1.57  2006/02/17 20:57:24  arno
% update history, new dataset etc...
%
% Revision 1.56  2006/02/07 19:21:09  arno
% nothing
%
% Revision 1.55  2006/02/07 19:17:51  arno
% checking dataset consistency
%
% Revision 1.54  2006/02/04 00:02:06  arno
% fixing history
%
% Revision 1.53  2006/02/03 17:54:38  arno
% typo
%
% Revision 1.52  2006/02/03 00:47:30  arno
% same
%
% Revision 1.51  2006/02/03 00:46:49  arno
% fix study call
%
% Revision 1.50  2006/02/03 00:19:30  arno
% debug last changes
%
% Revision 1.49  2006/02/03 00:07:10  arno
% typo
%
% Revision 1.48  2006/02/03 00:05:45  arno
% dealing with study structure
%
% Revision 1.47  2006/02/02 22:43:08  arno
% processing empty dataset
%
% Revision 1.46  2006/02/02 00:40:17  arno
% removing data from file in ALLEEG
%
% Revision 1.45  2006/02/02 00:30:29  arno
% important changes for processing multiple datasets
%
% Revision 1.44  2006/02/02 00:11:39  arno
% nothing
%
% Revision 1.43  2006/02/01 18:39:29  arno
% automatically save datasets ...
%
% Revision 1.42  2006/02/01 07:11:07  arno
% typo
%
% Revision 1.41  2006/02/01 07:09:48  arno
% retreiving from multiple datasets
%
% Revision 1.40  2006/02/01 06:53:34  arno
% new saving data etc...
%
% Revision 1.39  2006/02/01 06:33:42  arno
% revert version 1.37
%
% Revision 1.37  2006/02/01 00:45:19  arno
% newset
%
% Revision 1.36  2006/01/31 22:36:30  arno
% allowing no options to be set
%
% Revision 1.35  2006/01/31 20:21:50  arno
% eeglab options
%
% Revision 1.34  2006/01/31 00:32:02  arno
% fix problem for saved field
%
% Revision 1.33  2006/01/31 00:23:48  arno
% allow no previous dataset
%
% Revision 1.32  2006/01/30 22:50:16  arno
% new format
%
% Revision 1.31  2005/11/02 18:35:28  arno
% nothing
%
% Revision 1.30  2005/09/27 22:03:38  arno
% save multiple datasets
%
% Revision 1.29  2005/08/16 17:54:39  scott
% edit help message. EEG.changes_not_saved -> EEG.saved   -sm
%
% Revision 1.28  2005/08/08 18:43:33  arno
% do not erase filename and filepath
%
% Revision 1.27  2005/08/08 18:40:59  arno
% fix saving file
%
% Revision 1.26  2005/08/04 17:25:18  arno
% typo
%
% Revision 1.25  2005/08/04 16:28:55  arno
% modified -> changes to save
%
% Revision 1.24  2005/08/04 16:27:30  arno
% implement EEG.modified
%
% Revision 1.23  2005/08/04 15:36:48  arno
% remove option of keeping only 1 dataset
%
% Revision 1.22  2005/07/30 01:22:24  arno
% allowing to remove channels for multiple datasets
%
% Revision 1.21  2004/07/07 19:07:30  arno
% return empty if cancel
%
% Revision 1.20  2003/12/05 20:01:05  arno
% eeg_hist for history
%
% Revision 1.19  2003/12/05 00:48:28  arno
% saving setname
%
% Revision 1.18  2003/07/20 19:37:11  scott
% typo
%
% Revision 1.17  2003/05/30 16:48:31  arno
% removing debug message
%
% Revision 1.16  2003/05/12 15:35:22  arno
% updated output command
%
% Revision 1.15  2003/02/27 00:38:58  arno
% typo
%
% Revision 1.14  2003/02/03 20:01:58  arno
% debugging overwrite option call from the command line
%
% Revision 1.13  2003/02/03 01:46:42  scott
% header edits -sm
%
% Revision 1.12  2002/12/04 19:11:35  arno
% macOSX directory compatibility
%
% Revision 1.11  2002/11/14 17:40:55  arno
% comment -> description
%
% Revision 1.10  2002/10/23 15:04:00  arno
% isppc -> computer
%
% Revision 1.9  2002/10/15 23:28:07  arno
% debugging dataset save
%
% Revision 1.8  2002/10/15 17:07:19  arno
% drawnow
%
% Revision 1.7  2002/10/10 17:08:29  arno
% output command ';'
%
% Revision 1.6  2002/08/14 01:29:17  arno
% updating save
%
% Revision 1.5  2002/08/12 18:33:48  arno
% questdlg2
%
% Revision 1.4  2002/05/03 03:05:11  arno
% editing interface
%
% Revision 1.3  2002/05/03 02:59:44  arno
% debugging cancel
%
% Revision 1.2  2002/04/26 02:51:26  arno
% adding com parameter
%
% Revision 1.1  2002/04/26 02:46:37  arno
% Initial revision
%

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

if nargin < 3
   help pop_newset;
   return;
end;
CURRENTSET = OLDSET;
com = sprintf('[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, %s); ', vararg2str( { OLDSET varargin{:} } ));

[g varargin] = finputcheck(varargin, { ...
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
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    end;
    return;
elseif ~isempty(g.retrieve) % command line call
    % ***************************************************
    % case 4 -> single dataset, does not have to be saved, 
    %           retrieving another dataset
    % ***************************************************
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
elseif length(varargin) == 0 & length(EEG) == 1 % if several arguments, assign values 
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
    % case 5 -> single dataset, has to be saved, several dataset to retrieve
    % ***************************************************
    if length(g.retrieve) > 1 | ( g.study & ~isempty(g.retrieve)) % selecting several datasets or a study is present
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
    % ***************************************************
    % case 7 -> single dataset modified, study is absent, old copy has to
    %           be flush to disk or overwritten 
    % ***************************************************
    elseif ~saved & option_storedisk
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
        text_old = 'Some changes have not been saved. What do you want to do with the old dataset?';
    else
    % ***************************************************
    % case 9 -> single dataset not modified, study is absent, no constraint on saving
    % ***************************************************
        text_old       = 'What do you want to do with the old dataset (not modified since last saved)?';
        cb_saveold     = '';
        cb_overwrite   = 'Overwrite current dataset|New dataset';
    end;
        
    geometry    = { [1] [0.12 0.5 1 0.5] [0.12 0.5 1 0.5] [1] [1] [0.12 1.8 0.1 0.1] [0.12 0.5 1 0.5] };
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
            uilist = { uilist{:} ...
                 { 'Style', 'checkbox'  , 'string', '', 'callback', cb_loadold 'tag' 'cb_loadold' } ...
                 { 'Style', 'text',       'string', 'Reload copy from disk (will be done after optional saving above)' } {} {} };
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
                        EEG.saved = 'yes';
                        EEG = pop_saveset(EEG, [ filename ext ], filepath);
	 case 'saveold',    [filepath filename ext] = fileparts( args{ind+1} );
                        ALLEEG(OLDSET).saved = 'yes';
                        EEG = pop_saveset(ALLEEG(OLDSET), [ filename ext ], filepath);
	 case 'overwrite' , if strcmpi(args{ind+1}, 'on') | strcmpi(args{ind+1}, 'yes')
                            overWflag = 1; % so it can be done at the end
                        end;
	 otherwise, error(['pop_newset error: unrecognized key ''' args{ind} '''']); 
    end;
end;

% remove data from file if necessary
% ----------------------------------
if option_storedisk & ~isempty(ALLEEG) & OLDSET ~= 0
    ALLEEG(OLDSET) = update_datafield(ALLEEG(OLDSET));
end;

% moving/erasing/creating datasets
% --------------------------------
if ~isempty(g.retrieve)
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
        [ALLEEG, EEG] = eeg_store( ALLEEG, EEG, OLDSET);
    else
        if strcmpi(EEG.saved, 'yes')
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0); % 0 means that it is saved on disk
        else
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG);
        end;
    end;
end;        

com = sprintf('[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, %s); ', vararg2str( { OLDSET args{:} } ));
return;

function num = popask( text )
	 ButtonName=questdlg2( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end;

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
    EEG.icaact = [];
