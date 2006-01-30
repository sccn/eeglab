% pop_newset() - Edit/save EEG dataset structure information.
%
% Usage:
%   >> [ALLEEG EEG CURRENTSET] = pop_newset( ALLEEG, EEG, CURRENTSET, ...
%                                            'key', val,...);
% Inputs and outputs:
%   ALLEEG     - array of EEG dataset structures
%   EEG        - current dataset structure or structure array
%   CURRENTSET - index/indices of the current EEG dataset(s) in ALLEEG
%
% Optional inputs:
%   'setname'     - Name of the new dataset
%   'comments'    - ['string'] comments on the new dataset
%   'overwrite'   - ['on'|'off'] overwrite the old dataset
%   'savenew'     - ['filename'] filename to use to save the new dataset
%   'saveold'     - ['filename'] filename to use to save the old dataset
%   'retrieve'    - [index] retrieve the old dataset (ignore changes)
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

function [ALLEEG, EEG, CURRENTSET, com] = pop_newset( ALLEEG, EEG, CURRENTSET, varargin);

com = '';
if nargin < 3
   help pop_newset;
   return;
end;

if nargin < 4 & length(EEG) == 1 % if several arguments, assign values 
    % popup window parameters	
    % -----------------------
    comcomment = ['tmpuserdat = get(gcbf, ''userdata'');' ...
				  'tmpuserdat = pop_comments(tmpuserdat, ''Edit dataset comments'');' ...
				  'set(gcbf, ''userdata'', tmpuserdat); clear tmpuserdat;'];
	comsave1   = ['[tmpfile tmppath] = uiputfile(''*.set'', ''Enter filename''); drawnow;' ...
				  'if tmpfile ~= 0,' ...
				  '    set(findobj(''parent'', gcbf, ''tag'', ''saveedit1''), ''string'', fullfile(tmppath, tmpfile));' ...
				  'end;' ...
				  'clear tmpuserdat tmpfile tmppath;'];
    cb_save1   = [ 'set(findobj(gcbf, ''userdata'', ''saveedit1''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));' ];
    cb_save2   = [ 'set(findobj(gcbf, ''userdata'', ''saveedit2''), ''enable'', fastif(get(gcbo, ''value''), ''on'', ''off''));' ];
	comsave2   = ['[tmpfile tmppath] = uiputfile(''*.set'', ''Enter filename''); drawnow;' ...
				  'if tmpfile ~= 0,' ...
				  '    set(findobj(''parent'', gcbf, ''tag'', ''saveedit2''), ''string'', fullfile(tmppath, tmpfile));' ...
				  'end;' ...
				  'clear tmpuserdat tmpfile tmppath;'];
    enable_save2 = 'off';
    value_save2  = 0;
    cb_owrt      = '';
	userdat = EEG.comments;
    
    % status of parent dataset etc...
    eeg_optionsbackup;
    eeg_options;
    saved = 1;
    filename = '';
    if ~isempty(ALLEEG)
        if strcmpi(ALLEEG(CURRENTSET).saved, 'no')
            saved = 0;
            if ~isempty(ALLEEG(CURRENTSET).filename)
                filename = fullfile(ALLEEG(CURRENTSET).filepath, ALLEEG(CURRENTSET).filename);
            end;
        end;
    end;
    overwrite_or_save = 0;
    if ~saved & option_storedisk
        text_old = 'What do you want to do with the old dataset (some changes have not been saved)?';
        cb_save2     = [ 'if ~get(findobj(gcbf, ''tag'', ''cb_owrt''), ''value''),' ...
                       '   set(gcbo, ''value'', 1);' ...
                       '    warndlg2(strvcat(''Cannot unset the save checkbox!'','' '',' ...
                           '''The old dataset has to be saved because,'',' ...
                           '''according to your memory option, only one full '',' ...
                           '''dataset can be retained in memory. Thus, you either'',' ...
                           ''' have to save or overwrite the old dataset''));' ... 
                       'else, ' cb_save2 ...
                       'end;' ];
        cb_owrt      = [ 'if ~get(findobj(gcbf, ''tag'', ''cb_save2''), ''value''),' ...
                       '   set(gcbo, ''value'', 1);' ...
                       '    warndlg2(strvcat(''Cannot unset the overwrite checkbox!'','' '',' ...
                           '''The old dataset has to be saved because,'',' ...
                           '''according to your memory option, only one full '',' ...
                           '''dataset can be retained in memory. Thus, you either'',' ...
                           ''' have to save or overwrite the old dataset''));' ... 
                       'end;' ];
        enable_save2 = 'on';
        value_save2  = 1;
        overwrite_or_save = 1;
    elseif ~saved
        text_old = 'What do you want to do with the old dataset (some changes have not been saved)?';
    else
        text_old       = 'What do you want to do with the old dataset (NOT modified since last change)?';
        cb_save2       = 'set(gcbo, ''value'', 0);';
        cb_overwrite   = 'Overwrite current dataset|New dataset';
    end;
        
    geometry    = { [1] [0.12 0.5 1 0.5] [0.12 0.5 1 0.5] [1] [1] [0.12 1.8 0.1 0.1] [0.12 0.5 1 0.5] };
    geomvert    = [ ];
    uilist = { { 'style', 'text',       'string', 'What do you want to do with the new dataset?', 'fontweight', 'bold' } ...
         {} { 'Style', 'text',       'string', 'Name it:' } ...
		 { 'Style', 'edit',       'string', EEG.setname } ...
		 { 'Style', 'pushbutton', 'string', 'Edit description', 'tooltipstring', 'Modify comments of this new dataset', 'callback', comcomment } ...
         { 'Style', 'checkbox'  , 'string', '', 'callback', cb_save1 } ...
         { 'Style', 'text',       'string', 'Save it as file:' } ...
         { 'Style', 'edit',       'string', '', 'tag', 'saveedit1' 'userdata' 'saveedit1'  'enable' 'off' } ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', comsave1 'userdata' 'saveedit1' 'enable' 'off' } ...
         { } ...
         { 'style', 'text', 'string', text_old 'fontweight' 'bold' } ...
         { 'Style', 'checkbox'  , 'string', '' 'tag' 'cb_owrt' 'callback' cb_owrt } ...
         { 'Style', 'text'      , 'string', 'Overwrite it in memory (set=yes; unset=create a new dataset)' } {} ...
         { } ...
         { 'Style', 'checkbox'  , 'string', '', 'callback', cb_save2 'value' value_save2  'tag' 'cb_save2' } ...
         { 'Style', 'text'      , 'string', 'Save it as file:' } ...
         { 'Style', 'edit'      , 'string', filename, 'tag', 'saveedit2' 'userdata' 'saveedit2'  'enable' enable_save2 } ...
         { 'Style', 'pushbutton', 'string', 'Browse', 'callback', comsave2 'userdata' 'saveedit2'  'enable' enable_save2 } };
    
    % show GUI (do not return if old dataset has to saved or overwritten)
    % -------------------------------------------------------------------
    cont = 1;
    while cont
        [result userdat] = inputgui( 'geometry', geometry, 'uilist', uilist, 'helpcom', 'pophelp(''pop_newset'');', ...
                                     'title', fastif(isempty(EEG.data), 'Import dataset info -- pop_newset()', 'Edit dataset info -- pop_newset()'), ...
                                     'userdata', userdat, 'geomvert', geomvert);
        cont = 0;
        if length(result) ~= 0 & overwrite_or_save
            if result{5} & ~result{4} % save but not overwrite
                if isempty(result{6})
                    warndlg2(strvcat('Error: You must enter a name for the old dataset!',' ', ...
                           'The old dataset has to be saved because,', ...
                           'according to your memory option, only one full ', ...
                           'dataset can be retained in memory. Thus, you either', ...
                           'have to save or overwrite the old dataset')); ... 
                    cont = 1;
                end;
            end;
        end;
    end;
    
    % decode parameters
    % -----------------
    if length(result) == 0,
		args = { 'retrieve', CURRENTSET }; % cancel
	else 
		args = { 'setname', result{1} };
		if result{2} 
            if ~isempty(result{3}) 
                args = { args{:} 'savenew', result{3} };
            else
                disp('Warning: no file name given for new dataset (the dataset will not be saved on disk)');
            end;
		end;
		if ~strcmp(EEG.comments, userdat)
			args = { args{:} 'comments', userdat };
		end;
        if result{4}
			args = { args{:} 'overwrite' 'on' };
        end;
		if result{5} 
            if ~isempty(result{6}) 
                args = { args{:} 'saveold', result{6} };
            else
                disp('Warning: no file name given for old dataset (the dataset will not be saved on disk)');
            end;
		end;
    end;
elseif length(EEG) > 1
    % processing multiple datasets
    % ----------------------------
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET ); % it is possible to undo the operation here
    com = '[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET );';
    return;
else
    % no interactive inputs
    args = varargin;
end;

% assigning values
% ----------------
overWflag    = 0;
EEG.saved = 'no';
for ind = 1:2:length(args)
    switch lower(args{ind})
	 case 'setname'   , EEG.setname = args{ind+1}; EEG = eeg_hist(EEG, [ 'EEG.setname=''' EEG.setname ''';' ]);
	 case 'comments'  , EEG.comments = args{ind+1};
	 case 'retrieve'  , if ~isempty(ALLEEG) 
                            EEG = eeg_retrieve(ALLEEG, args{ind+1}); 
                        else
                            EEG = eeg_emptyset;
                        end;
                        overWflag = 1; com = ''; return;
	 case { 'save' 'savenew' }, [filepath filename ext] = fileparts( args{ind+1} );
                        EEG.saved = 'yes';
                        EEG = pop_saveset(EEG, [ filename ext ], filepath);
	 case 'saveold',    [filepath filename ext] = fileparts( args{ind+1} );
                        ALLEEG(CURRENTSET).saved = 'yes';
                        EEG = pop_saveset(ALLEEG(CURRENTSET), [ filename ext ], filepath);
	 case 'overwrite' , if strcmpi(args{ind+1}, 'on') | strcmpi(args{ind+1}, 'yes')
                            overWflag = 1; % so it can be done at the end
                        end;
	 otherwise, error(['pop_newset error: unrecognized key ''' args{ind} '''']); 
    end;
end;
if overWflag
	[ALLEEG, EEG] = eeg_store( ALLEEG, EEG, CURRENTSET);
else
	[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0); % 0 means that it is saved on disk
end;
	
% generate the output command
% ---------------------------
com = sprintf( '[ALLEEG EEG %s] = pop_newset(ALLEEG, EEG, %s, %s);', inputname(3), inputname(3), vararg2str(args));
return;

function num = popask( text )
	 ButtonName=questdlg2( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end;
