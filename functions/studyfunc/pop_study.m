%
% pop_study()  - create a new STUDY set structure defining a group of related EEG datasets. 
%                      The STUDY set also contains information about each of the datasets: 
%                      the subject code, subject group, exp. condition and session. This can be 
%                      provided interactively in a pop-up window or be automatically filled 
%                      in by the function. Defaults: Assume a different subject for each 
%                      dataset and only one condition; leave subject group and session fields 
%                      empty. Additional STUDY information about the STUDY name, task and 
%                      miscellaneous notes can also be saved in the STUDY structure.
% Usage:  
%                        >> [ STUDY ALLEEG] = pop_study; % specify sets interactively
%                        >> [ STUDY ALLEEG] = pop_study(ALLEEG); % included the loaded
%                                                                % EEG datasets
% Optional Inputs:
%   ALLEEG               - vector of EEG dataset structures to be included in the STUDY. 
%
% Outputs:
%   STUDY                - a new STUDY set comprising some or all of the datasets in 
%                          ALLEEG, plus other information about the experiments. 
%   ALLEEG               - an updated ALLEEG structure including the STUDY datasets. 
%
% Graphic interface buttons:
%  "STUDY set name"      - [edit box] a name for the STUDY structure {default: ''}   
%  "STUDY set task name" - [edit box] name for the task performed by the subject {default: ''}   
%  "STUDY set notes"     - [edit box] notes about the experiment, the datasets, the STUDY, 
%                          or anything else to store with the rest of the STUDY information 
%                          {default: ''}   
%  "subject"             - [edit box] subject code associated with the dataset. If no 
%                          subject code is provided, each dataset will assumed to be from 
%                          a different subject {default: 'S1', 'S2', ..., 'Sn'}
%  "session"             - [edit box] dataset session. If no session information is
%                          provided, all datasets that belong to one subject are assumed to
%                          have been recorded within one session {default: []}
%  "condition"           - [edit box] dataset condition. If no condition code is provided,
%                          all datasets are assumed to be from the same condition {default:[]}
%  "group"               - [edit box] the subject group the dataset belongs to. If no group
%                          is provided, all subjects and datasets are assumed to belong to
%                          the same group. {default: []}
%  "Save this STUDY set to disk file" - [check box] If checked, save the new STUDY set
%                          structure to disk. If no filename is provided, a window will 
%                          pop up to ask for the file name. 
%
%  See also  create_study(), pop_loadstudy(), pop_preclust(), eeg_preclust(), pop_clust()
%
% Authors: Arnaud Delorme, Hilit Serby, Scott Makeig, SCCN, INC, UCSD, July 22, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, July 22, 2005, arno@sccn.ucsd.edu
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

% Coding notes: Useful information on functions and global variables used.

% $Log: not supported by cvs2svn $
% Revision 1.28  2006/02/23 18:56:10  scott
% help msg -sm
%
% Revision 1.27  2006/02/23 00:43:05  arno
% trying to retain selection in listbox
%
% Revision 1.26  2006/02/23 00:31:05  arno
% better message
%
% Revision 1.25  2006/02/23 00:16:52  arno
% fix component selection
%
% Revision 1.24  2006/02/23 00:06:20  arno
% call to editstudy
%
% Revision 1.23  2006/02/23 00:04:54  arno
% nothing
%
% Revision 1.22  2006/02/23 00:01:45  arno
% dipole selection
%
% Revision 1.21  2006/02/22 21:50:40  arno
% fixing delclust option
%
% Revision 1.20  2006/02/09 22:44:17  arno
% text edit
%
% Revision 1.19  2006/02/09 21:49:15  arno
% adding rmclust options
%
% Revision 1.18  2006/02/09 21:39:28  arno
% automatically update other datasets with same session and subject
%
% Revision 1.17  2006/02/09 20:07:09  arno
% OK for creating new dataset, consistency etc...
%
% Revision 1.16  2006/02/09 00:40:35  arno
% adding new button to delete cluster information
%
% Revision 1.15  2006/02/08 23:23:40  arno
% allow comps input
%
% Revision 1.14  2006/02/08 23:15:59  arno
% don't know
%
% Revision 1.13  2006/02/03 21:47:18  arno
% fixing text
%
% Revision 1.12  2006/02/03 20:48:35  arno
% typo
%
% Revision 1.11  2006/02/03 20:47:03  arno
% adding savedat option
%
% Revision 1.10  2006/02/03 20:36:41  arno
% saving log
%

function [STUDY, ALLEEG, com]  = pop_study(STUDY, ALLEEG, varargin)

com = '';

if nargin < 1
    help pop_study;
    return;
end

% type of call (gui, script or internal)
% --------------------------------------
mode = 'internal_command';
if ~isstr(STUDY) %intial settings
    mode = 'script';
    if nargin > 2
        for index = 1:length(varargin)
            if isstr(varargin{index})
                if strcmpi(varargin{index}, 'gui')
                    mode = 'gui';
                end;
            end;
        end;
    end;
end;

if isempty(STUDY)
    newstudy       = 1;
    STUDY.name     = '';
    STUDY.task     = '';
    STUDY.notes    = '';
    STUDY.filename = '';
    STUDY.cluster  = [];
else
    newstudy       = 0;
end;

if strcmpi(mode, 'script') % script mode
    [STUDY ALLEEG] = editstudy(STUDY, ALLEEG, varargin{:});
    return;
elseif strcmpi(mode, 'gui') % GUI mode
    % show warning if necessary
    % -------------------------
    if isreal(ALLEEG)
        if ALLEEG == 0
            res = questdlg2( strvcat('Datasets currently loaded will be removed from EEGLAB memory.', ...
                                     'Are you sure you want to pursue and ignore currenlty loaded datasets?'), ...
                                     'Discard datasets loaded in EEGLAB?', 'Cancel', 'Yes', 'Yes');
            if strcmpi(res, 'cancel'), return; end;
        end;
        ALLEEG = [];
    end;
    
    % set initial datasetinfo
    % -----------------------
    if isfield(STUDY, 'datasetinfo')
        datasetinfo = STUDY.datasetinfo;
        different = 0;
        for k = 1:length(ALLEEG)
            if ~strcmpi(datasetinfo(k).filename, ALLEEG(k).filename), different = 1; break; end;
            if ~strcmpi(datasetinfo(k).subject,   ALLEEG(k).subject),   different = 1; break; end;
            if ~strcmpi(datasetinfo(k).condition, ALLEEG(k).condition), different = 1; break; end;
            if ~strcmpi(char(datasetinfo(k).group), char(ALLEEG(k).group)),     different = 1; break; end;              
            if datasetinfo(k).session ~= ALLEEG(k).session,             different = 1; break; end;
        end
        if different
            info = 'from_STUDY_different_from_ALLEEG';
        else
            info = 'from_STUDY';
        end;
        if ~isfield(datasetinfo, 'comps');
            datasetinfo(1).comps = [];
        end;
    else
        info = 'from_ALLEEG';
        if length(ALLEEG) > 0
            datasetinfo(length(ALLEEG)).filename  = [];
            datasetinfo(length(ALLEEG)).filepath  = [];
            datasetinfo(length(ALLEEG)).subject   = [];
            datasetinfo(length(ALLEEG)).session   = [];
            datasetinfo(length(ALLEEG)).condition = [];
            datasetinfo(length(ALLEEG)).group     = [];     
            for k = 1:length(ALLEEG)
                datasetinfo(k).filename  = ALLEEG(k).filename;   
                datasetinfo(k).filepath  = ALLEEG(k).filepath;   
                datasetinfo(k).subject   = ALLEEG(k).subject;
                datasetinfo(k).session   = ALLEEG(k).session;
                datasetinfo(k).condition = ALLEEG(k).condition;
                datasetinfo(k).group     = ALLEEG(k).group;                    
            end
            if ~isfield(datasetinfo, 'comps');
                datasetinfo(1).comps = [];
            end;
        else
            datasetinfo = [];
        end;
    end;

    nextpage = 'pop_study(''nextpage'', gcbf);';
    prevpage = 'pop_study(''prevpage'', gcbf);';
    delset   = 'pop_study(''clear'',    gcbf, get(gcbo, ''userdata''));';
    loadset  = 'pop_study(''load'',     gcbf, get(guiind, ''userdata''), get(guiind, ''string'')); clear guiind;';
    loadsetedit  = [ 'guiing = gcbo;' loadset ];
    subcom    = 'pop_study(''subject''  , gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    sescom    = 'pop_study(''session''  , gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    condcom   = 'pop_study(''condition'', gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    grpcom    = 'pop_study(''group''    , gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    grpcom    = 'pop_study(''component'', gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    cb_del    = 'pop_study(''delclust'' , gcbf, ''showwarning'');';
    cb_dipole = 'pop_study(''dipselect'', gcbf, ''showwarning'');';

	browsestudy = [ '[filename, filepath] = uiputfile2(''*.study'', ''Use exsiting STUDY set to import dataset information -- pop_study()''); ' ... 
                      'set(findobj(''parent'', gcbf, ''tag'', ''usestudy_file''), ''string'', [filepath filename]);' ];
	saveSTUDY = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''save''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
	browsesave = [ '[filename, filepath] = uiputfile2(''*.study'', ''Save STUDY with .study extension -- pop_clust()''); ' ... 
                   'if filename ~= 0,' ...
                   '   set(findobj(''parent'', gcbf, ''tag'', ''studyfile''), ''string'', [filepath filename]);' ...
                   'end;' ...
                   'clear filename filepath;' ];
	
    texthead = fastif(newstudy, 'Create a new STUDY set', 'Edit STUDY set information');
	guispec = { ...
        {'style' 'text' 'string' texthead 'FontWeight' 'Bold' 'FontSize' 11 'HorizontalAlignment' 'center'} ...
        {'style' 'text' 'string' 'STUDY set name:' } { 'style' 'edit' 'string' STUDY.name 'tag' 'study_name' } ...
        {'style' 'text' 'string' 'STUDY set task name:' } { 'style' 'edit' 'string' STUDY.task 'tag' 'study_task' } ...
        {'style' 'text' 'string' 'STUDY set notes:' } { 'style' 'edit' 'string' STUDY.notes 'tag' 'study_notes' } {}...
        {} ...
        {'style' 'text' 'string' 'dataset filename' 'userdata' 'addt'} {'style' 'text' 'string' 'browse' 'userdata' 'addt'} ...
        {'style' 'text' 'string' 'subject'    'userdata' 'addt'} ...
        {'style' 'text' 'string' 'session'    'userdata' 'addt'} ...
        {'style' 'text' 'string' 'condition'  'userdata' 'addt'} ...
        {'style' 'text' 'string' 'group'      'userdata' 'addt'} ...
        {'style' 'pushbutton' 'string' 'components' 'userdata' 'addt' 'callback' cb_dipole } ...
        {} };
	guigeom = { [1] [1 2] [1 2] [1 2] [1] [0.2 1.05 0.35 0.4 0.35 0.6 0.4 0.6 0.3]};
	
    % create edit boxes
    % -----------------
    for index = 1:10
        guigeom = { guigeom{:} [0.2 1 0.2 0.5 0.2 0.5 0.5 0.5 0.3] };
        select_com = ['[inputname, inputpath] = uigetfile2(''*.set*;*.SET*'', ''Choose dataset to add to STUDY -- pop_study()'');'...
                      'if inputname ~= 0,' ...
                      '   guiind = findobj(''parent'', gcbf, ''tag'', ''set' int2str(index) ''');' ...
                      '   set( guiind,''string'', [inputpath inputname]);' ...
                          loadset ...
                      'end;'];
        numstr = int2str(index);
        guispec = { guispec{:} ...
        {'style' 'text'       'string' numstr  'tag' [ 'num'   int2str(index) ] 'userdata' index } ...
        {'style' 'edit'       'string' ''      'tag' [ 'set'   int2str(index) ] 'userdata' index 'callback' loadsetedit} ...
        {'style' 'pushbutton' 'string' '...'   'tag' [ 'brw'   int2str(index) ] 'userdata' index 'Callback' select_com} ...
        {'style' 'edit'       'string' ''      'tag' [ 'sub'   int2str(index) ] 'userdata' index 'Callback' subcom} ...
        {'style' 'edit'       'string' ''      'tag' [ 'sess'  int2str(index) ] 'userdata' index 'Callback' sescom} ...
        {'style' 'edit'       'string' ''      'tag' [ 'cond'  int2str(index) ] 'userdata' index 'Callback' condcom} ...
        {'style' 'edit'       'string' ''      'tag' [ 'group' int2str(index) ] 'userdata' index 'Callback' grpcom} ...
        {'style' 'pushbutton' 'string' 'All'   'tag' [ 'comps' int2str(index) ] 'userdata' index 'Callback' grpcom} ...
        {'style' 'pushbutton' 'string' 'CLear' 'tag' [ 'clear' int2str(index) ] 'userdata' index 'callback' delset} };
    end;
    
    if strcmpi(info, 'from_STUDY_different_from_ALLEEG')
        text1    = 'The dataset info (condition, group, ...) is different from info stored in study. [set] Overwrite dataset info.';
        value_cb = 0;
    else
        text1    = 'Update dataset info (unset=local to study). If datasets are stored on disk, they will be overwritten.';
        value_cb = 1;
    end;
    guispec = { guispec{:} ...
                {'style' 'text'       'string'  'Important note: Removed datasets will not be saved prior to being delete from EEGLAB memory' } ...
                {} ...
                {'style' 'pushbutton' 'string'  '<'      'Callback' prevpage 'userdata' 'addt'} ...
                {'style' 'text'       'string'  'Page 1' 'tag' 'page' 'horizontalalignment' 'center' } ... 
                {'style' 'pushbutton' 'string'  '>'      'Callback' nextpage 'userdata' 'addt'} {} ...
                {} ...
                {'style' 'checkbox'   'value'   value_cb 'tag' 'copy_to_dataset' } ...
                {'style' 'text'       'string'  text1 } ...
                {'style' 'checkbox'   'value'   0        'tag' 'delclust' 'callback' cb_del } ...
                {'style' 'text'       'string'  'Delete cluster information (allow loading new datasets, set new starting components for clustering, ...)' } ...
                {'style' 'text'       'string'  'Save this STUDY set to disk file'} ...
                {'style' 'edit'       'string'  ''       'tag' 'studyfile'                        'userdata' 'save'} ...
                {'style' 'pushbutton' 'string'  '...'    'tag' 'browsesave' 'Callback' browsesave 'userdata' 'save'} {} };
	guigeom = { guigeom{:} [1] [1 0.2 0.3 0.2 1] [1] [0.14 3] [0.14 3] [1 1.5 0.3] [1]};

    if ~isempty(STUDY.filename)
        guispec{end-3} = {'style' 'checkbox' 'string' '' 'value' 1 'tag' 'studyfile' };
        guispec{end-2} = {'style' 'text'     'string' 'Resave study. Uncheck and use menu File > Save study as to save under a new filename'};
        guispec{end-1} = {};
        guigeom{end-1} = [0.08 1.5 0.1];
    end;
	
    fig_arg{1} = ALLEEG;      % datasets         
    fig_arg{2} = datasetinfo; % datasetinfo
    fig_arg{3} = 1;           % page
    fig_arg{4} = {};          % all commands
    fig_arg{5} = (length(STUDY.cluster) > 1); % are cluster present
    fig_arg{6} = STUDY; % are cluster present
    
    % generate GUI
    % ------------
    optiongui = { 'geometry', guigeom, ...
                  'uilist'  , guispec, ...
                  'helpcom' , 'pophelp(''pop_study'')', ...
                  'title'   , 'Create a new STUDY set -- pop_study()', ...
                  'userdata', fig_arg, ...
                  'eval'    , 'pop_study(''delclust'', gcf); pop_study(''redraw'', gcf);' };
	[result, userdat2, strhalt, outstruct] = inputgui( 'mode', 'noclose', optiongui{:});
    if isempty(result), return; end;
    if ~isempty(get(0, 'currentfigure')) currentfig = gcf; end;
    
    while test_wrong_parameters(currentfig)
    	[result, userdat2, strhalt, outstruct] = inputgui( 'mode', currentfig, optiongui{:});
        if isempty(result), return; end;
    end;
    close(currentfig);

    % convert GUI selection to options
    % --------------------------------
    allcom = simplifycom(userdat2{4});
    options = {};
    if ~strcmpi(result{1}, STUDY.name ), options = { options{:} 'name'        result{1} }; end;
    if ~strcmpi(result{2}, STUDY.task ), options = { options{:} 'task'        result{2} }; end;
    if ~strcmpi(result{3}, STUDY.notes), options = { options{:} 'notes'       result{3} }; end;
    if ~isempty(result{4}),              options = { options{:} 'commands'    allcom    }; end;
    if isnumeric(outstruct(1).studyfile)
        if outstruct(1).studyfile == 1,  options = { options{:} 'resave'      'on' }; end;
    else
        if ~isempty(outstruct(1).studyfile), options = { options{:} 'filename' outstruct(1).studyfile }; end;
    end;
    if outstruct(1).copy_to_dataset == 1
         options = { options{:} 'updatedat' 'on' };
         eeglab_options;
         if option_storedisk
             options = { options{:} 'savedat' 'on' };
         end;
    else options = { options{:} 'updatedat' 'off' };
    end;
    if outstruct(1).delclust == 1
        options = { options{:} 'rmclust' 'on' };
    end;
    
    % run command and create history
    % ------------------------------
    [STUDY ALLEEG] = editstudy(STUDY, ALLEEG, options{:});
    com = sprintf( '[%s %s] = pop_study( %s, %s, %s );', inputname(1), inputname(2), ...
                                                         inputname(1), inputname(2), vararg2str(options(3:end)) );
    
else % internal command
    
    com = STUDY;
    hdl = ALLEEG; %figure handle

    % userdata info
    % -------------
    userdat     = get(hdl, 'userdata');    
    ALLEEG      = userdat{1};
    datasetinfo = userdat{2};
    page        = userdat{3};
    allcom      = userdat{4};
    clusterpresent = userdat{5};
    STUDY       = userdat{6};

    switch  com
        case 'subject'
            guiindex  = varargin{1};
            realindex = guiindex+(page-1)*10;
            
            datasetinfo(realindex).subject   = varargin{2};
            if get(findobj(hdl, 'tag', 'copy_to_dataset'), 'value')
                ALLEEG(realindex).subject    = varargin{2};
            end;
            allcom = { allcom{:} { 'index' realindex 'subject' varargin{2} } };
            userdat{1} = ALLEEG;
            userdat{2} = datasetinfo;
            userdat{4} = allcom;
            set(hdl, 'userdata', userdat);            

        case 'session'
            guiindex  = varargin{1};
            realindex = guiindex+(page-1)*10;
            datasetinfo(realindex).session   = str2num(varargin{2});
            if get(findobj(hdl, 'tag', 'copy_to_dataset'), 'value')
                ALLEEG(realindex).session    = str2num(varargin{2});
            end;
            allcom = { allcom{:} { 'index' realindex 'session' str2num(varargin{2}) } };
            userdat{1} = ALLEEG;
            userdat{2} = datasetinfo;
            userdat{4} = allcom;
            set(hdl, 'userdata', userdat);            

        case 'group'
            guiindex  = varargin{1};
            realindex = guiindex+(page-1)*10;
            datasetinfo(realindex).group   = varargin{2};
            if get(findobj(hdl, 'tag', 'copy_to_dataset'), 'value')
                ALLEEG(realindex).group    = varargin{2};
            end;
            allcom = { allcom{:} { 'index' realindex 'group' varargin{2} } };
            userdat{1} = ALLEEG;
            userdat{2} = datasetinfo;
            userdat{4} = allcom;
            set(hdl, 'userdata', userdat);            

        case 'condition'
            guiindex  = varargin{1};
            realindex = guiindex+(page-1)*10;
            datasetinfo(realindex).condition   = varargin{2};
            if get(findobj(hdl, 'tag', 'copy_to_dataset'), 'value')
                ALLEEG(realindex).conditon     = varargin{2};
            end;
            allcom = { allcom{:} { 'index' realindex 'condition' varargin{2} } };
            userdat{1} = ALLEEG;
            userdat{2} = datasetinfo;
            userdat{4} = allcom;
            set(hdl, 'userdata', userdat);            
            
        case 'dipselect'
            STUDY.datasetinfo = datasetinfo;
            
            res = inputdlg2( { strvcat('Enter threshold max residual variance in % to pre-select components', ...
                                       'Note this will delete all existing clusters') }, ...
                             'Pre-select components', 1, { '15' } );
            
            STUDY = editstudy(STUDY, ALLEEG, 'commands', { 'dipselect' str2num(res{1})/100 'return' });
            allcom = { allcom{:} { 'dipselect' str2num(res{1})/100 } };
            datasetinfo   = STUDY.datasetinfo;
            
            userdat{2} = datasetinfo;
            userdat{4} = allcom;
            set(hdl, 'userdata', userdat);            
            pop_study('redraw', hdl);

        case 'component'
            guiindex  = varargin{1};
            realindex = guiindex+(page-1)*10;
            
            for index = 1:size(ALLEEG(realindex).icaweights,1)
                complist{index} = [ 'IC ' int2str(index) ];
            end;
            [tmps,tmpv] = listdlg2('PromptString', 'Select components', 'SelectionMode', ...
                                    'multiple', 'ListString', strvcat(complist), 'value', datasetinfo(realindex).comps);
            if tmpv ~= 0 % no cancel                
                
                % find other subjects with the same session
                % -----------------------------------------
                for index = 1:length(datasetinfo)
                    if realindex == index | (strcmpi(datasetinfo(index).subject, datasetinfo(realindex).subject) & ...
                                ~isempty(datasetinfo(index).subject) & ...
                                isequal( datasetinfo(index).session, datasetinfo(realindex).session ) )
                        datasetinfo(index).comps = tmps;
                        allcom = { allcom{:} { 'index' index 'comps' tmps } };
                        set(findobj('tag', [ 'comps' int2str(index) ]), ...
                            'string', formatbut(tmps), 'horizontalalignment', 'left');
                    end;
                end;
                
            end;
            userdat{2} = datasetinfo;
            userdat{4} = allcom;
            set(hdl, 'userdata', userdat);            
            pop_study('redraw', hdl);

        case 'clear'
            guiindex  = varargin{1};
            realindex = guiindex+(page-1)*10;
            
            datasetinfo(realindex).filename  = '';   
            datasetinfo(realindex).filepath  = '';   
            datasetinfo(realindex).subject   = '';
            datasetinfo(realindex).session   = [];
            datasetinfo(realindex).condition = '';
            datasetinfo(realindex).group     = '';                    
            datasetinfo(realindex).comps     = [];                    

            allcom = { allcom{:} { 'remove' realindex } };
            userdat{1} = ALLEEG;
            userdat{2} = datasetinfo;
            userdat{4} = allcom;
            set(hdl, 'userdata', userdat);
            pop_study('redraw', hdl);

        case 'nextpage'
            userdat{3} = page+1;
            set(hdl, 'userdata', userdat);
            pop_study('redraw', hdl);

        case 'prevpage'
            userdat{3} = max(1,page-1);
            set(hdl, 'userdata', userdat);
            pop_study('redraw', hdl);

        case 'load'
            guiindex  = varargin{1};
            filename  = varargin{2};
            realindex = guiindex+(page-1)*10;

            % load dataset
            % ------------
            TMPEEG = pop_loadset('filename', filename, 'loadmode', 'info');
            ALLEEG = eeg_store(ALLEEG, eeg_checkset(TMPEEG), realindex);
            
            % update datasetinfo structure
            % ----------------------------
            datasetinfo(realindex).filename  = ALLEEG(realindex).filename;   
            datasetinfo(realindex).filepath  = ALLEEG(realindex).filepath;   
            datasetinfo(realindex).subject   = ALLEEG(realindex).subject;
            datasetinfo(realindex).session   = ALLEEG(realindex).session;
            datasetinfo(realindex).condition = ALLEEG(realindex).condition;
            datasetinfo(realindex).group     = ALLEEG(realindex).group;                    
            datasetinfo(realindex).comps     = [];                    

            allcom = { allcom{:} { 'index' realindex 'load' filename } };
            userdat{1} = ALLEEG;
            userdat{2} = datasetinfo;
            userdat{4} = allcom;
            set(hdl, 'userdata', userdat);
            pop_study('redraw', hdl);

         case 'delclust'
            if clusterpresent
                if ~get(findobj(hdl, 'tag', 'delclust'), 'value')
                    for k = 1:10
                        set(findobj('parent', hdl, 'tag',['set'   num2str(k)]), 'style', 'text');
                        set(findobj('parent', hdl, 'tag',['comps' num2str(k)]), 'enable', 'off');
                        set(findobj('parent', hdl, 'tag',['sess'  num2str(k)]), 'enable', 'off');
                        set(findobj('parent', hdl, 'tag',['brw'   num2str(k)]), 'enable', 'off');
                    end;
                else
                    for k = 1:10
                        set(findobj('parent', hdl, 'tag',['set'   num2str(k)]), 'style', 'edit');
                        set(findobj('parent', hdl, 'tag',['comps' num2str(k)]), 'enable', 'on');
                        set(findobj('parent', hdl, 'tag',['sess'  num2str(k)]), 'enable', 'on');
                        set(findobj('parent', hdl, 'tag',['brw'   num2str(k)]), 'enable', 'on');
                    end;
                end;
            else 
                set(findobj(hdl, 'tag', 'delclust'), 'value', 0)
                if nargin > 2
                    warndlg2('No cluster present');
                end;
            end;
         
        case 'redraw'
            for k = 1:10
                kk = k+(page-1)*10; % real index
                if kk > length(datasetinfo)
                    set(findobj('parent', hdl, 'tag',['num' num2str(k)]), 'string', int2str(kk));
                    set(findobj('parent', hdl, 'tag',['set' num2str(k)]), 'string', '');
                    set(findobj('parent', hdl, 'tag',['sub' num2str(k)]), 'string','');
                    set(findobj('parent', hdl, 'tag',['sess' num2str(k)]), 'string','');
                    set(findobj('parent', hdl, 'tag',['cond' num2str(k)]), 'string','');
                    set(findobj('parent', hdl, 'tag',['comps' num2str(k)]), 'string','');
                    set(findobj('parent', hdl, 'tag',['group' num2str(k)]), 'string','');
                else
                    set(findobj('parent', hdl, 'tag',['num' num2str(k)]), 'string', int2str(kk));
                    set(findobj('parent', hdl, 'tag',['set' num2str(k)]), 'string', fullfile(datasetinfo(kk).filepath, datasetinfo(kk).filename));
                    set(findobj('parent', hdl, 'tag',['sub' num2str(k)]), 'string', datasetinfo(kk).subject);
                    set(findobj('parent', hdl, 'tag',['sess' num2str(k)]), 'string', int2str(datasetinfo(kk).session));
                    set(findobj('parent', hdl, 'tag',['cond' num2str(k)]), 'string', datasetinfo(kk).condition);
                    set(findobj('parent', hdl, 'tag',['comps' num2str(k)]), 'string', formatbut(datasetinfo(kk).comps));
                    set(findobj('parent', hdl, 'tag',['group' num2str(k)]), 'string', datasetinfo(kk).group);
                end;
            end
            if page<10
                 pagestr =  [ ' Page ' int2str(page) ];
            else pagestr =  [ 'Page ' int2str(page) ];
            end;
            set(findobj('parent', hdl, 'tag','page'), 'string', pagestr );
    end
end;

% remove empty elements in allcom
% -------------------------------
function allcom = simplifycom(allcom);

    for index = length(allcom):-1:2
        if strcmpi(allcom{index}{1}, 'index') & strcmpi(allcom{index-1}{1}, 'index')
            if allcom{index}{2} == allcom{index-1}{2}
                allcom{index-1}(end+1:end+2) = allcom{index}(3:4);
                allcom(index) = [];
            end;
        end;
    end;
    
% test for wrong parameters
% -------------------------
function bool = test_wrong_parameters(hdl)
    userdat     = get(hdl, 'userdata');    
    datasetinfo = userdat{2};

    bool = 0;
    for index = 1:length(datasetinfo)
        if ~isempty(datasetinfo(index).filename)
            if isempty(datasetinfo(index).subject) & bool == 0
                bool = 1; warndlg2('All dataset must have a subject name or code', 'Error');
            end;
        end;
    end;
    
    anysession   = any(~cellfun('isempty', { datasetinfo.session }));
    allsession   = all(~cellfun('isempty', { datasetinfo.session }));
    anycondition = any(~cellfun('isempty', { datasetinfo.condition }));
    allcondition = all(~cellfun('isempty', { datasetinfo.condition }));
    anygroup     = any(~cellfun('isempty', { datasetinfo.group }));
    allgroup     = all(~cellfun('isempty', { datasetinfo.group }));

    if anygroup & ~allgroup
         bool = 1; warndlg2('If one dataset has a group, they must all have one', 'Error');
    end;
    if anycondition & ~allcondition
         bool = 1; warndlg2('If one dataset has a condition, they must all have one', 'Error');
    end;
    if anysession & ~allsession
         bool = 1; warndlg2('If one dataset has a session index, they must all have one', 'Error');
    end;
    
function strbut = formatbut(complist)
    if isempty(complist)
        strbut = 'All';
    else
        if length(complist) > 5, strbut = [ ' ' int2str(complist(1:4)) ' ...' ];
        else                     strbut = [ ' ' int2str(complist) ];
        end;
    end;
