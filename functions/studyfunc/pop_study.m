% pop_study()  - create a new STUDY set structure defining a group of related EEG 
%                      datasets. the STUDY set also contains information about each of the
%                      datasets: the subject, group, condition and session. This can be 
%                      provided interactively in a pop-up window or be automatically filled 
%                      in by the function. Default: assume a different subject for each 
%                      dataset and only one condition, with dataset group and session fields 
%                      empty. Additional STUDY information about the STUDY name, task and 
%                      miscellaneous notes can also be saved in the STUDY structure.
% Usage:  
%                        >> [ STUDY ALLEEG] = pop_study; % specify sets interactively
%                        >> [ STUDY ALLEEG] = pop_study(ALLEEG); % included loaded
%                                                                      % EEG datasets
% Optional Inputs:
%   ALLEEG               - vector of EEG dataset structures to be included in the STUDY. 
%
% Outputs:
%   STUDY                - a new STUDY set comprising some or all of the datasets in 
%                          ALLEEG, plus other information on the experiments. 
%   ALLEEG               - an updated ALLEEG structure including the STUDY datasets. 
%
% Graphic interface buttons:
%  "STUDY set name"      - [edit box] a name for the STUDY structure {default: ''}   
%  "STUDY set task name" - [edit box] Name of the task performed by the subject {default: ''}   
%  "STUDY set notes"     - [edit box] Notes about the experiment, the datasets, the STUDY, 
%                          or anything else to keep with the rest of the STUDY information 
%                          {default: ''}   
%  "Include the datasets loaded in eeglab" - [check box] If checked and optional input ALLEEG 
%                          is provided, its datasets are added to the STUDY. The datasets will 
%                          be available in the "# datasets in STUDY" list box, where they can 
%                          be edited. 
%  "Include datasets and information from another STUDY" - [check box] If checked, an existing 
%                          STUDY set (full path) filename must be provided in the associated 
%                          edit box. The dataset information (file name/path, subject, 
%                          condition, etc.) will be copied from the existing STUDY set. 
%                          Those datasets will be added to datasets that were imported from 
%                          an ALLEEG EEG structure (see option above) or to any manually 
%                          added datasets (see option below).  
%  "Include these datasets" - [check box] If checked, datasets and dataset information can be 
%                          manually added to the STUDY set. The dataset full path and filename
%                          must be provided or selected using he bowse button... You may then
%                          provide subject, session, condition, and subject group information
%                          for the dataset. If you leave those fields empty, default values
%                          will be used. If you do provide this information for one dataset, 
%                          you will need to enter it for all datasets you add to the STUDY.
%  "subject"             - [edit box] The subject associated with the dataset. If no subject
%                          information is provided, each dataset will assumed to be from a
%                          different subject {defaults: 'S1', 'S2', ..., 'Sn'}
%  "session"             - [edit box] The dataset session. If no session information is
%                          provided, all datasets that belong to a subject are assumed to
%                          have been recorded in one session {default: []}
%  "condition"           - [edit box] The dataset condition. If no condition is provided,
%                          all datasets are assumed to have the same condition {default:[]}
%  "group"               - [edit box] the subject group the dataset belongs to. If no group
%                          is provided all datasets and subjects are assumed to belong to
%                          the same group. {default: []}
%  "Add these and more"  - [pushbutton] Add the entered datasets to the STUDY and clear the
%                          input fields so that additional datasets can be added. The added
%                          datasets will appear in the "# datasets in STUDY" list box.
%  "Add these"           - [pushbutton] add the datasets entered to the STUDY. The added
%                          datasets will appear in the "# datasets in STUDY" list box and
%                          no further dataset entry will be possible. 
%  "# datasets in STUDY" - [list box] A list of the dataset file names to add to the
%                          STUDY, both those entered manually and those in the ALLEEG 
%                          structure. Any of these datasets may be deleted from the list 
%                          using the 'Delete' button, or their information (subject, 
%                          session, etc.) may be edited.
%  "Save this STUDY set to disk file" - [check box] If checked, save the new STUDY set
%                          structure to disk. If no filename is provided, a window will 
%                          pop up to ask for the file name. 
%
%  See also  create_study(), pop_loadstudy(), pop_preclust(), eeg_preclust(), pop_clust()
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, July 22, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, July 22, 2005, hilit@sccn.ucsd.edu
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
    STUDY.name  = '';
    STUDY.task  = '';
    STUDY.notes = '';
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
    datasetinfo.filename  = [];
    datasetinfo.subject   = [];
    datasetinfo.session   = [];
    datasetinfo.condition = [];
    datasetinfo.group     = [];     
    for k = 1:length(ALLEEG)
        datasetinfo(k).filename  = fullfile(ALLEEG(k).filepath, ALLEEG(k).filename);   
        datasetinfo(k).subject   = ALLEEG(k).subject;
        datasetinfo(k).session   = ALLEEG(k).session;
        datasetinfo(k).condition = ALLEEG(k).condition;
        datasetinfo(k).group     = ALLEEG(k).group;                    
    end

    nextpage = 'pop_study(''nextpage'', gcbf);';
    prevpage = 'pop_study(''prevpage'', gcbf);';
    delset   = 'pop_study(''clear'',    gcbf, get(gcbo, ''userdata''));';
    loadset  = 'pop_study(''load'',     gcbf, get(guiind, ''userdata''), get(guiind, ''string'')); clear guiind;';
    loadsetedit  = [ 'guiing = gcbo;' loadset ];
    subcom    = 'pop_study(''subject''  , gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    sescom    = 'pop_study(''session''  , gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    condcom   = 'pop_study(''condition'', gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    grpcom    = 'pop_study(''group''    , gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    
	browsestudy = [ '[filename, filepath] = uiputfile2(''*.study'', ''Use exsiting STUDY set to import dataset information -- pop_study()''); ' ... 
                      'set(findobj(''parent'', gcbf, ''tag'', ''usestudy_file''), ''string'', [filepath filename]);' ];
	saveSTUDY = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''save''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
	browsesave = [ '[filename, filepath] = uiputfile2(''*.study'', ''Save STUDY with .study extension -- pop_clust()''); ' ... 
                   'if filename ~= 0,' ...
                   '   set(findobj(''parent'', gcbf, ''tag'', ''studyfile''), ''string'', [filepath filename]);' ...
                   'end;' ...
                   'clear filename filepath;' ];
	
	guispec = { ...
    {'style' 'text' 'string' 'Create a new STUDY set' 'FontWeight' 'Bold' 'FontSize' 11 'HorizontalAlignment' 'center'} ...
	{'style' 'text' 'string' 'STUDY set name:' } { 'style' 'edit' 'string' STUDY.name 'tag' 'study_name' } ...
	{'style' 'text' 'string' 'STUDY set task name:' } { 'style' 'edit' 'string' STUDY.task 'tag' 'study_task' } ...
	{'style' 'text' 'string' 'STUDY set notes:' } { 'style' 'edit' 'string' STUDY.notes 'tag' 'study_notes' } {}...
	{'style' 'text' 'string' 'datasets (ALLEEG)' 'userdata' 'addt'} {'style' 'text' 'string' 'browse' 'userdata' 'addt'} ...
	{'style' 'text' 'string' 'subject' 'userdata' 'addt'} {'style' 'text' 'string' 'session' 'userdata' 'addt'} ...
	{'style' 'text' 'string' 'condition' 'userdata' 'addt'} {'style' 'text' 'string' 'group'  'userdata' 'addt'} {} };
	guigeom = { [1] [1 2] [1 2] [1 2] [1] [1 0.5 0.5 0.5 0.5 0.5 0.5]};
	
    % create edit boxes
    % -----------------
    for index = 1:10
        guigeom = { guigeom{:} [1 0.5 0.5 0.5 0.5 0.5 0.5] };
        select_com = ['[inputname, inputpath] = uigetfile2(''*.set*;*.SET*'', ''Choose dataset to add to STUDY -- pop_study()'');'...
                      'if inputname ~= 0,' ...
                      '   guiind = findobj(''parent'', gcbf, ''tag'', ''set' int2str(index) ''');' ...
                      '   set( guiind,''string'', [inputpath inputname]);' ...
                          loadset ...
                      'end;'];
        guispec = { guispec{:} ...
        {'style' 'edit'       'string' ''      'tag' [ 'set'   int2str(index) ] 'userdata' index 'callback' loadsetedit} ...
        {'style' 'pushbutton' 'string' '...'                                    'userdata' index 'Callback' select_com} ...
        {'style' 'edit'       'string' ''      'tag' [ 'sub'   int2str(index) ] 'userdata' index 'Callback' subcom} ...
        {'style' 'edit'       'string' ''      'tag' [ 'sess'  int2str(index) ] 'userdata' index 'Callback' sescom} ...
        {'style' 'edit'       'string' ''      'tag' [ 'cond'  int2str(index) ] 'userdata' index 'Callback' condcom} ...
        {'style' 'edit'       'string' ''      'tag' [ 'group' int2str(index) ] 'userdata' index 'Callback' grpcom} ...
        {'style' 'pushbutton' 'string' 'CLear' 'tag' [ 'clear' int2str(index) ] 'userdata' index 'callback' delset} };
    end;
    guispec = { guispec{:} ...
    {'style' 'pushbutton' 'string' '<'      'Callback' prevpage 'userdata' 'addt'} {} ...
    {'style' 'text'       'string' 'Page 1' 'tag' 'page' } {} ... 
	{'style' 'pushbutton' 'string' '>'      'Callback' nextpage 'userdata' 'addt'} ...
    {'style' 'checkbox'   'value' 1 'tag' 'copy_to_dataset' } ...
    {'style' 'text'       'string' 'Update dataset info (unset=local to study)' } ...
    {'style' 'checkbox'   'value' 0 'tag' 'save_dataset' } ...
    {'style' 'text'       'string' 'Resave modified datasets' } ...
    {'style' 'text'       'string' 'Important note: Removed datasets will not be saved prior to being delete from EEGLAB memory' } ...
    { } ... 
    {'style' 'text'       'string' 'Save this STUDY set to disk file'} ...
	{'style' 'edit'       'string' ''       'tag' 'studyfile' 'userdata' 'save'} ...
	{'style' 'pushbutton' 'string' '...'    'tag' 'browsesave' 'Callback' browsesave 'userdata' 'save'} {} };
	guigeom = { guigeom{:} [0.3 1 0.3 1 0.3] [0.2 3] [0.2 3] [1] [1] [1 1.5 0.3] [1]};
	
    fig_arg{1} = ALLEEG;      % datasets         
    fig_arg{2} = datasetinfo; % datasetinfo
    fig_arg{3} = 1;           % page
    fig_arg{4} = {};          % all commands
    
    % generate GUI
    % ------------
    optiongui = { 'geometry', guigeom, ...
                  'uilist'  , guispec, ...
                  'helpcom' , 'pophelp(''pop_study'')', ...
                  'title'   , 'Create a new STUDY set -- pop_study()', ...
                  'userdata', fig_arg, ...
                  'eval'    , 'pop_study(''redraw'', gcf);' };
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
    if ~isempty(ALLEEG)
         options = { 'ALLEEG' ALLEEG };
    else options = {};
    end;
    if ~strcmpi(result{1}, STUDY.name ), options = { options{:} 'name'        result{1} }; end;
    if ~strcmpi(result{2}, STUDY.task ), options = { options{:} 'task'        result{2} }; end;
    if ~strcmpi(result{3}, STUDY.notes), options = { options{:} 'notes'       result{3} }; end;
    if ~isempty(result{4}), options = { options{:} 'commands'    allcom    }; end;
    if ~isempty( outstruct(1).studyfile )
        options = { options{:} 'filename' { outstruct(1).studyfile } };
    end;
    if outstruct(1).save_dataset == 1
         options = { options{:} 'savedat' 'on' };
    else options = { options{:} 'savedat' 'off' };
    end;
    if outstruct(1).copy_to_dataset == 1
         options = { options{:} 'updatedat' 'on' };
    else options = { options{:} 'updatedat' 'off' };
    end;
    
    % run command and create history
    % ------------------------------
    [STUDY ALLEEG] = editstudy(options{:});
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

        case 'clear'
            guiindex  = varargin{1};
            realindex = guiindex+(page-1)*10;
            
            datasetinfo(realindex).filename  = '';   
            datasetinfo(realindex).subject   = '';
            datasetinfo(realindex).session   = [];
            datasetinfo(realindex).condition = '';
            datasetinfo(realindex).group     = '';                    

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
            datasetinfo(realindex).filename  = fullfile(ALLEEG(realindex).filepath, ALLEEG(realindex).filename);   
            datasetinfo(realindex).subject   = ALLEEG(realindex).subject;
            datasetinfo(realindex).session   = ALLEEG(realindex).session;
            datasetinfo(realindex).condition = ALLEEG(realindex).condition;
            datasetinfo(realindex).group     = ALLEEG(realindex).group;                    

            allcom = { allcom{:} { 'index' realindex 'load' filename } };
            userdat{1} = ALLEEG;
            userdat{2} = datasetinfo;
            userdat{4} = allcom;
            set(hdl, 'userdata', userdat);
            pop_study('redraw', hdl);

        case 'redraw'
            for k = 1:10
                kk = k+(page-1)*10; % real index
                if kk > length(datasetinfo)
                    set(findobj('parent', hdl, 'tag',['set' num2str(k)]), 'string', '');
                    set(findobj('parent', hdl, 'tag',['sub' num2str(k)]), 'string','');
                    set(findobj('parent', hdl, 'tag',['sess' num2str(k)]), 'string','');
                    set(findobj('parent', hdl, 'tag',['cond' num2str(k)]), 'string','');
                    set(findobj('parent', hdl, 'tag',['group' num2str(k)]), 'string','');
                else
                    set(findobj('parent', hdl, 'tag',['set' num2str(k)]), 'string', datasetinfo(kk).filename);
                    set(findobj('parent', hdl, 'tag',['sub' num2str(k)]), 'string',datasetinfo(kk).subject);
                    set(findobj('parent', hdl, 'tag',['sess' num2str(k)]), 'string',int2str(datasetinfo(kk).session));
                    set(findobj('parent', hdl, 'tag',['cond' num2str(k)]), 'string',datasetinfo(kk).condition);
                    set(findobj('parent', hdl, 'tag',['group' num2str(k)]), 'string',datasetinfo(kk).group);
                end;
            end
            set(findobj('parent', hdl, 'tag','page'), 'string', [ 'Page ' int2str(page) ] );
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
    anycondition = all(~cellfun('isempty', { datasetinfo.condition }));
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