% pop_study()  - create a new STUDY set structure defining a group of related EEG datasets. 
%                The STUDY set also contains information about each of the datasets: the 
%                subject code, subject group, experimental condition, and session. This can 
%                be provided interactively in a pop-up window or be automatically filled 
%                in by the function. Defaults: Assume a different subject for each 
%                dataset and only one condition; leave subject group and session fields 
%                empty. Additional STUDY information about the STUDY name, task and 
%                miscellaneous notes can also be saved in the STUDY structure.
% Usage:  
%  >> [ STUDY ALLEEG ] = pop_study([],[], 'gui', 'on'); % create new study interactively
%  >> [ STUDY ALLEEG ] = pop_study(STUDY, ALLEEG, 'gui', 'on'); % edit study interactively
%  >> [ STUDY ALLEEG ] = pop_study(STUDY, ALLEEG, 'key', 'val', ...); % edit study
%
% Optional Inputs:
%   STUDY                - existing study structure. 
%   ALLEEG               - vector of EEG dataset structures to be included in the STUDY. 
%
% Optional Inputs:
%   All "'key', 'val'" inputs of std_editset() may be used.
%
% Outputs:
%   STUDY                - new STUDY set comprising some or all of the datasets in 
%                          ALLEEG, plus other information about the experiments. 
%   ALLEEG               - an updated ALLEEG structure including the STUDY datasets. 
%
% Graphic interface buttons:
%  "STUDY set name"      - [edit box] name for the STUDY structure {default: ''}   
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
%  "run"                 - [edit box] dataset run. If no run information is
%                          provided, all datasets that belong to one subject are assumed to
%                          have been recorded within one run {default: []}
%  "condition"           - [edit box] dataset condition. If no condition code is provided,
%                          all datasets are assumed to be from the same condition {default:[]}
%  "group"               - [edit box] the subject group the dataset belongs to. If no group
%                          is provided, all subjects and datasets are assumed to belong to
%                          the same group. {default: []}
%  "Save this STUDY set to disk file" - [check box] If checked, save the new STUDY set
%                          structure to disk. If no filename is provided, a window will 
%                          pop up to ask for it. 
%
% See also: std_editset, pop_loadstudy(), pop_preclust(), pop_clust()
%
% Authors: Arnaud Delorme, Hilit Serby, Scott Makeig, SCCN, INC, UCSD, July 22, 2005

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, July 22, 2005, arno@sccn.ucsd.edu
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
if ~ischar(STUDY) %initial settings
    mode = 'script';
    if nargin > 2
        for index = 1:length(varargin)
            if ischar(varargin{index})
                if strcmpi(varargin{index}, 'gui')
                    mode = 'gui';
                end
            end
        end
    end
end

if isempty(STUDY)
    newstudy       = 1;
    STUDY.name     = '';
    STUDY.task     = '';
    STUDY.notes    = '';
    STUDY.filename = '';
    STUDY.cluster  = [];
    STUDY.history = 'STUDY = [];';
else
    newstudy       = 0;
end

if strcmpi(mode, 'script') % script mode
    [STUDY, ALLEEG] = std_editset(STUDY, ALLEEG, varargin{:});
    return;
elseif strcmpi(mode, 'gui') % GUI mode

    % check events
    fieldList = { 'condition' 'group' 'subject' 'session' 'run' };
    if ~isempty(ALLEEG) && isstruct(ALLEEG) && isfield(ALLEEG, 'event')
        for iField = 1:length(fieldList)
            if any(cellfun(@(x)isfield(x, fieldList{iField}), {ALLEEG.event}))
                fprintf(2, 'Rename field "%s" in datasets'' event structure as it conflicts with the STUDY field bearing the same name\n', fieldList{iField});
            end
        end
    end
    
    % show warning if necessary
    % -------------------------
    if isreal(ALLEEG)
        if ALLEEG == 0
            res = questdlg2( strvcat('Datasets currently loaded will be removed from EEGLAB memory.', ...
                                     'Are you sure you want to continue?'), ...
                                     'Discard loaded EEGLAB datasets?', 'Cancel', 'Yes', 'Yes');
            if strcmpi(res, 'cancel'), return; end
        end
        ALLEEG = [];
        alleegEmpty = 1;
    else
        alleegEmpty = 0;
    end
    
    % set initial datasetinfo
    % -----------------------
    if isfield(STUDY, 'datasetinfo')
        datasetinfo = STUDY.datasetinfo;
        different = 0;
        for k = 1:length(ALLEEG)
            if ~strcmpi(datasetinfo(k).filename, ALLEEG(k).filename), different = 1; break; end
            if ~strcmpi(datasetinfo(k).subject,   ALLEEG(k).subject),   different = 1; break; end
            if ~strcmpi(datasetinfo(k).condition, ALLEEG(k).condition), different = 1; break; end
            if ~strcmpi(char(datasetinfo(k).group), char(ALLEEG(k).group)),     different = 1; break; end       
            if ~isequal(datasetinfo(k).session, ALLEEG(k).session),             different = 1; break; end
            if ~isequal(datasetinfo(k).run, ALLEEG(k).run),                     different = 1; break; end
        end
        if different
            info = 'from_STUDY_different_from_ALLEEG';
        else
            info = 'from_STUDY';
        end
        if ~isfield(datasetinfo, 'comps');
            datasetinfo(1).comps = [];
        end
    else
        info = 'from_ALLEEG';
        if length(ALLEEG) > 0
            datasetinfo(length(ALLEEG)).filename  = [];
            datasetinfo(length(ALLEEG)).filepath  = [];
            datasetinfo(length(ALLEEG)).subject   = [];
            datasetinfo(length(ALLEEG)).session   = [];
            datasetinfo(length(ALLEEG)).run       = [];
            datasetinfo(length(ALLEEG)).condition = [];
            datasetinfo(length(ALLEEG)).group     = [];     
            for k = 1:length(ALLEEG)
                datasetinfo(k).filename  = ALLEEG(k).filename;   
                datasetinfo(k).filepath  = ALLEEG(k).filepath;   
                datasetinfo(k).subject   = ALLEEG(k).subject;
                datasetinfo(k).session   = ALLEEG(k).session;
                datasetinfo(k).run       = ALLEEG(k).run;
                datasetinfo(k).condition = ALLEEG(k).condition;
                datasetinfo(k).group     = ALLEEG(k).group;                    
            end
            if ~isfield(datasetinfo, 'comps');
                datasetinfo(1).comps = [];
            end
        else
            datasetinfo = [];
        end
    end

    nextpage    = 'pop_study(''nextpage'', gcbf);';
    prevpage    = 'pop_study(''prevpage'', gcbf);';
    delset      = 'pop_study(''clear'',    gcbf, get(gcbo, ''userdata''));';
    loadset     = 'pop_study(''load'',     gcbf, get(guiind, ''userdata''), get(guiind, ''string'')); clear guiind;';
    loadsetedit = [ 'guiind = gcbo;' loadset ];
    subcom      = 'pop_study(''subject''  , gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    sescom      = 'pop_study(''session''  , gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    runcom      = 'pop_study(''run''      , gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    condcom     = 'pop_study(''condition'', gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    grpcom      = 'pop_study(''group''    , gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    compcom     = 'pop_study(''component'', gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    cb_del      = 'pop_study(''delclust'' , gcbf, ''showwarning'');';
    cb_dipole   = 'pop_study(''dipselect'', gcbf, ''showwarning'');';

	browsestudy = [ '[filename, filepath] = uiputfile2(''*.study'', ''Use existing STUDY set to import dataset information -- pop_study()''); ' ... 
                      'set(findobj(''parent'', gcbf, ''tag'', ''usestudy_file''), ''string'', [filepath filename]);' ];
	saveSTUDY = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''save''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
	browsesave = [ '[filename, filepath] = uiputfile2(''*.study'', ''Save STUDY with .study extension -- pop_clust()''); ' ... 
                   'if filename ~= 0,' ...
                   '   set(findobj(''parent'', gcbf, ''tag'', ''studyfile''), ''string'', [filepath filename]);' ...
                   'end;' ...
                   'clear filename filepath;' ];
	
    texthead = fastif(newstudy, 'Create a new STUDY set', 'Edit STUDY set information - remember to save changes');
	guispec = { ...
        {'style' 'text' 'string' texthead 'FontWeight' 'Bold' 'HorizontalAlignment' 'center'} ...
        {} {'style' 'text' 'string' 'STUDY set name:' } { 'style' 'edit' 'string' char(STUDY.name) 'tag' 'study_name' } ...
        {} {'style' 'text' 'string' 'STUDY set task name:' } { 'style' 'edit' 'string' char(STUDY.task) 'tag' 'study_task' } ...
        {} {'style' 'text' 'string' 'STUDY set notes:' } { 'style' 'edit' 'string' char(STUDY.notes) 'tag' 'study_notes' } {}...
        {} ...
        {'style' 'text' 'string' 'dataset filename' 'userdata' 'addt'} {'style' 'text' 'string' 'browse' 'userdata' 'addt'} ...
        {'style' 'text' 'string' 'subject'    'userdata' 'addt'} ...
        {'style' 'text' 'string' 'session'    'userdata' 'addt'} ...
        {'style' 'text' 'string' 'run'        'userdata' 'addt'} ...
        {'style' 'text' 'string' 'condition'  'userdata' 'addt'} ...
        {'style' 'text' 'string' 'group'      'userdata' 'addt'} ...
        {'style' 'pushbutton' 'string' 'Select by r.v.' 'userdata' 'addt' 'callback' cb_dipole } ...
        {} };
    geomUnit  = [0.2 1 3.5];
    geomUnit2 = [0.2 1.05 0.35 0.4 0.35 0.25 0.6 0.4 0.6 0.3];
	guigeom = { [1] geomUnit geomUnit geomUnit [1] geomUnit2 };
	
    % create edit boxes
    % -----------------
    for index = 1:10
        guigeom = { guigeom{:} geomUnit2 };
        select_com = ['[inputname, inputpath] = uigetfile2(''*.set;*.SET'', ''Choose dataset to add to STUDY -- pop_study()'');'...
                      'if inputname ~= 0,' ...
                      '   guiind = findobj(''parent'', gcbf, ''tag'', ''set' int2str(index) ''');' ...
                      '   set( guiind,''string'', [inputpath inputname]);' ...
                          loadset ...
                      'end; clear inputname inputpath;'];
        numstr = int2str(index);
        guispec = { guispec{:}, ...
        {'style' 'text'       'string' numstr  'tag' [ 'num'   int2str(index) ] 'userdata' index }, ...
        {'style' 'edit'       'string' ''      'tag' [ 'set'   int2str(index) ] 'userdata' index 'callback' loadsetedit}, ...
        {'style' 'pushbutton' 'string' '...'   'tag' [ 'brw'   int2str(index) ] 'userdata' index 'Callback' select_com}, ...
        {'style' 'edit'       'string' ''      'tag' [ 'sub'   int2str(index) ] 'userdata' index 'Callback' subcom}, ...
        {'style' 'edit'       'string' ''      'tag' [ 'sess'  int2str(index) ] 'userdata' index 'Callback' sescom}, ...
        {'style' 'edit'       'string' ''      'tag' [ 'run'   int2str(index) ] 'userdata' index 'Callback' runcom}, ...
        {'style' 'edit'       'string' ''      'tag' [ 'cond'  int2str(index) ] 'userdata' index 'Callback' condcom}, ...
        {'style' 'edit'       'string' ''      'tag' [ 'group' int2str(index) ] 'userdata' index 'Callback' grpcom}, ...
        {'style' 'pushbutton' 'string' 'All comp.'   'tag' [ 'comps' int2str(index) ] 'userdata' index 'Callback' compcom}, ...
        {'style' 'pushbutton' 'string' 'CLear' 'tag' [ 'clear' int2str(index) ] 'userdata' index 'callback' delset} };
    end
    
    if strcmpi(info, 'from_STUDY_different_from_ALLEEG')
        text1    = 'Dataset info (condition, group, ...) differs from study info. [set] = Overwrite dataset info for each dataset on disk.';
        value_cb = 0;
    else
        text1    = 'Update dataset info - datasets stored on disk will be overwritten (unset = Keep study info separate).';
        value_cb = 1;
    end
    guispec = { guispec{:}, ...
                {'style' 'text'       'string'  'Important note: Removed datasets will not be saved before being deleted from EEGLAB memory' }, ...
                {}, ...
                {'style' 'pushbutton' 'string'  '<'      'Callback' prevpage 'userdata' 'addt'}, ...
                {'style' 'text'       'string'  'Page 1' 'tag' 'page' 'horizontalalignment' 'center' }, ... 
                {'style' 'pushbutton' 'string'  '>'      'Callback' nextpage 'userdata' 'addt'}, {}, ...
                {}, ...
                {'style' 'checkbox'   'value'   value_cb 'tag' 'copy_to_dataset' }, ...
                {'style' 'text'       'string'  text1 }, ...
                {'style' 'checkbox'   'value'   0        'tag' 'delclust' 'callback' cb_del }, ...
                {'style' 'text'       'string'  'Delete cluster information (to allow loading new datasets, set new components for clustering, etc.)' } };
	guigeom = { guigeom{:} [1] [1 0.2 0.3 0.2 1] [1] [0.14 3] [0.14 3] };

%     if ~isempty(STUDY.filename)
%         guispec{end-3} = {'style' 'checkbox' 'string' '' 'value' 0 'tag' 'studyfile' };
%         guispec{end-2} = {'style' 'text'     'string' 'Re-save STUDY. Uncheck and use menu File > Save study as to save under a new filename'};
%         guispec(end-1) = [];
%         guigeom{end-1} = [0.14 3];
%     end
	
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
	[result, userdat2, strhalt, outstruct, instruct] = inputgui( 'mode', 'noclose', optiongui{:});
    if isempty(result), return; end
    if ~isempty(get(0, 'currentfigure')) currentfig = gcf; end
    
    while test_wrong_parameters(currentfig)
    	[result, userdat2, strhalt, outstruct] = inputgui( 'mode', currentfig, optiongui{:});
        if isempty(result), return; end
    end
    close(currentfig);

    % convert GUI selection to options
    % --------------------------------
    allcom = simplifycom(userdat2{4});
    options = {};
    if ~strcmpi(result{1}, STUDY.name ), options = { options{:} 'name'        result{1} }; end
    if ~strcmpi(result{2}, STUDY.task ), options = { options{:} 'task'        result{2} }; end
    if ~strcmpi(result{3}, STUDY.notes), options = { options{:} 'notes'       result{3} }; end
    if ~isempty(allcom),                 options = { options{:} 'commands'    allcom    }; end
%     if isnumeric(outstruct(1).studyfile)
%         if outstruct(1).studyfile == 1,  options = { options{:} 'resave'      'on' }; end
%     else
%         if ~isempty(outstruct(1).studyfile), options = { options{:} 'filename' outstruct(1).studyfile }; end
%     end
    if outstruct(1).copy_to_dataset == 1
         options = { options{:} 'updatedat' 'on' };
         eeglab_options;
         if option_storedisk
             options = { options{:} 'savedat' 'on' };
         end
    else options = { options{:} 'updatedat' 'off' };
    end
    
    if outstruct(1).delclust == 1
        options = { options{:} 'rmclust' 'on' };
    else
        options = { options{:} 'rmclust' 'off' };
    end
    
    % ---
    if ~isequal(outstruct, instruct) && (outstruct(1).delclust ~= 1) % notice that isequal is sensitive to fields order. isequaln isn't backward compatible 
        strfields = fieldnames(outstruct);
        for i = 1:length(strfields)
            strdiff(i) = strcmp(outstruct.(strfields{i}),instruct.(strfields{i}));
        end
        % If the information of the STUDY differ, then remove information from clusters
        % Ignoring ('STUDY set name','STUDY set task','STUDY set notes')
        if any(~strdiff(3:end-2)) 
            options{find(strcmp(options,'rmclust'))+1} = 'on';
        end
    end
    % ---
    
    % check channel labels
    % --------------------
    TMPALLEEG   = userdat2{1};
    if isfield(TMPALLEEG, 'chanlocs')
        allchans = { TMPALLEEG.chanlocs };
        if any(cellfun('isempty', allchans))
            txt = strvcat('Some datasets do not have channel labels. Do you wish to generate', ...
                          'channel labels automatically for all datasets ("1" for channel 1,', ...
                          '"2" for channel 2, ...). Datasets will be overwritten on disk.', ...
                          'If you abort, the STUDY will not be created.');
            res = questdlg2(txt, 'Dataset format problem', 'Yes', 'No, abort', 'Yes');
            if strcmpi(res, 'yes'), options = { options{:} 'addchannellabels' 'on' 'savedat' 'on'}; 
            else return;
            end
        end
    end
        
    % run command and create history
    % ------------------------------
    if alleegEmpty
        com = sprintf( '[STUDY ALLEEG] = std_editset( STUDY, [], %s );\n[STUDY ALLEEG] = std_checkset(STUDY, ALLEEG);', vararg2str(options) );
    else
        com = sprintf( '[STUDY ALLEEG] = std_editset( STUDY, ALLEEG, %s );\n[STUDY ALLEEG] = std_checkset(STUDY, ALLEEG);', vararg2str(options) );
    end
    [STUDY, ALLEEG] = std_editset(STUDY, ALLEEG, options{:});

    % check sessions
    % --------------
    allSessions = [ STUDY.datasetinfo.session ];
    if ~isempty(allSessions)
        if ~ismember(1, allSessions)
            fprintf(2, 'Sessions are usually numbered starting at 1; Do not use session to store custom information\n');
        end
    end
    
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
            end
            allcom = { allcom{:}, { 'index' realindex 'subject' varargin{2} } };
            userdat{1} = ALLEEG;
            userdat{2} = datasetinfo;
            userdat{4} = allcom;
            set(hdl, 'userdata', userdat);            

        case 'run'
            guiindex  = varargin{1};
            realindex = guiindex+(page-1)*10;
            datasetinfo(realindex).run = str2num(varargin{2});
            if get(findobj(hdl, 'tag', 'copy_to_dataset'), 'value')
                ALLEEG(realindex).run  = str2num(varargin{2});
            end
            allcom = { allcom{:}, { 'index' realindex 'run' str2num(varargin{2}) } };
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
            end
            allcom = { allcom{:}, { 'index' realindex 'session' str2num(varargin{2}) } };
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
            end
            allcom = { allcom{:}, { 'index' realindex 'group' varargin{2} } };
            userdat{1} = ALLEEG;
            userdat{2} = datasetinfo;
            userdat{4} = allcom;
            set(hdl, 'userdata', userdat);            

        case 'condition'
            guiindex  = varargin{1};
            realindex = guiindex+(page-1)*10;
            datasetinfo(realindex).condition   = varargin{2};
            if get(findobj(hdl, 'tag', 'copy_to_dataset'), 'value')
                ALLEEG(realindex).condition     = varargin{2};
            end
            allcom = { allcom{:}, { 'index' realindex 'condition' varargin{2} } };
            userdat{1} = ALLEEG;
            userdat{2} = datasetinfo;
            userdat{4} = allcom;
            set(hdl, 'userdata', userdat);            
            
        case 'dipselect'
            STUDY.datasetinfo = datasetinfo;
            
            res = inputdlg2_with_checkbox( { strvcat('Enter maximum residual (topo map - dipole proj.) var. (in %)', ...
                                       'NOTE: This will delete any existing component clusters!') }, ...
                                       'pop_study():  Pre-select components', 1, { '15' },'pop_study' );
           
            if isempty(res), return; end
            if res{2} == 1
                STUDY = std_editset(STUDY, ALLEEG, 'commands', { 'inbrain', 'on', 'dipselect' str2num(res{1})/100 'return' });
                allcom = { allcom{:}, { 'inbrain', 'on', 'dipselect' str2num(res{1})/100 } };
            else
                STUDY = std_editset(STUDY, ALLEEG, 'commands', { 'inbrain', 'off','dipselect' str2num(res{1})/100 'return' });
                allcom = { allcom{:}, { 'inbrain', 'off', 'dipselect' str2num(res{1})/100 } };
            end
                
            datasetinfo   = STUDY.datasetinfo;
            
            userdat{2} = datasetinfo;
            userdat{4} = allcom;
            set(hdl, 'userdata', userdat);   
            set(findobj(hdl, 'tag', 'delclust'), 'value', 1);
            pop_study('delclust', hdl);
            pop_study('redraw', hdl);

        case 'component'
            guiindex  = varargin{1};
            realindex = guiindex+(page-1)*10;
            
            for index = 1:size(ALLEEG(realindex).icaweights,1)
                complist{index} = [ 'IC ' int2str(index) ];
            end
            [tmps,tmpv] = listdlg2('PromptString', 'Select components', 'SelectionMode', ...
                                    'multiple', 'ListString', strvcat(complist), 'initialvalue', datasetinfo(realindex).comps);
            if tmpv ~= 0 % no cancel                
                
                % find other subjects with the same session
                % -----------------------------------------
                for index = 1:length(datasetinfo)
                    if realindex == index || (strcmpi(datasetinfo(index).subject, datasetinfo(realindex).subject) && ...
                                ~isempty(datasetinfo(index).subject) && ...
                                isequal( datasetinfo(index).session, datasetinfo(realindex).session ) && ...
                                isequal( datasetinfo(index).run, datasetinfo(realindex).run ) )
                        datasetinfo(index).comps = tmps;
                        allcom = { allcom{:}, { 'index' index 'comps' tmps } };
                        set(findobj('tag', [ 'comps' int2str(index) ]), ...
                            'string', formatbut(tmps), 'horizontalalignment', 'left');
                    end
                end
                
            end
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
            datasetinfo(realindex).run       = [];
            datasetinfo(realindex).condition = '';
            datasetinfo(realindex).group     = '';                    
            datasetinfo(realindex).comps     = [];                    

            allcom = { allcom{:}, { 'remove' realindex } };
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
            datasetinfo(realindex).run       = ALLEEG(realindex).run;
            datasetinfo(realindex).condition = ALLEEG(realindex).condition;
            datasetinfo(realindex).group     = ALLEEG(realindex).group;                    
            datasetinfo(realindex).comps     = [];                    

            allcom = { allcom{:}, { 'index' realindex 'load' filename } };
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
                    end
                else
                    for k = 1:10
                        set(findobj('parent', hdl, 'tag',['set'   num2str(k)]), 'style', 'edit');
                        set(findobj('parent', hdl, 'tag',['comps' num2str(k)]), 'enable', 'on');
                        set(findobj('parent', hdl, 'tag',['sess'  num2str(k)]), 'enable', 'on');
                        set(findobj('parent', hdl, 'tag',['brw'   num2str(k)]), 'enable', 'on');
                    end
                end
            else 
                set(findobj(hdl, 'tag', 'delclust'), 'value', 0)
                if nargin > 2
                    warndlg2('No cluster present');
                end
            end
         
        case 'redraw'
            for k = 1:10
                kk = k+(page-1)*10; % real index
                if kk > length(datasetinfo)
                    set(findobj('parent', hdl, 'tag',['num' num2str(k)]), 'string', int2str(kk));
                    set(findobj('parent', hdl, 'tag',['set' num2str(k)]), 'string', '');
                    set(findobj('parent', hdl, 'tag',['sub' num2str(k)]), 'string','');
                    set(findobj('parent', hdl, 'tag',['sess' num2str(k)]), 'string','');
                    set(findobj('parent', hdl, 'tag',['run'  num2str(k)]), 'string','');
                    set(findobj('parent', hdl, 'tag',['cond' num2str(k)]), 'string','');
                    set(findobj('parent', hdl, 'tag',['comps' num2str(k)]), 'string','');
                    set(findobj('parent', hdl, 'tag',['group' num2str(k)]), 'string','');
                else
                    set(findobj('parent', hdl, 'tag',['num' num2str(k)]), 'string', int2str(kk));
                    set(findobj('parent', hdl, 'tag',['set' num2str(k)]), 'string', fullfile(char(datasetinfo(kk).filepath), char(datasetinfo(kk).filename)));
                    set(findobj('parent', hdl, 'tag',['sub' num2str(k)]), 'string', char(datasetinfo(kk).subject));
                    set(findobj('parent', hdl, 'tag',['sess' num2str(k)]), 'string', int2str(datasetinfo(kk).session));
                    set(findobj('parent', hdl, 'tag',['run' num2str(k)]),  'string', int2str(datasetinfo(kk).run));
                    set(findobj('parent', hdl, 'tag',['cond' num2str(k)]), 'string', char(datasetinfo(kk).condition));
                    set(findobj('parent', hdl, 'tag',['comps' num2str(k)]), 'string', formatbut(datasetinfo(kk).comps));
                    set(findobj('parent', hdl, 'tag',['group' num2str(k)]), 'string', char(datasetinfo(kk).group));
                end
            end
            if page<10
                 pagestr =  [ ' Page ' int2str(page) ];
            else pagestr =  [ 'Page ' int2str(page) ];
            end
            set(findobj('parent', hdl, 'tag','page'), 'string', pagestr );
    end
end

% remove empty elements in allcom
% -------------------------------
function allcom = simplifycom(allcom);

    for index = length(allcom)-1:-1:1
        if strcmpi(allcom{index}{1}, 'index') && strcmpi(allcom{index+1}{1}, 'index')
            if allcom{index}{2} == allcom{index+1}{2} % same dataset index
                allcom{index}(end+1:end+length(allcom{index+1})-2) = allcom{index+1}(3:end);
                allcom(index+1) = [];
            end
        end
    end
    
% test for wrong parameters
% -------------------------
function bool = test_wrong_parameters(hdl)
    userdat     = get(hdl, 'userdata');    
    datasetinfo = userdat{2};
    datastrinfo = userdat{1};

    bool = 0;
    for index = 1:length(datasetinfo)
        if ~isempty(datasetinfo(index).filename)
            if isempty(datasetinfo(index).subject) && bool == 0
                bool = 1; warndlg2('All datasets must have a subject name or code', 'Error');
            end
        end
    end
    
    nonempty     = cellfun('isempty', { datasetinfo.filename });
    anysession   = any(~cellfun('isempty', { datasetinfo(nonempty).session }));
    allsession   = all(~cellfun('isempty', { datasetinfo(nonempty).session }));
    anyrun       = any(~cellfun('isempty', { datastrinfo(nonempty).run}));
    allrun       = all(~cellfun('isempty', { datastrinfo(nonempty).run}));
    anycondition = any(~cellfun('isempty', { datasetinfo(nonempty).condition }));
    allcondition = all(~cellfun('isempty', { datasetinfo(nonempty).condition }));
    anygroup     = any(~cellfun('isempty', { datasetinfo(nonempty).group }));
    allgroup     = all(~cellfun('isempty', { datasetinfo(nonempty).group }));
    anydipfit    = any(~cellfun('isempty', { datastrinfo(nonempty).dipfit}));
    alldipfit    = all(~cellfun('isempty', { datastrinfo(nonempty).dipfit}));

    if anygroup && ~allgroup
         bool = 1; warndlg2('If one dataset has a group label, they must all have one', 'Error');
    end
    if anycondition && ~allcondition
         bool = 1; warndlg2('If one dataset has a condition label, they must all have one', 'Error');
    end
    if anysession && ~allsession
         bool = 1; warndlg2('If one dataset has a session index, they must all have one', 'Error');
    end
    if anyrun && ~allrun
         bool = 1; warndlg2('If one dataset has a session index, they must all have one', 'Error');
    end
    if anydipfit && ~alldipfit
         bool = 1; warndlg2('Dipole''s data across datasets is not uniform');
    end
function strbut = formatbut(complist)
    if isempty(complist)
        strbut = 'All comp.';
    else
        if length(complist) > 3, strbut = [ 'Comp.: ' int2str(complist(1:2)) ' ...' ];
        else                     strbut = [ 'Comp.: ' int2str(complist) ];
        end
    end

    
%---------------------- helper functions -------------------------------------    
    
function [result] = inputdlg2_with_checkbox(Prompt,Title,LineNo,DefAns,funcname);

if nargin < 4
   help inputdlg2;
   return;
end;	
if nargin < 5
	funcname = '';
end
	
if length(Prompt) ~= length(DefAns)
	error('inputdlg2: prompt and default answer cell array must have the same size');
end

geometry = {};
listgui = {};

% determine if vertical or horizontal
% -----------------------------------
geomvert = [];
for index = 1:length(Prompt)
	geomvert = [geomvert size(Prompt{index},1) 1];  % default is vertical geometry
end
if all(geomvert == 1) && length(Prompt) > 1
	geomvert = []; % horizontal
end

for index = 1:length(Prompt)
	if ~isempty(geomvert) % vertical
		geometry = { geometry{:} [ 1] [1 ]};
	else
		geometry = { geometry{:} [ 1 0.6 ]};
	end
	listgui = { listgui{:} { 'Style', 'text', 'string', Prompt{index}}  ...
				{ 'Style', 'edit', 'string', DefAns{index} } { 'Style', 'checkbox', 'string','Keep only in-brain dipoles (requires Fieldtrip extension).','value',1 }  };
end
geometry = [1 1 1];geomvert = [2 1 1];
result = inputgui(geometry, listgui, ['pophelp(''' funcname ''');'], Title, [], 'normal', geomvert);
