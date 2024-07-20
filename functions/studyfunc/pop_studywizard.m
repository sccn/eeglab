% POP_STUDYWIZARD - browse for dataset to create an EEGLAB STUDY
%
% Usage:  
%  >> [ STUDY ALLEEG ] = pop_studywizard; % create new study interactively
%  >> [ STUDY ALLEEG ] = pop_studywizard( 'key', 'val', ...); % edit study
%
% Optional Inputs:
%   All "'key', 'val'" inputs of STD_EDITSET may be used.
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
%  "dataset filename/browse" - Load dataset from specified filename. It may be 
%                          an EEGLAB dataset or any file supported by File-IO (if 
%                          it is installed from the plugin manager).
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
%  "task"                - [edit box] task for the file. {default: []}
%
% See also: POP_STUDY, EEG_IMPORT
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, July 22, 2024

% Copyright (C) Arnaud Delorme, 2024
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

function [STUDY, ALLEEG, com]  = pop_studywizard(STUDY, ALLEEG, varargin)

com = '';

if nargin == 0
    STUDY.name     = '';
    STUDY.task     = '';
    STUDY.notes    = '';
    STUDY.filename = '';
    STUDY.cluster  = [];
    STUDY.history = 'STUDY = [];';
    ALLEEG = [];

    % check events
    fieldList = { 'condition' 'group' 'subject' 'session' 'run' 'task' };
    
    % set initial datasetinfo
    % -----------------------
    datasetinfo = [];
    nextpage    = 'pop_studywizard(''nextpage'', gcbf);';
    prevpage    = 'pop_studywizard(''prevpage'', gcbf);';
    delset      = 'pop_studywizard(''clear'',    gcbf, get(gcbo, ''userdata''));';
    loadset     = 'pop_studywizard(''load'',     gcbf, get(guiind, ''userdata''), get(guiind, ''string'')); clear guiind;';
    loadsetedit = [ 'guiind = gcbo;' loadset ];
    subcom      = 'pop_studywizard(''subject''  , gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    sescom      = 'pop_studywizard(''session''  , gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    runcom      = 'pop_studywizard(''run''      , gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    condcom     = 'pop_studywizard(''condition'', gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    grpcom      = 'pop_studywizard(''group''    , gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
    taskcom     = 'pop_studywizard(''task''     , gcbf, get(gcbo, ''userdata''), get(gcbo, ''string''));';
	
	guispec = { ...
        {'style' 'text' 'string' 'Select datasets' 'FontWeight' 'Bold' 'HorizontalAlignment' 'center'} ...
        {},{}, ...
        {'style' 'text' 'string' 'Raw dataset filename' 'userdata' 'addt'} {'style' 'text' 'string' 'browse' 'userdata' 'addt'} ...
        {'style' 'text' 'string' 'subject'    'userdata' 'addt'} ...
        {'style' 'text' 'string' 'session'    'userdata' 'addt'} ...
        {'style' 'text' 'string' 'run'        'userdata' 'addt'} ...
        {'style' 'text' 'string' 'condition'  'userdata' 'addt'} ...
        {'style' 'text' 'string' 'group'      'userdata' 'addt'} ...
        {'style' 'text' 'string' 'task'       'userdata' 'addt'} ...
        {} };
    geomUnit  = [0.2 1 3.5];
    geomUnit2 = [0.2 1.05 0.35 0.4 0.35 0.25 0.6 0.4 0.6 0.3];
	guigeom = { 1 1 geomUnit2 };
	
    % create edit boxes
    % -----------------
    for index = 1:10
        guigeom = { guigeom{:} geomUnit2 };
        select_com = ['[inputname, inputpath] = uigetfile2(''*'', ''Choose dataset to add -- pop_studywizard()'');'...
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
        {'style' 'edit'       'string' ''      'tag' [ 'task' int2str(index) ]  'userdata' index 'Callback' taskcom}, ...
        {'style' 'pushbutton' 'string' 'CLear' 'tag' [ 'clear' int2str(index) ] 'userdata' index 'callback' delset} };
    end
    
    guispec = { guispec{:}, ...
                {},{}, ...
                {'style' 'pushbutton' 'string'  '<'      'Callback' prevpage 'userdata' 'addt'}, ...
                {'style' 'text'       'string'  'Page 1' 'tag' 'page' 'horizontalalignment' 'center' }, ... 
                {'style' 'pushbutton' 'string'  '>'      'Callback' nextpage 'userdata' 'addt'} {} };
    guigeom = { guigeom{:} 1 [1 0.2 0.3 0.2 1] };

    fig_arg{1} = ALLEEG;      % datasets         
    fig_arg{2} = datasetinfo; % datasetinfo
    fig_arg{3} = 1;           % page
    fig_arg{4} = {};          % all commands
    fig_arg{5} = (length(STUDY.cluster) > 1); % are cluster present
    fig_arg{6} = STUDY;       % are cluster present
    fig_arg{7} = [];          % tags
    fig_arg{8} = 'eeg';       % modality
    
    % generate GUI
    % ------------
    optiongui = { 'geometry', guigeom, ...
                  'uilist'  , guispec, ...
                  'helpcom' , 'pophelp(''pop_studywizard'')', ...
                  'title'   , 'Create a new STUDY set -- pop_studywizard()', ...
                  'userdata', fig_arg };
	[~, ~, ~, ~, instruct, alltags] = inputgui( 'mode', 'plot', optiongui{:}); % plot only to get tags

    % save tags
    currentfig = gcf;
    fig_arg{7} = alltags; 
    set(currentfig, 'userdata', fig_arg);
    
    % update GUI and wait for user input
    pop_studywizard('redraw', currentfig);
    [result, userdat2, ~, outstruct] = inputgui( 'mode', currentfig, optiongui{:});
    if isempty(result), return; end

    while test_wrong_parameters(currentfig)
    	[result, userdat2, ~, outstruct] = inputgui( 'mode', currentfig, optiongui{:});
        if isempty(result), return; end
    end
    close(currentfig);

    % convert GUI selection to options
    % --------------------------------
    allcom = simplifycom(userdat2{4});
    options = {};
    if ~isempty(allcom),                 options = { options{:} 'commands'    allcom    }; end
    options = { options{:} 'updatedat' 'on' };
    
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
    com = sprintf( '[STUDY ALLEEG] = std_editset( STUDY, [], %s );\n[STUDY ALLEEG] = std_checkset(STUDY, ALLEEG);', vararg2str(options) );
    [STUDY, ALLEEG] = std_editset(STUDY, ALLEEG, options{:});
    STUDY.etc.modality = userdat2{9};

    % check sessions
    % --------------
    allSessions = [ STUDY.datasetinfo.session ];
    if ~isempty(allSessions)
        if ~ismember(1, allSessions)
            fprintf(2, 'Sessions are usually numbered starting at 1; Do not use session to store custom information\n');
        end
    end
    
elseif nargin > 1 && isa(ALLEEG, 'matlab.ui.Figure')
    
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
    tags        = userdat{7};

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

        case 'task'
            guiindex  = varargin{1};
            realindex = guiindex+(page-1)*10;
            datasetinfo(realindex).task   = varargin{2};
            if get(findobj(hdl, 'tag', 'copy_to_dataset'), 'value')
                ALLEEG(realindex).task    = varargin{2};
            end
            allcom = { allcom{:}, { 'index' realindex 'task' varargin{2} } };
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
            datasetinfo(realindex).task      = '';                    
            datasetinfo(realindex).comps     = [];                    

            allcom = { allcom{:}, { 'remove' realindex } };
            userdat{1} = ALLEEG;
            userdat{2} = datasetinfo;
            userdat{4} = allcom;
            set(hdl, 'userdata', userdat);
            pop_studywizard('redraw', hdl);

        case 'nextpage'
            userdat{3} = page+1;
            set(hdl, 'userdata', userdat);
            pop_studywizard('redraw', hdl);

        case 'prevpage'
            userdat{3} = max(1,page-1);
            set(hdl, 'userdata', userdat);
            pop_studywizard('redraw', hdl);

        case 'load'
            guiindex  = varargin{1};
            filename  = varargin{2};
            realindex = guiindex+(page-1)*10;

            % load dataset
            % ------------
            [pathTmp, fileTmp, ext] = fileparts(filename);
            if isequal(ext, '.set')
                modality = 'eeg';
                TMPEEG = pop_loadset('filename', filename, 'loadmode', 'info');
            else
                fprintf('Converting file to EEGLAB format...\n')
                [TMPEEG, modality] = eeg_import(filename);
                filename = fullfile(pathTmp, [ fileTmp '_bids_tmp_.set' ]);
                pop_saveset(TMPEEG, 'filename', filename);
                TMPEEG = pop_loadset('filename', filename, 'loadmode', 'info');
            end
            ALLEEG = eeg_store(ALLEEG, eeg_checkset(TMPEEG), realindex);
            
            % update datasetinfo structure
            % ----------------------------
            if ~isfield(ALLEEG, 'task'), ALLEEG(realindex).task = ''; end
            datasetinfo(realindex).filename  = ALLEEG(realindex).filename;   
            datasetinfo(realindex).filepath  = ALLEEG(realindex).filepath;   
            datasetinfo(realindex).subject   = ALLEEG(realindex).subject;
            datasetinfo(realindex).session   = ALLEEG(realindex).session;
            datasetinfo(realindex).run       = ALLEEG(realindex).run;
            datasetinfo(realindex).condition = ALLEEG(realindex).condition;
            datasetinfo(realindex).task      = ALLEEG(realindex).task;
            datasetinfo(realindex).group     = ALLEEG(realindex).group;                    
            datasetinfo(realindex).comps     = [];                    

            allcom = { allcom{:}, { 'index' realindex 'load' filename } };
            userdat{1} = ALLEEG;
            userdat{2} = datasetinfo;
            userdat{4} = allcom;
            userdat{9} = modality;
            set(hdl, 'userdata', userdat);
            pop_studywizard('redraw', hdl);
         
        case 'redraw'
            for k = 1:10
                kk = k+(page-1)*10; % real index
                if kk > length(datasetinfo)
                    set(tags.(['num' num2str(k)]), 'string', int2str(kk));
                    set(tags.(['set' num2str(k)]), 'string', '');
                    set(tags.(['sub' num2str(k)]), 'string','');
                    set(tags.(['sess' num2str(k)]), 'string','');
                    set(tags.(['run'  num2str(k)]), 'string','');
                    set(tags.(['cond' num2str(k)]), 'string','');
                    set(tags.(['task' num2str(k)]), 'string','');
                    set(tags.(['group' num2str(k)]), 'string','');
                else
                    set(tags.(['num' num2str(k)]), 'string', int2str(kk));
                    set(tags.(['set' num2str(k)]), 'string', fullfile(char(datasetinfo(kk).filepath), char(datasetinfo(kk).filename)));
                    set(tags.(['sub' num2str(k)]), 'string', char(datasetinfo(kk).subject));
                    set(tags.(['sess' num2str(k)]), 'string', int2str(datasetinfo(kk).session));
                    set(tags.(['run' num2str(k)]),  'string', int2str(datasetinfo(kk).run));
                    set(tags.(['cond' num2str(k)]), 'string', char(datasetinfo(kk).condition));
                    set(tags.(['task' num2str(k)]), 'string', char(datasetinfo(kk).task));
                    set(tags.(['group' num2str(k)]), 'string', char(datasetinfo(kk).group));
                end
            end
            if page<10
                 pagestr =  [ ' Page ' int2str(page) ];
            else pagestr =  [ 'Page ' int2str(page) ];
            end
            set(tags.('page'), 'string', pagestr );
    end
else
    % script execution
    [STUDY, ALLEEG] = std_editset(STUDY, ALLEEG, varargin{:});
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
    anytask      = any(~cellfun('isempty', { datasetinfo(nonempty).task }));
    alltask      = all(~cellfun('isempty', { datasetinfo(nonempty).task }));
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
    if anytask && ~alltask
         bool = 1; warndlg2('If one dataset has a task, they must all have one', 'Error');
    end
    if anydipfit && ~alldipfit
         bool = 1; warndlg2('Dipole''s data across datasets is not uniform');
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
