% pop_studydesign() - create a STUDY design structure.
%
% Usage: 
%   >> [STUDY, ALLEEG] = pop_studydesign(STUDY, ALLEEG, key1, val1, ...);  
%
% Inputs:
%   STUDY      - EEGLAB STUDY set
%   ALLEEG     - vector of the EEG datasets included in the STUDY structure 
%
% Optional inputs:
%
% Authors: Arnaud Delorme, April 2010

% Copyright (C) Arnaud Delorme & Scott Makeig, SCCN/INC/UCSD, October 11, 2004, smakeig@ucsd.edu
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

function [STUDY allcom] = pop_studydesign(STUDY, ALLEEG, designind, varargin);  
    
allcom = '';
if nargin < 2
    help pop_studydesign;
    return;
end

if ~ischar(STUDY) && isfield(STUDY, 'currentdesign') && STUDY.currentdesign > length(STUDY.design)
    STUDY.currentdesign = length(STUDY.design);
end

if nargin < 3 && ~ischar(STUDY)
    
    %% create GUI
    [ usrdat.factors usrdat.factorvals usrdat.factsubj usrdat.pairing] = std_getindvar(STUDY, 'both', 1);
    
    usrdat.factors     = { 'None' usrdat.factors{:} };
    usrdat.factorvals  = { {}     usrdat.factorvals{:} };
    usrdat.factsubj    = { {}     usrdat.factsubj{:} };
    usrdat.subjects    = STUDY.subject;
    usrdat.datasetinfo = STUDY.datasetinfo;
    usrdat.design      = STUDY.design;
    usrdat.filepath    = STUDY.filepath;
    for ind = 1:length(usrdat.design)
        for iVar = length(usrdat.design(ind).variable):-1:1
            if isempty(usrdat.design(ind).variable(iVar).label)
                usrdat.design(ind).variable(iVar) = [];
            end
        end
    end
    
    % check numerical variables
    for ind1 = 1:length(usrdat.factors)
        if all(cellfun(@isnumeric, usrdat.factorvals{ind1}))
             usrdat.numerical(ind1) = 1;
        else usrdat.numerical(ind1) = 0;
        end
    end
        
    cb_selectsubj   = 'pop_studydesign(''selectsubj'', gcbf);';
    cb_setsubj      = 'pop_studydesign(''setsubj'', gcbf);';
    cb_rename       = 'pop_studydesign(''rename'', gcbf);';
    cb_add          = 'pop_studydesign(''add'', gcbf);';
    cb_del          = 'pop_studydesign(''del'', gcbf);';
    cb_listboxfact1 = 'pop_studydesign(''selectfact'', gcf, 0);';
    cb_listboxfact2 = 'pop_studydesign(''selectfact'', gcf, 1);';
    cb_lbval        = 'pop_studydesign(''updategui'', gcbf);';
    cb_selectdesign = 'pop_studydesign(''selectdesign'', gcbf);';
    cb_selectdata   = 'pop_studydesign(''selectdatatrials'', gcbf);';
    cb_selectfolder = 'pop_studydesign(''selectfolder'', gcbf);';
    cb_setfolder    = 'pop_studydesign(''updategui'', gcbf);';
    cb_newvar       = 'pop_studydesign(''newvar'', gcbf);';
    cb_editvar      = 'pop_studydesign(''editvar'', gcbf);';
    cb_delvar       = 'pop_studydesign(''delvar'', gcbf);';
    cb_list         = 'pop_studydesign(''list'', gcbf);';
    cb_plotdmat     = 'pop_studydesign(''plotdmat'', gcbf);';
    cb_importgvar   = 'pop_studydesign(''importgvar'', gcbf);';
    
    icadefs;
%               { 'style' 'text'       'string' 'Subjects' 'fontweight' 'bold' } ...
%               { 'style' 'listbox'    'string' usrdat.subjects 'tag' 'lbsubj' 'min' 0 'max' 2 'value' 1 'callback' cb_selectsubj } ...
    uilist = { ...
               { 'style' 'text'       'string' 'Include these subjects (default: all)' 'fontweight' 'bold' } ...
               { 'style' 'edit'       'string' '' 'tag' 'subjects' 'callback' cb_setsubj } ...
               { 'style' 'pushbutton' 'string' '...' 'callback' cb_selectsubj } ...
               { 'style' 'text'       'string' 'Design name' 'fontweight' 'bold' } ...
               { 'style' 'listbox'    'string' { usrdat.design.name } 'tag' 'listboxdesign' 'callback' cb_selectdesign 'value' STUDY.currentdesign } ...
               { 'style' 'pushbutton' 'string' 'New'    'callback' cb_add } ...
               { 'style' 'pushbutton' 'string' 'Rename' 'callback' cb_rename } ...
               { 'style' 'pushbutton' 'string' 'Delete' 'callback' cb_del } ...
               ...
               { 'panel' 'title' 'Edit the independent variables for this design' 'fontweight' 'bold' 'backgroundcolor' GUIBACKCOLOR } ...
               { 'style' 'text'       'string' '' 'fontweight' 'bold' } ...
               { 'style' 'pushbutton' 'string' 'New' 'callback' cb_newvar } ...
               { 'style' 'pushbutton' 'string' 'Edit'   'callback' cb_editvar } ...
               { 'style' 'pushbutton' 'string' 'Delete' 'callback' cb_delvar } ...
               { 'style' 'pushbutton' 'string' 'List factors' 'callback'  cb_list } ...
               { 'style' 'listbox'    'string' ''  'tag' 'lbfact0' 'value' 2 } ...
               { 'style' 'checkbox'   'string' ' Re-save STUDY file (if saved previously)' 'tag' 'chk_save' 'value' 1 }  };
%               { 'style' 'checkbox'   'string' 'Delete all pre-computed datafiles' 'tag' 'chk_del' 'callback' cb_lbval } };
%               { 'style' 'checkbox'  'string' 'Paired statistics' 'tag' 'lbpair0' 'callback' cb_lbval } ...
%               { 'style' 'checkbox'  'string' 'Paired statistics' 'tag' 'lbpair1' 'callback' cb_lbval } ...
%               { 'style' 'pushbutton' 'string' 'Plot'   'callback' cb_plotdmat} ...

    ht = 11.5;
    h2 = 1.1;
    h1 = 1.5;
    geometry = { {3 ht [1    1]    [2 1] } ...
                 {3 ht [2.55 1]    [1 1] } ...
                 {3 ht [3.6  1]    [0.4 1] } ...
                 ...
                 {3 ht [1    h1+1]    [1 1] } ...
                 {3 ht [1    h1+2]    [3 2.6] } ...
                 {3 ht [2.85 h1+1.15] [0.4  1] } ... % Edit
                 {3 ht [3.15 h1+1.15] [0.5  1] } ... % Raname
                 {3 ht [3.55 h1+1.15] [0.45 1] } ... % Delete
                 ...
                 {3 ht [1    h2+5.4] [3 5.1] } ...
                 {3 ht [1.07 h2+6  ] [1 1] } ...      % Indep variables 
                 {3 ht [2.20 h2+6.2] [0.4  1] } ...    % New
                 {3 ht [2.50 h2+6.2] [0.4  1] } ...    % Edit
                 {3 ht [2.80 h2+6.2] [0.45 1] } ...   % Delete
                 {3 ht [3.15 h2+6.2] [0.75 1] } ...   % List var
                 {3 ht [1.07 h2+7  ] [2.85 3] } ...   % listbox variable
                 {3 ht [1    h2+11 ] [3 1] } ...
                 };
%                 {3 ht [3.3  6+h2]   [0.5 1] } ...    % Subject text
%                 {3 ht [3.3  7+h2  ] [0.7 3] } ...    % listbox subject  
% {3 ht [1    9.8+h2] [3 1] } ...

    for i = 1:length(geometry), geometry{i}{3} = geometry{i}{3}-1; end;            
    streval = [ 'pop_studydesign(''selectdesign'', gcf);' ];    
    [tmp usrdat tmp2 result] = inputgui('uilist', uilist, 'title', 'Edit STUDY design -- pop_studydesign()', 'helpbut', 'Web help', 'helpcom',  'web(''http://sccn.ucsd.edu/wiki/Chapter_03:_Working_with_STUDY_designs'', ''-browser'')', 'geom', geometry, 'userdata', usrdat, 'eval', streval);
    if isempty(tmp), return; end
    
    % call std_makedesign
    % -------------------
    des    = usrdat.design;
    allcom = '';
    if length(des) < length(STUDY.design)
        for index = length(STUDY.design):-1:length(des)+1
            fprintf('Deleting STUDY design %d\n', index);
            com    = 'STUDY.design(index) = [];'; eval(com);
            allcom = [ allcom 10 com ];
        end
    end
    for index = 1:length(des)
        tmpdes  = des(index);
        if ~isfield(tmpdes.variable, 'vartype'), tmpdes.variable(1).vartype = []; end; 
        if index > length(STUDY.design) || ~isequal(STUDY.design(index), tmpdes)
            fprintf('Updating/creating STUDY design %d\n', index);
            
            if isfield(tmpdes, 'variable')
                for iVar = 1:length(tmpdes.variable)
                    if isempty(tmpdes.variable(iVar).vartype)
                        tmpdes.variable(iVar).vartype = 'categorical';
                    end
                end
            end
            
            [STUDY com] = std_makedesign(STUDY, ALLEEG, index, tmpdes);
            allcom = [ allcom 10 com ];
        else
            fprintf('STUDY design %d not modified\n', index);
        end
    end
    if result.listboxdesign ~= STUDY.currentdesign
        fprintf('Selecting STUDY design %d\n', result.listboxdesign);
        com = sprintf('STUDY = std_selectdesign(STUDY, ALLEEG, %d);', result.listboxdesign); eval(com);
        allcom = [ allcom 10 com ];
    end
    if result.chk_save == 1 && ~isempty(STUDY.filename)
        fprintf('Resaving STUDY\n');
        [STUDY ALLEEG com] = pop_savestudy(STUDY, ALLEEG, 'savemode', 'resave');
        allcom = [ allcom 10 com ];
    end
    if ~isempty(allcom), allcom(1) = []; end
    
elseif ischar(STUDY)
    com = STUDY;
    fig = ALLEEG;
    usrdat = get(fig, 'userdata');
    datinfo  = usrdat.datasetinfo;
    des      = usrdat.design;
    subjects = usrdat.subjects;
    filepath = usrdat.filepath;
    val = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
    
    switch com
        % summary of callbacks
        % case 'add', Add new study design        
        % case 'del', Delete study design
        % case 'rename', Rename study design
        %
        % case 'selectdesign', select a specific design
        % case 'updategui', update the study information (whenever the
        %                      user click on a button 
        %
        % case 'selectdatatrials', new GUI to select specific dataset and trials
        % case 'selectdatatrialssel', % select in the GUI above
        % case 'selectdatatrialsadd', % add new selection in the GUI above
        
        case 'selectsubj'
            strSubj = get(findobj(fig, 'tag', 'subjects'), 'string');
            if ~isempty(strSubj)
                res = str2cell(strSubj);
                selected = eeg_chaninds(struct('labels', subjects), res);
            else
                selected = 1:length(subjects);
            end
            [inds, strSubj, cellSubj] = pop_chansel(struct('labels', subjects), 'select', selected);
            if ~isempty(inds)
                if length(inds) == length(subjects)
                    strSubj = '';
                end
                set(findobj(fig, 'tag', 'subjects'), 'string', strSubj);
                disp('Warning: setting subjects for all designs');
                for iDes = 1:length(des)
                    des(iDes).cases.value = cellSubj;
                end
            end
            
        case 'setsubj'
            strSubj = get(findobj(fig, 'tag', 'subjects'), 'string');
            res = str2cell(strSubj);
            if ~isempty(setdiff(res, subjects))
                set(findobj(fig, 'tag', 'subjects'), 'string', '');
                warndlg2('Unknown subject(s) - check syntax');
                res = subjects';
            end
            for iDes = 1:length(des)
                des(iDes).cases.value = res;
            end
        
        case 'add' % Add new study design
            des(end+1) = des(val);
            des(end).variable = [];
            des(end).limo     = [];
            des(end).name = sprintf('Design %d', length(des));
            set(findobj(fig, 'tag', 'listboxdesign'), 'string', { des.name }, 'value', length(des));
            
        case 'del' % Delete study design
            val = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            if val == 1
                warndlg2('The first STUDY design cannot be removed, only modified');
                return;
            end
            des(val) = [];
            
        case 'rename' % Rename study design
            val        = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            strs       = get(findobj(fig, 'tag', 'listboxdesign'), 'string');
            result     = inputdlg2( { 'Study design name:                                                                    ' }, ...
                                      'Rename Study Design', 1,  { strs{val} }, 'pop_studydesign');
            if isempty(result), return; end
            des(val).name  = result{1};
                      
        case 'updategui' % update the study information (whenever the user click on a button)
            val = min(val, length(des));
            set(findobj(fig, 'tag', 'listboxdesign'), 'string', { des.name }, 'value', val );
            set(findobj(fig, 'tag', 'edit_storedir'), 'string', des(val).filepath);
            set(findobj(fig, 'tag', 'edit_selectdattrials'),  'string', vararg2str( des(val).include ));
            
            % update subjects
            if length(des(val).cases.value) ~= length(subjects)
                tmp = des(val).cases.value;
                tmp = cellfun(@(x) [ x ' ' ], tmp, 'uniformoutput', false);
                subjectStr = [tmp{:}];
                set(findobj(fig, 'tag', 'subjects'), 'string', subjectStr);
            end
            
            % categorical var
            curVal = {};
            for iVar = 1:length(des(val).variable)
                if ~isempty(des(val).variable(iVar).value) && strcmpi(des(val).variable(iVar).vartype, 'categorical')
                    valStr = '';
                    valStrCell = encodevals(des(val).variable(iVar).value);
                    for iVal = 1:length(valStrCell)
                        valStr = [ valStr sprintf('%s - ', valStrCell{iVal}) ];
                    end
                    valStr(end-2:end) = [];
                    strCond = sprintf('Categorical variable: %s - Values (%s)', des(val).variable(iVar).label, valStr);
                else
                    strCond = sprintf('Continuous variable: %s', des(val).variable(iVar).label);
                end
                curVal{end+1} = strCond;

            end
            curValLbfact = get(findobj(fig, 'tag', 'lbfact0'), 'value');
            valVar = min(curValLbfact, length(des(val).variable(iVar)));
            if isequal(valVar, 0), valVar = []; end
            if isempty(valVar) && ~isempty(curVal), valVar = 1; end
            set(findobj(fig, 'tag', 'lbfact0'), 'string', curVal, 'value', valVar);
            return;
            
        case 'selectdatatrials' % select specific dataset and trials
            usrdat.parent = fig;
            cb_sel = 'pop_studydesign(''selectdatatrialssel'',gcbf);';
            cb_add = 'pop_studydesign(''selectdatatrialsadd'',gcbf);';
            uilist = { { 'style' 'text'    'string' strvcat('Press ''Add'' to add data', 'selection. Multiple variables', 'are combined using AND.') } ...
                       { 'style' 'text'    'string' 'Select data based on variable', 'fontweight' 'bold' } ...
                       { 'style' 'listbox' 'string' usrdat.factors  'tag' 'lbfact2' 'callback' cb_sel 'value' 1 } ...
                       { 'style' 'text'    'string' 'Select data based on value(s)' 'fontweight' 'bold' } ...
                       { 'style' 'listbox' 'string' encodevals(usrdat.factorvals{1}) 'tag' 'lbval2' 'min' 0 'max' 2} };
            cb_renamehelp = [ 'set(findobj(gcf, ''tag'', ''help''), ''string'', ''Add'');' ...
                              'set(findobj(gcf, ''tag'', ''cancel''), ''string'', ''Erase'', ''callback'', ''set(findobj(''''tag'''', ''''edit_selectdattrials''''), ''''string'''', '''''''');'');' ...
                              'set(findobj(gcf, ''tag'', ''ok''), ''string'', ''Close'');' ];
            usrdat.fig = fig;
            inputgui('uilist', uilist, 'geometry', { [1] [1] [1] [1] [1] }, 'geomvert', [2 1 2.5 1 2.5], ...
                'helpcom', cb_add, 'userdata', usrdat, 'eval', cb_renamehelp);
            return;
            
        case 'selectdatatrialssel' % select in the GUI above
            val1  = get(findobj(fig, 'tag', 'lbfact2'), 'value');
            valfact = [1:length(usrdat.factorvals{val1})];
            tmpval = get(findobj(fig, 'tag', 'lbval2'), 'value');
            if max(tmpval) > max(valfact)
                set(findobj(fig, 'tag', 'lbval2'), 'value', valfact, 'string', encodevals(usrdat.factorvals{val1}));
            else
                set(findobj(fig, 'tag', 'lbval2'), 'string', encodevals(usrdat.factorvals{val1}), 'value', valfact);
            end
            return;

        case 'selectdatatrialsadd' % Add button in the GUI above
            val1    = get(findobj(fig, 'tag', 'lbfact2'), 'value');
            val2    = get(findobj(fig, 'tag', 'lbval2') , 'value');
            %close(fig);
            fig     = usrdat.parent;
            val     = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            des(val).include{end+1} = { usrdat.factors{val1} usrdat.factorvals{val1}(val2) };
            
        case 'updateinclude' % Add button in the GUI above
            try
                des(val).include = eval( get(findobj(fig, 'tag', 'lbfact2'), 'string') );
            catch,
                warndlg2('Syntax error');
            end
            
        case 'selectfolder'
            res = uigetdir;
            if ~isempty(findstr(filepath, res)) && findstr(filepath, res) == 1 && ~isequal(filepath, res)
                res = res(length(filepath)+2:end);
            end
            if res(1) == 0, return; end
            des(val).filepath = res;
            
        case 'list'
            val1   = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            val2   = get(findobj(fig, 'tag', 'lbfact0'), 'value');
            STD.datasetinfo = datinfo;
            STD.design = des(val1);
            STD = std_addvarlevel(STD, 1);
            pop_listfactors(STD.design);

        case 'delvar'
            val1   = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            val2   = get(findobj(fig, 'tag', 'lbfact0'), 'value');
            des(val1).variable(val2) = [];
            
        case 'newvar'
            val    = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            [tmpVar tmpVarList cat] = pop_addindepvar(usrdat);
            if ~isempty(tmpVar) && ~strcmp(tmpVar,'None')
                des(val).variable(end+1).label = tmpVar;
                des(val).variable(end  ).value = tmpVarList; % empty for cont var
                if cat, des(val).variable(end).vartype = 'categorical'; else des(val).variable(end).vartype = 'continuous'; end
            end
            % update var list
            [ usrdat.factors usrdat.factorvals usrdat.factsubj usrdat.pairing] = std_getindvar(struct('design', des, 'datasetinfo', datinfo), 'both', 1);
            usrdat.factors     = { 'None' usrdat.factors{:} };
            usrdat.factorvals  = { {}     usrdat.factorvals{:} };
            usrdat.factsubj    = { {}     usrdat.factsubj{:} };
            
        case 'editvar'
            val    = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            val2   = get(findobj(fig, 'tag', 'lbfact0'), 'value');
            tmpval = des(val).variable(val2).value;
            if strcmpi(des(val).variable(val2).vartype, 'continuous'), tmpval = []; end
            [tmpVar tmpVarList cat] = pop_addindepvar(usrdat, [], des(val).variable(val2).label, tmpval);
            if ~isempty(tmpVar) && ~strcmp(tmpVar,'None')
                des(val).variable(val2).label = tmpVar;
                des(val).variable(val2).value = tmpVarList; % empty for cont var
                if cat, des(val).variable(val2).vartype = 'categorical'; else des(val).variable(val2).vartype = 'continuous'; end
            end
            % update var list
            [ usrdat.factors usrdat.factorvals usrdat.factsubj usrdat.pairing] = std_getindvar(struct('design', des, 'datasetinfo', datinfo), 'both', 1);
            usrdat.factors     = { 'None' usrdat.factors{:} };
            usrdat.factorvals  = { {}     usrdat.factorvals{:} };
            usrdat.factsubj    = { {}     usrdat.factsubj{:} };
            
        case 'plotdmat'
            val    = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            std_plotdmat(usrdat.design(val),usrdat.datasetinfo);
            
        case 'importgvar'
            val    = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            usrdat = pop_importgroupvar(usrdat,val);
            
    end

    usrdat.design = des;
    set(fig, 'userdata', usrdat);
    pop_studydesign( 'updategui', fig);
    
    if strcmpi(com, 'newvar')
        set(findobj(fig, 'tag', 'lbfact0'), 'value', length(des(val).variable));
    end
end

function res = str2cell(strSubj)
    res = textscan(strSubj, '%s');
    res = res{1};
    for iRes = 1:length(res)
        res{iRes}(res{iRes} == '''') = [];
    end
    res = res';

function res = strmatchmult(a, b)
    if isempty(b), res = []; return; end
    res = zeros(1,length(a));
    for index = 1:length(a)
        tmpi = std_indvarmatch(a{index}, b);
        res(index) = tmpi(1); % in case there is a duplicate
    end
    %[tmp ind] = mysetdiff(b, a);
    %res = setdiff_bc([1:length(b)], ind);

function cellarray = mysort(cellarray)
    return; % was crashing for combinations of selection
            % also there is no reason the order should be different
    if ~isempty(cellarray) && ischar(cellarray{1})
        cellarray = sort(cellarray);
    end

function [cellout inds ] = mysetdiff(cell1, cell2)
    if (~isempty(cell1) && ischar(cell1{1})) || (~isempty(cell2) && ischar(cell2{1}))
         [ cellout inds ] = setdiff_bc(cell1, cell2);
    else [ cellout inds ] = setdiff_bc([ cell1{:} ], [ cell2{:} ]);
         cellout = mattocell(cellout);
    end

% encode string an numerical values for list boxes
function cellout = encodevals(cellin)
    if isempty(cellin) 
        cellout = {};
    elseif ~iscell(cellin)
        cellout = { num2str(cellin) };
    elseif ischar(cellin{1}) || iscell(cellin{1})
        for index = 1:length(cellin)
            if ischar(cellin{index})
                cellout{index} = cellin{index};
            else
                cellout{index} =  cellin{index}{1};
                for indcell = 2:length(cellin{index})
                    cellout{index} = [ cellout{index} ' & ' cellin{index}{indcell} ];
                end
            end
        end
    else
        for index = 1:length(cellin)
            if length(cellin{index}) == 1
                cellout{index} = num2str(cellin{index});
            else
                cellout{index} =  num2str(cellin{index}(1));
                for indcell = 2:length(cellin{index})
                    cellout{index} = [ cellout{index} ' & ' num2str(cellin{index}(indcell)) ];
                end
            end
        end
    end
