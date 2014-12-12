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

function [STUDY allcom] = pop_studydesign(STUDY, ALLEEG, designind, varargin);  
    
allcom = '';
if nargin < 2
    help pop_studydesign;
    return;
end;

if nargin < 3 && ~isstr(STUDY)
    
    %% create GUI
    [ usrdat.factors usrdat.factorvals usrdat.factsubj] = std_getindvar(STUDY, 'both', 1);
    
    usrdat.factors     = { 'None' usrdat.factors{:} };
    usrdat.factorvals  = { {}     usrdat.factorvals{:} };
    usrdat.factsubj    = { {}     usrdat.factsubj{:} };
    usrdat.subjects    = STUDY.subject;
    usrdat.datasetinfo = STUDY.datasetinfo;
    usrdat.design      = STUDY.design;
    usrdat.filepath    = STUDY.filepath;
    for ind = 1:length(usrdat.design)
        usrdat.design(ind).deletepreviousfiles = 0;
        for iVar = length(usrdat.design(ind).variable):-1:1
            if isempty(usrdat.design(ind).variable(iVar).label)
                usrdat.design(ind).variable(iVar) = [];
            end;
        end;
    end;
    
    % check numerical variables
    for ind1 = 1:length(usrdat.factors)
        if all(cellfun(@isnumeric, usrdat.factorvals{ind1}))
             usrdat.numerical(ind1) = 1;
        else usrdat.numerical(ind1) = 0;
        end;
    end;
        
    cb_rename = 'pop_studydesign(''rename'', gcbf);';
    cb_add    = 'pop_studydesign(''add'', gcbf);';
    cb_del    = 'pop_studydesign(''del'', gcbf);';
    cb_listboxfact1 = 'pop_studydesign(''selectfact'', gcf, 0);';
    cb_listboxfact2 = 'pop_studydesign(''selectfact'', gcf, 1);';
    cb_selectsubj   = 'pop_studydesign(''selectsubj'', gcbf);';
    cb_combinevals1 = 'pop_studydesign(''combinevals'', gcbf, 0);';
    cb_combinevals2 = 'pop_studydesign(''combinevals'', gcbf, 1);';
    cb_lbval        = 'pop_studydesign(''updategui'', gcbf);';
    cb_selectdesign = 'pop_studydesign(''selectdesign'', gcbf);';
    cb_selectdata   = 'pop_studydesign(''selectdatatrials'', gcbf);';
    cb_selectfolder = 'pop_studydesign(''selectfolder'', gcbf);';
    cb_setfolder    = 'pop_studydesign(''updategui'', gcbf);';
    cb_newvar       = 'pop_studydesign(''newvar'', gcbf);';
    cb_editvar      = 'pop_studydesign(''editvar'', gcbf);';
    cb_delvar       = 'pop_studydesign(''delvar'', gcbf);';
    cb_plotdmat     = 'pop_studydesign(''plotdmat'', gcbf);';
    
    uilist = { { 'style' 'text'       'string' 'Select STUDY design' 'fontweight' 'bold' } ...
               { 'style' 'listbox'    'string' { usrdat.design.name } 'tag' 'listboxdesign' 'callback' cb_selectdesign 'value' STUDY.currentdesign } ...
               { 'style' 'pushbutton' 'string' 'Add design'    'callback' cb_add } ...
               { 'style' 'pushbutton' 'string' 'Rename design' 'callback' cb_rename } ...
               { 'style' 'pushbutton' 'string' 'Delete design' 'callback' cb_del } ...
               { 'style' 'text'       'string' 'Subjects' 'fontweight' 'bold' } ...
               { 'style' 'text'       'string' 'Subject variables' 'fontweight' 'bold' } ...
               { 'style' 'pushbutton' 'string' 'New'    'callback' cb_newvar } ...
               { 'style' 'pushbutton' 'string' 'Edit'   'callback' cb_editvar } ...
               { 'style' 'pushbutton' 'string' 'Delete' 'callback' cb_delvar } ...
               { 'style' 'pushbutton' 'string' 'Plot'   'callback' cb_plotdmat} ...
               { 'style' 'listbox'    'string' usrdat.subjects 'tag' 'lbsubj' 'min' 0 'max' 2 'value' 1 'callback' cb_selectsubj } ...
               { 'style' 'listbox'    'string' ''  'tag' 'lbfact0' 'value' 2 } ...
               { 'style' 'listbox'    'string' ''  'tag' 'lbfact1' 'value' 1 } ... %usrdat.factors
               { 'style' 'text'       'string' 'Group variables' 'fontweight' 'bold' } ...
               { 'style' 'pushbutton' 'string' 'New' } ...
               { 'style' 'pushbutton' 'string' 'Edit' } ...
               { 'style' 'pushbutton' 'string' 'Delete' } ...
               { 'style' 'pushbutton' 'string' 'Plot'   'callback' cb_plotdmat} ...
               { 'style' 'pushbutton' 'string' 'Use only specific datasets/trials' 'callback' cb_selectdata } ...
               { 'style' 'edit'       'string' '' 'tag' 'edit_selectdattrials' 'callback' cb_lbval } ...
               { 'style' 'text'       'string' 'Store pre-computed files in folder' } ...
               { 'style' 'edit'       'string' '' 'tag' 'edit_storedir' 'callback' cb_setfolder } ...
               { 'style' 'pushbutton' 'string' '...' 'callback' cb_selectfolder } ...
               { 'style' 'checkbox'   'string' 'Delete all pre-computed datafiles associated with this specific STUDY design' 'tag' 'chk_del' 'callback' cb_lbval } ...
               { 'style' 'checkbox'   'string' 'Save the STUDY' 'tag' 'chk_save' 'value' 1 } };
%               { 'style' 'checkbox'  'string' 'Paired statistics' 'tag' 'lbpair0' 'callback' cb_lbval } ...
%               { 'style' 'checkbox'  'string' 'Paired statistics' 'tag' 'lbpair1' 'callback' cb_lbval } ...

    geometry = { {3 16 [1 1] [2 1] } ...
                 {3 16 [1 2] [2 2.6] } ...
                 {3 16 [3 2.0] [1 1] } ...
                 {3 16 [3 2.8] [1 1] } ...
                 {3 16 [3 3.6] [1 1] } ...
                 ...
                 {3 16 [1 5] [1 1] } ...
                 {3 16 [2 5] [1 1] } ...
                 {3 16 [2.80 5.2] [0.35 1] } ... % New
                 {3 16 [3.05 5.2] [0.35 1] } ... % Edit
                 {3 16 [3.30 5.2] [0.45 1] } ... % Remove
                 {3 16 [3.65 5.2] [0.35 1] } ... % Plot
                 {3 16 [1 6 ] [0.9 7.5] } ...
                 {3 16 [2 6 ] [2 3] } ...
                 {3 16 [2 10.5] [2 3] } ...
                 {3 16 [2 9.5] [1 1] } ...
                 {3 16 [2.80 9.7] [0.35 1] } ... % New
                 {3 16 [3.05 9.7] [0.35 1] } ... % Edit
                 {3 16 [3.30 9.7] [0.45 1] } ... % Remove
                 {3 16 [3.65 9.7] [0.35 1] } ... % Plot
                 {3 16 [1    14] [1.3 1] } ...
                 {3 16 [2.25 14] [1.7 1] } ...
                 {3 16 [1    15] [1.3 1] } ...
                 {3 16 [2.25 15] [1.2 1] } ...
                 {3 16 [3.45 15] [0.55 1] } ...
                 {3 16 [1 16] [3 1] } ...
                 {3 16 [1 16.8] [3 1] } ...
                 };

    for i = 1:length(geometry), geometry{i}{3} = geometry{i}{3}-1; end;            
    streval = [ 'pop_studydesign(''selectdesign'', gcf);' ];    
    [tmp usrdat tmp2 result] = inputgui('uilist', uilist, 'title', 'Edit STUDY design -- pop_studydesign()', 'helpbut', 'Web help', 'helpcom',  'web(''http://sccn.ucsd.edu/wiki/Chapter_03:_Working_with_STUDY_designs'', ''-browser'')', 'geom', geometry, 'userdata', usrdat, 'eval', streval);
    if isempty(tmp), return; end;
    
    % call std_makedesign
    % -------------------
    des    = usrdat.design;
    allcom = '';
    if length(des) < length(STUDY.design)
        for index = length(des)+1:length(STUDY.design)
            fprintf('Deleting STUDY design %d\n', index);
            com    = 'STUDY.design(index).name = '';'; eval(com);
            allcom = [ allcom 10 com ];
        end;
    end;
    for index = 1:length(des)
        tmpdes  = rmfield(des(index), 'deletepreviousfiles');
        rmfiles = fastif(des(index).deletepreviousfiles, 'limited', 'off');
        if index > length(STUDY.design) || ~isequal(STUDY.design(index), tmpdes) || strcmpi(rmfiles, 'on')
            fprintf('Updating/creating STUDY design %d\n', index);
            
            % test if file exist and issue warning
            if length(STUDY.design) >= index && isfield('cell', STUDY.design) && ~isempty(STUDY.design(index).cell) && ...
                ~isempty(dir([ STUDY.design(index).cell(1).filebase '.*' ])) &&  strcmpi(rmfiles, 'off')
                if ~isequal(tmpdes.variable(1).label, STUDY.design(index).variable(1).label) || ...
                     ~isequal(tmpdes.variable(2).label, STUDY.design(index).variable(2).label) || ...
                      ~isequal(tmpdes.include, STUDY.design(index).include) || ...
                       ~isequal(tmpdes.variable(1).value, STUDY.design(index).variable(1).value) || ...
                        ~isequal(tmpdes.cases.value, STUDY.design(index).cases.value) || ...
                         ~isequal(tmpdes.variable(2).value, STUDY.design(index).variable(2).value)
                    res = questdlg2(strvcat([ 'Precomputed data files exist for design ' int2str(index) '.' ], ' ', ...
                                       'Modifying this design without deleting the associated files', ...
                                       'might mean that they will stay on disk and will be unusable'), ...
                                       'STUDY design warning', 'Abort', 'Continue', 'Continue');
                    if strcmpi(res, 'Abort'), return; end;
                end;
            end;
            [STUDY com] = std_makedesign(STUDY, ALLEEG, index, tmpdes, 'delfiles', rmfiles);
            allcom = [ allcom 10 com ];
        else
            fprintf('STUDY design %d not modified\n', index);
        end;
    end;
    if result.listboxdesign ~= STUDY.currentdesign
        fprintf('Selecting STUDY design %d\n', result.listboxdesign);
        com = sprintf('STUDY = std_selectdesign(STUDY, ALLEEG, %d);', result.listboxdesign); eval(com);
        allcom = [ allcom 10 com ];
    end;
    if result.chk_save == 1
        fprintf('Resaving STUDY\n');
        [STUDY ALLEEG com] = pop_savestudy(STUDY, ALLEEG, 'savemode', 'resave');
        allcom = [ allcom 10 com ];
    end;
    if ~isempty(allcom), allcom(1) = []; end;
    
elseif isstr(STUDY)
    com = STUDY;
    fig = ALLEEG;
    usrdat = get(fig, 'userdata');
    datinfo  = usrdat.datasetinfo;
    des      = usrdat.design;
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
        % case 'selectfact', select a specific ind. var. (update value listboxes)
        % case 'combinevals', combine values in value listboxes
        % case 'selectsubj', select specific subjects
        %
        % case 'selectdatatrials', new GUI to select specific dataset and trials
        % case 'selectdatatrialssel', % select in the GUI above
        % case 'selectdatatrialsadd', % add new selection in the GUI above
        
        case 'add', % Add new study design
            des(end+1) = des(val);
            des(end).name = sprintf('Design %d', length(des));
            set(findobj(fig, 'tag', 'listboxdesign'), 'string', { des.name }, 'value', length(des));
            
        case 'del', % Delete study design
            val = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            if val == 1
                warndlg2('The first STUDY design cannot be removed, only modified');
                return;
            end;
            des(val) = [];
            
        case 'rename', % Rename study design
            val        = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            strs       = get(findobj(fig, 'tag', 'listboxdesign'), 'string');
            result     = inputdlg2( { 'Study design name:                                                                    ' }, ...
                                      'Rename Study Design', 1,  { strs{val} }, 'pop_studydesign');
            if isempty(result), return; end;
            des(val).name  = result{1};
                      
        case 'updategui', % update the study information (whenever the user click on a button)
            val = min(val, length(des));
            set(findobj(fig, 'tag', 'listboxdesign'), 'string', { des.name }, 'value', val );
            set(findobj(fig, 'tag', 'chk_del'), 'value', des(val).deletepreviousfiles);
            set(findobj(fig, 'tag', 'edit_storedir'), 'string', des(val).filepath);
            set(findobj(fig, 'tag', 'edit_selectdattrials'),  'string', vararg2str( des(val).include ));
            
            % update subjects
            cellSubj = get(findobj(fig, 'tag', 'lbsubj'), 'string');
            [tmp valSubj] = intersect(cellSubj, des(val).cases.value);
            set(findobj(fig, 'tag', 'lbsubj'), 'value', valSubj);
            
            % categorical var
            curVal = {};
            for iVar = 1:length(des(val).variable)
                if ~isempty(des(val).variable(iVar).value)
                    valStr = '';
                    valStrCell = encodevals(des(val).variable(iVar).value);
                    for iVal = 1:length(valStrCell)
                        valStr = [ valStr sprintf('%s - ', valStrCell{iVal}) ];
                    end;
                    valStr(end-2:end) = [];
                    strCond = sprintf('Categorical variable: %s - Values (%s)', des(val).variable(iVar).label, valStr);
                else
                    strCond = sprintf('Continuous variable: %s', des(val).variable(iVar).label);
                end;
                curVal{end+1} = strCond;

            end;
            curValLbfact = get(findobj(fig, 'tag', 'lbfact0'), 'value');
            valVar = min(curValLbfact, length(des(val).variable(iVar)));
            if isequal(valVar, 0), valVar = []; end;
            if isempty(valVar) && ~isempty(curVal), valVar = 1; end;
            set(findobj(fig, 'tag', 'lbfact0'), 'string', curVal, 'value', valVar);
            return;
            
        case 'selectdatatrials', % select specific dataset and trials
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
            
        case 'selectdatatrialssel', % select in the GUI above
            val1  = get(findobj(fig, 'tag', 'lbfact2'), 'value');
            valfact = [1:length(usrdat.factorvals{val1})];
            tmpval = get(findobj(fig, 'tag', 'lbval2'), 'value');
            if max(tmpval) > max(valfact)
                set(findobj(fig, 'tag', 'lbval2'), 'value', valfact, 'string', encodevals(usrdat.factorvals{val1}));
            else
                set(findobj(fig, 'tag', 'lbval2'), 'string', encodevals(usrdat.factorvals{val1}), 'value', valfact);
            end;
            return;

        case 'selectdatatrialsadd', % Add button in the GUI above
            val1    = get(findobj(fig, 'tag', 'lbfact2'), 'value');
            val2    = get(findobj(fig, 'tag', 'lbval2') , 'value');
            %close(fig);
            fig     = usrdat.parent;
            val     = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            des(val).include{end+1} = { usrdat.factors{val1} usrdat.factorvals{val1}(val2) };
                        
        case 'selectsubj', % select specific subjects
            valSubj  = get(findobj(fig, 'tag', 'lbsubj'), 'value');
            cellSubj = get(findobj(fig, 'tag', 'lbsubj'), 'string');
            des(val).cases.value = cellSubj(valSubj);
            
        case 'updateinclude', % Add button in the GUI above
            try
                des(val).include = eval( get(findobj(fig, 'tag', 'lbfact2'), 'string') );
            catch,
                warndlg2('Syntax error');
            end;
            
        case 'selectfolder',
            res = uigetdir;
            if ~isempty(findstr(filepath, res)) && findstr(filepath, res) == 1 && ~isequal(filepath, res)
                res = res(length(filepath)+2:end);
            end;
            if res(1) == 0, return; end;
            des(val).filepath = res;
            
        case 'delvar'
            val1   = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            val2   = get(findobj(fig, 'tag', 'lbfact0'), 'value');
            des(val1).variable(val2) = [];
            
        case 'newvar'
            val    = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            [tmpVar tmpVarList] = pop_addindepvar(usrdat);
            if ~isempty(tmpVar)
                des(val).variable(end+1).label = tmpVar;
                des(val).variable(end  ).value = tmpVarList; % empty for cont var
            end;
            
        case 'editvar'
            val    = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            val2   = get(findobj(fig, 'tag', 'lbfact0'), 'value');
            [tmpVar tmpVarList] = pop_addindepvar(usrdat, [], des(val).variable(val2).label, des(val).variable(val2).value);
            if ~isempty(tmpVar)
                des(val).variable(val2).label = tmpVar;
                des(val).variable(val2).value = tmpVarList; % empty for cont var
            end;
        case 'plotdmat'
            val    = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            std_plotdmat(usrdat,val);
    end;

    usrdat.design = des;
    set(fig, 'userdata', usrdat);
    pop_studydesign( 'updategui', fig);
    
    if strcmpi(com, 'newvar')
        set(findobj(fig, 'tag', 'lbfact0'), 'value', length(des(val).variable));
    end;
end;

function res = strmatchmult(a, b);
    if isempty(b), res = []; return; end;
    res = zeros(1,length(a));
    for index = 1:length(a)
        tmpi = std_indvarmatch(a{index}, b);
        res(index) = tmpi(1); % in case there is a duplicate
    end;
    %[tmp ind] = mysetdiff(b, a);
    %res = setdiff_bc([1:length(b)], ind);

function cellarray = mysort(cellarray)
    return; % was crashing for combinations of selection
            % also there is no reason the order should be different
    if ~isempty(cellarray) && isstr(cellarray{1})
        cellarray = sort(cellarray);
    end;

function [cellout inds ] = mysetdiff(cell1, cell2);
    if (~isempty(cell1) && isstr(cell1{1})) || (~isempty(cell2) && isstr(cell2{1}))
         [ cellout inds ] = setdiff_bc(cell1, cell2);
    else [ cellout inds ] = setdiff_bc([ cell1{:} ], [ cell2{:} ]);
         cellout = mattocell(cellout);
    end;

% encode string an numerical values for list boxes
function cellout = encodevals(cellin)
    if isempty(cellin) 
        cellout = {};
    elseif ~iscell(cellin)
        cellout = { num2str(cellin) };
    elseif ischar(cellin{1}) || iscell(cellin{1})
        for index = 1:length(cellin)
            if isstr(cellin{index})
                cellout{index} = cellin{index};
            else
                cellout{index} =  cellin{index}{1};
                for indcell = 2:length(cellin{index})
                    cellout{index} = [ cellout{index} ' & ' cellin{index}{indcell} ];
                end;
            end;
        end;
    else
        for index = 1:length(cellin)
            if length(cellin{index}) == 1
                cellout{index} = num2str(cellin{index});
            else
                cellout{index} =  num2str(cellin{index}(1));
                for indcell = 2:length(cellin{index})
                    cellout{index} = [ cellout{index} ' & ' num2str(cellin{index}(indcell)) ];
                end;
            end;
        end;
    end;
