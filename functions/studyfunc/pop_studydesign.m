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
    end;
    
    % build menu
    popupselectsubj = { 'Select all subjects' };
    for ind1 = 1:length(usrdat.factors)
        if ~isempty(usrdat.factsubj{ind1})
            if any(cellfun(@length, usrdat.factsubj{ind1}) ~= length(usrdat.subjects))
                for ind2 = 1:length(usrdat.factorvals{ind1})
                    if ~iscell(usrdat.factorvals{ind1}{ind2}) % not a combined value
                        tmpval = encodevals(usrdat.factorvals{ind1}(ind2));
                        popupselectsubj{end+1} = [ num2str(usrdat.factors{ind1}) ' - ' tmpval{1} ];
                    end;
                end;
            end;
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
    cb_lbval        = 'pop_studydesign(''updatedesign'', gcbf);';
    cb_selectdesign = 'pop_studydesign(''selectdesign'', gcbf);';
    cb_selectdata   = 'pop_studydesign(''selectdatatrials'', gcbf);';
    cb_selectfolder = 'pop_studydesign(''selectfolder'', gcbf);';
    cb_setfolder    = 'pop_studydesign(''updatedesign'', gcbf);';
    uilist = { { 'style' 'text'       'string' 'Select STUDY design' 'fontweight' 'bold' } ...
               { 'style' 'listbox'    'string' { usrdat.design.name } 'tag' 'listboxdesign' 'callback' cb_selectdesign 'value' STUDY.currentdesign } ...
               { 'style' 'pushbutton' 'string' 'Add design'    'callback' cb_add } ...
               { 'style' 'pushbutton' 'string' 'Rename design' 'callback' cb_rename } ...
               { 'style' 'pushbutton' 'string' 'Delete design' 'callback' cb_del } ...
               { 'style' 'text'       'string' 'Subjects' 'fontweight' 'bold' } ...
               { 'style' 'text'       'string' 'Independent variable 1   ' 'fontweight' 'bold' } ...
               { 'style' 'text'       'string' 'Independent variable 2   ' 'fontweight' 'bold' } ...
               { 'style' 'listbox'    'string' usrdat.subjects 'tag' 'lbsubj' 'min' 0 'max' 2 'value' 1 'callback' cb_lbval } ...
               { 'style' 'listbox'    'string' usrdat.factors  'tag' 'lbfact0' 'callback' cb_listboxfact1 'value' 2 } ...
               { 'style' 'listbox'    'string' usrdat.factors  'tag' 'lbfact1' 'callback' cb_listboxfact2 'value' 1 } ...
               { 'style' 'text'       'string' 'Ind. var. 1 values ' } ...
               { 'style' 'text'       'string' 'Ind. var. 2 values' } ...
               { 'style' 'listbox'    'string' '' 'tag' 'lbval0' 'min' 0 'max' 2 'callback' cb_lbval } ...
               { 'style' 'listbox'    'string' '' 'tag' 'lbval1' 'min' 0 'max' 2 'callback' cb_lbval } ...
               { 'style' 'popupmenu'  'string' popupselectsubj 'tag' 'popupselect' 'callback' cb_selectsubj } ...
               { 'style' 'pushbutton' 'string' 'Combine selected values' 'tag' 'combine1' 'callback' cb_combinevals1 } ...
               { 'style' 'pushbutton' 'string' 'Combine selected values' 'tag' 'combine2' 'callback' cb_combinevals2 } ...
               { 'style' 'popupmenu'  'string' 'Paired statistics|Unpaired statistics' 'tag' 'lbpair0' 'callback' cb_lbval } ...
               { 'style' 'popupmenu'  'string' 'Paired statistics|Unpaired statistics' 'tag' 'lbpair1' 'callback' cb_lbval } ...
               { 'style' 'pushbutton' 'string' 'Use only specific datasets/trials' 'callback' cb_selectdata } ...
               { 'style' 'edit'       'string' '' 'tag' 'edit_selectdattrials' 'callback' cb_lbval } ...
               { 'style' 'text'       'string' 'Store pre-computed files in folder' } ...
               { 'style' 'edit'       'string' '' 'tag' 'edit_storedir' 'callback' cb_setfolder } ...
               { 'style' 'pushbutton' 'string' '...' 'callback' cb_selectfolder } ...
               { 'style' 'checkbox'   'string' 'Delete all pre-computed datafiles associated with this specific STUDY design' 'tag' 'chk_del' 'callback' cb_lbval } ...
               { 'style' 'checkbox'   'string' 'Save the STUDY' 'tag' 'chk_save' 'value' 1 } };
%               { 'style' 'checkbox'  'string' 'Paired statistics' 'tag' 'lbpair0' 'callback' cb_lbval } ...
%               { 'style' 'checkbox'  'string' 'Paired statistics' 'tag' 'lbpair1' 'callback' cb_lbval } ...

    geometry = { {3 18 [1 1] [2 1] } ...
                 {3 18 [1 2] [2 3] } ...
                 {3 18 [3 2] [1 1] } ...
                 {3 18 [3 3] [1 1] } ...
                 {3 18 [3 4] [1 1] } ...
                 {3 18 [1 5] [1 1] } ...
                 {3 18 [2 5] [1 1] } ...
                 {3 18 [3 5] [1 1] } ...
                 {3 18 [1 6] [1 8] } ...
                 {3 18 [2 6] [0.98 3] } ...
                 {3 18 [3 6] [0.98 3] } ...
                 {3 18 [2 9] [1 1] } ...
                 {3 18 [3 9] [1 1] } ...
                 {3 18 [2 10] [0.98 3] } ...
                 {3 18 [3 10] [0.98 3] } ...
                 {3 18 [1 14] [1 1] } ...
                 {3 18 [2 13] [1 1] } ...
                 {3 18 [3 13] [1 1] } ...
                 {3 18 [2 14] [1 1] } ...
                 {3 18 [3 14] [1 1] } ...
                 {3 18 [1    15.5] [1.3 1] } ...
                 {3 18 [2.25 15.5] [1.7 1] } ...
                 {3 18 [1    16.5] [1.3 1] } ...
                 {3 18 [2.25 16.5] [1.2 1] } ...
                 {3 18 [3.45 16.5] [0.55 1] } ...
                 {3 18 [1 17.5] [3 1] } ...
                 {3 18 [1 19] [3 1] } ...
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
    
    switch com
        % summary of callbacks
        % case 'add', Add new study design        
        % case 'del', Delete study design
        % case 'rename', Rename study design
        %
        % case 'selectdesign', select a specific design
        % case 'updatedesign', update the study information (whenever the
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
            inde  = find( cellfun(@isempty,{ des.name }));
            if isempty(inde), inde = length(des)+1; end;
            des(inde(1)) = des(1);
            des(inde(1)).name = sprintf('Design %d', inde(1));
            set(findobj(fig, 'tag', 'listboxdesign'), 'string', { des.name }, 'value', inde(1));
            usrdat.design = des;
        	set(fig, 'userdata', usrdat);
            pop_studydesign( 'selectstudy', fig);
            return;
            
        case 'del', % Delete study design
            val = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            if val == 1
                warndlg2('The first STUDY design cannot be removed, only modified');
                return;
            end;
            des(val).name = '';
            set(findobj(fig, 'tag', 'listboxdesign'), 'value', 1, 'string', { des.name } );
            usrdat.design = des;
        	set(fig, 'userdata', usrdat);
            pop_studydesign( 'selectstudy', fig);
            return;
            
        case 'rename', % Rename study design
            val        = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            strs       = get(findobj(fig, 'tag', 'listboxdesign'), 'string');
            result     = inputdlg2( { 'Study design name:                                                                    ' }, ...
                                      'Rename Study Design', 1,  { strs{val} }, 'pop_studydesign');
            if isempty(result), return; end;
            des(val).name  = result{1};
            set(findobj(fig, 'tag', 'listboxdesign'), 'string', { des.name } );
            
        case 'selectdesign', % select a specific design
            val  = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            if isempty(des(val).name)
                set(findobj(fig, 'tag', 'lbfact0'), 'string', '', 'value', 1);
                set(findobj(fig, 'tag', 'lbfact1'), 'string', '', 'value', 1);
                set(findobj(fig, 'tag', 'lbval0') , 'string', '', 'value', 1);
                set(findobj(fig, 'tag', 'lbval1') , 'string', '', 'value', 1);
                set(findobj(fig, 'tag', 'lbpair0'), 'value', 1);
                set(findobj(fig, 'tag', 'lbpair1'), 'value', 1);
                set(findobj(fig, 'tag', 'lbsubj') , 'value' , 1, 'string', '');
                set(findobj(fig, 'tag', 'chk_del'), 'value', des(val).deletepreviousfiles );
                set(findobj(fig, 'tag', 'edit_selectdattrials'), 'string', '' );
                % do not change file path
                return;
            end;
            val1 = strmatch(des(val).variable(1).label, usrdat.factors, 'exact'); if isempty(val1), val1 = 1; end;
            val2 = strmatch(des(val).variable(2).label, usrdat.factors, 'exact'); if isempty(val2), val2 = 1; end;
            set(findobj(fig, 'tag', 'lbfact0'), 'string', usrdat.factors, 'value', val1);
            set(findobj(fig, 'tag', 'lbfact1'), 'string', usrdat.factors, 'value', val2);
            valfact1 = strmatchmult(des(val).variable(1).value, usrdat.factorvals{val1});
            valfact2 = strmatchmult(des(val).variable(2).value, usrdat.factorvals{val2});
            if isempty(valfact1), listboxtop1 = 1; else listboxtop1 = valfact1(1); end;
            if isempty(valfact2), listboxtop2 = 1; else listboxtop2 = valfact2(1); end;
            set(findobj(fig, 'tag', 'lbval0'), 'string', encodevals(usrdat.factorvals{val1}), 'value', valfact1, 'listboxtop', listboxtop1);
            set(findobj(fig, 'tag', 'lbval1'), 'string', encodevals(usrdat.factorvals{val2}), 'value', valfact2, 'listboxtop', listboxtop2);
            valsubj = strmatchmult(des(val).cases.value, usrdat.subjects);
            set(findobj(fig, 'tag', 'lbsubj'), 'string', usrdat.subjects, 'value', valsubj);
            if isempty(des(val).include), str = ''; else str = vararg2str(des(val).include); end;
            set(findobj(fig, 'tag', 'chk_del'), 'value', des(val).deletepreviousfiles );
            set(findobj(fig, 'tag', 'edit_selectdattrials'), 'string', str );
            set(findobj(fig, 'tag', 'popupselect'), 'value', 1 );
            set(findobj(fig, 'tag', 'lbpair0'), 'value', fastif(isequal(des(val).variable(1).pairing,'on'),1,2));
            set(findobj(fig, 'tag', 'lbpair1'), 'value', fastif(isequal(des(val).variable(2).pairing,'on'),1,2));
            if ~isfield(des, 'filepath') || isempty(des(val).filepath), des(val).filepath = ''; end;
            set(findobj(fig, 'tag', 'edit_storedir'), 'string', des(val).filepath);
            
        case 'updatedesign', % update the study information (whenever the user click on a button)
            val    = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            val1   = get(findobj(fig, 'tag', 'lbfact0'), 'value');
            val2   = get(findobj(fig, 'tag', 'lbfact1'), 'value');
            valf1  = get(findobj(fig, 'tag', 'lbval0'),  'value');
            valf2  = get(findobj(fig, 'tag', 'lbval1'),  'value');
            valp1  = get(findobj(fig, 'tag', 'lbpair0'), 'value');
            valp2  = get(findobj(fig, 'tag', 'lbpair1'), 'value');
            vals   = get(findobj(fig, 'tag', 'lbsubj'),  'value');
            valchk = get(findobj(fig, 'tag', 'chk_del'), 'value');
            filep  = get(findobj(fig, 'tag', 'edit_storedir'), 'string');
            strs   = get(findobj(fig, 'tag', 'edit_selectdattrials'),  'string');
            valpaired = { 'on' 'off' };
            
            if ~strcmpi(des(val).variable(1).label, usrdat.factors{val1})
                des(val).variable(1).label = usrdat.factors{val1};
                des(val).variable(1).value = usrdat.factorvals{val1}(valf1);
            end;
            if ~strcmpi(des(val).variable(2).label, usrdat.factors{val2})
                des(val).variable(2).label = usrdat.factors{val2};
                des(val).variable(2).value = usrdat.factorvals{val2}(valf2);
            end;
            if ~isequal(mysort(des(val).variable(1).value), mysort(usrdat.factorvals{val1}(valf1)))
                des(val).variable(1).value = usrdat.factorvals{val1}(valf1);
            end;
            if ~isequal(mysort(des(val).variable(2).value), mysort(usrdat.factorvals{val2}(valf2)))
                des(val).variable(2).value = usrdat.factorvals{val2}(valf2);
            end;
            if ~isequal(mysort(des(val).cases.value), mysort(usrdat.subjects(vals)))
                des(val).cases.value = usrdat.subjects(vals);
            end;
            if ~isequal(des(val).variable(1).pairing, valpaired{valp1})
                des(val).variable(1).pairing = valpaired{valp1};
            end;
            if ~isequal(des(val).variable(2).pairing, valpaired{valp2})
                des(val).variable(2).pairing = valpaired{valp2};
            end;
            if ~isequal(des(val).deletepreviousfiles, valchk)
                des(val).deletepreviousfiles = 1;
            end;
            if ~isfield(des, 'filepath') || ~isequal(des(val).filepath, filep)
                des(val).filepath = filep;
            end;
            if ~isequal(des(val).include, strs)
                try,
                    des(val).include = eval( [ '{' strs '}' ]);
                catch,
                    disp('Error while decoding list of parameters');
                    des(val).include = {};
                    set(findobj(fig, 'tag', 'edit_selectdattrials'),  'string', '');
                end;
            end;
            
        case 'selectfact', % select a specific ind. var. (update value listboxes)
            factval = designind;
            val   = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            val1  = get(findobj(fig, 'tag', [ 'lbfact' num2str(factval) ]), 'value');
            val2  = get(findobj(fig, 'tag', [ 'lbfact' num2str(~factval) ]), 'value');
%             if val1 == val2 && val1 ~= 1
%                 warndlg2('Cannot select twice the same independent variable');
%                 val1 = 1;
%                 set(findobj(fig, 'tag', [ 'lbfact' num2str(factval) ]), 'value', val1);
%             end;
            valfact = [1:length(usrdat.factorvals{val1})];
            set(findobj(fig, 'tag', ['lbval' num2str(factval) ]), 'string', encodevals(usrdat.factorvals{val1}), 'value', valfact, 'listboxtop', 1);
            pop_studydesign('updatedesign', fig);
            return;
            
        case 'combinevals', %  combine values in value listboxes
            factval = designind;
            val1    = get(findobj(fig, 'tag', [ 'lbfact' num2str(factval) ]), 'value');
            vals    = get(findobj(fig, 'tag', [ 'lbval'  num2str(factval) ]), 'value');
            strs    = get(findobj(fig, 'tag', [ 'lbval'  num2str(factval) ]), 'string');
            if length(vals) == 1
                warndlg2('You need to select several values to combine them');
                return;
            end;
            if ~iscell(usrdat.factorvals{val1})
                warndlg2('Cannot combine values from numerical variables');
                return;
            end;
            % combine values for string and integers
            if isstr(usrdat.factorvals{val1}{1}) || iscell(usrdat.factorvals{val1}{1})
                tmpcell = {};
                for indCell = vals(:)'
                    if iscell(usrdat.factorvals{val1}{indCell})
                        tmpcell = { tmpcell{:} usrdat.factorvals{val1}{indCell}{:} };
                    else
                        tmpcell = { tmpcell{:} usrdat.factorvals{val1}{indCell} };
                    end;
                end;
                usrdat.factorvals{val1}{end+1} = unique_bc(tmpcell);
            else
                usrdat.factorvals{val1}{end+1} = unique_bc([ usrdat.factorvals{val1}{vals} ]);
            end;
            set(findobj(fig, 'tag', ['lbval' num2str(factval) ]), 'string', encodevals(usrdat.factorvals{val1}));
            
        case 'selectsubj', % select specific subjects
            val = get(findobj(fig, 'tag', 'popupselect'), 'value');
            str = get(findobj(fig, 'tag', 'popupselect'), 'string');
            str = str{val};
            if val == 1, 
                subjbox = get(findobj(fig, 'tag', 'lbsubj'), 'string');
                set(findobj(fig, 'tag', 'lbsubj'), 'value', [1:length(subjbox)]);
                return;
            end;
            indunders = findstr( ' - ', str);
            factor    = str(1:indunders-1);
            factorval = str(indunders+3:end);
            
            % select subjects
            eval( [ 'allsetvals = { datinfo.' factor '};' ]);
            indset = strmatch(factorval, allsetvals, 'exact');
            subjects = unique_bc( { datinfo(indset).subject } );
            
            % change the subject listbox
            val     = get(findobj(fig, 'tag', 'popupselect'), 'value');
            subjbox = get(findobj(fig, 'tag', 'lbsubj'), 'string');
            indsubj = [];
            for ind = 1:length(subjects);
                indsubj(ind) = strmatch(subjects{ind}, subjbox, 'exact');
            end;
            set(findobj(fig, 'tag', 'lbsubj'), 'value', indsubj);
            set(findobj(fig, 'tag', 'popupselect'), 'value', 1);
            pop_studydesign('updatedesign', fig);
            return;
            %set(findobj(get(gcbf, ''userdata''), ''tag'', ''edit_selectdattrials''
            
        case 'selectdatatrials', % select specific dataset and trials
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
            pop_studydesign('updatedesign', fig);
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

        case 'selectfolder',
            res = uigetdir;
            if ~isempty(findstr(filepath, res)) && findstr(filepath, res) == 1 && ~isequal(filepath, res)
                res = res(length(filepath)+2:end);
            end;
            if res(1) == 0, return; end;
            set(findobj(fig, 'tag', 'edit_storedir'), 'string', res);
            pop_studydesign('updatedesign', fig);
            return;
            
        case 'selectdatatrialsadd', % Add button in the GUI above
            val1    = get(findobj(fig, 'tag', 'lbfact2'), 'value');
            val2    = get(findobj(fig, 'tag', 'lbval2') , 'value');
            objedit = findobj(usrdat.fig, 'tag', 'edit_selectdattrials');
            str     = get(objedit, 'string');
            if ~isempty(str), str = [ str ',' ]; end;
            set(objedit, 'string', [ str vararg2str( { usrdat.factors{val1} usrdat.factorvals{val1}(val2) }) ]);
            return;
            
    end;
    usrdat.design = des;
	set(fig, 'userdata', usrdat);
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
