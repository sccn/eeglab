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
    for ind = 1:length(usrdat.design)
        usrdat.design(ind).deletepreviousfiles = 0;
    end;
    
    % build menu
    popupselectsubj = { 'Select all subjects' };
    for ind1 = 1:length(usrdat.factors)
        if ~isempty(usrdat.factsubj{ind1})
            if any(cellfun(@length, usrdat.factsubj{ind1}) ~= length(usrdat.subjects))
                for ind2 = 1:length(usrdat.factorvals{ind1})
                    popupselectsubj{end+1} = [ usrdat.factors{ind1} ' - ' usrdat.factorvals{ind1}{ind2} ];
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
               { 'style' 'checkbox'   'string' 'Delete all datafiles associated with this STUDY design' 'tag' 'chk_del' 'callback' cb_lbval } ...
               { 'style' 'checkbox'   'string' 'Save the STUDY' 'tag' 'chk_save' 'value' 1 } };
%               { 'style' 'checkbox'  'string' 'Paired statistics' 'tag' 'lbpair0' 'callback' cb_lbval } ...
%               { 'style' 'checkbox'  'string' 'Paired statistics' 'tag' 'lbpair1' 'callback' cb_lbval } ...

    geometry = { {3 17 [1 1] [2 1] } ...
                 {3 17 [1 2] [2 3] } ...
                 {3 17 [3 2] [1 1] } ...
                 {3 17 [3 3] [1 1] } ...
                 {3 17 [3 4] [1 1] } ...
                 {3 17 [1 5] [1 1] } ...
                 {3 17 [2 5] [1 1] } ...
                 {3 17 [3 5] [1 1] } ...
                 {3 17 [1 6] [1 8] } ...
                 {3 17 [2 6] [0.98 3] } ...
                 {3 17 [3 6] [0.98 3] } ...
                 {3 17 [2 9] [1 1] } ...
                 {3 17 [3 9] [1 1] } ...
                 {3 17 [2 10] [0.98 3] } ...
                 {3 17 [3 10] [0.98 3] } ...
                 {3 17 [1 14] [1 1] } ...
                 {3 17 [2 13] [1 1] } ...
                 {3 17 [3 13] [1 1] } ...
                 {3 17 [2 14] [1 1] } ...
                 {3 17 [3 14] [1 1] } ...
                 {3 17 [1 15.5] [1.3 1] } ...
                 {3 17 [2.25 15.5] [1.7 1] } ...
                 {3 17 [1 16.5] [3 1] } ...
                 {3 17 [1 18] [3 1] } ...
                 };

    for i = 1:length(geometry), geometry{i}{3} = geometry{i}{3}-1; end;            
    streval = [ 'pop_studydesign(''selectdesign'', gcf);' ];    
    [tmp usrdat tmp2 result] = inputgui('uilist', uilist, 'geom', geometry, 'userdata', usrdat, 'eval', streval);
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
        rmfiles = fastif(des(index).deletepreviousfiles, 'on', 'off');
        if index > length(STUDY.design) || ~isequal(STUDY.design(index), tmpdes) || strcmpi(rmfiles, 'on')
            fprintf('Updating/creating STUDY design %d\n', index);
            
            % test if file exist and issue warning
            if ~isempty(dir([ STUDY.design(index).setinfo(1).filebase '.*' ])) &&  strcmpi(rmfiles, 'off')
                res = questdlg2(strvcat([ 'Precomputed data files exist for design ' int2str(index) '.' ], ' ', ...
                                       'Modifying this design without deleting the associated files', ...
                                       'might mean that they will stay on disk and will be unusable'), ...
                                       'STUDY design warning', 'Abord', 'Continue', 'Continue');
                if strcmpi(res, 'Abord'), return; end;
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
    datinfo = usrdat.datasetinfo;
    des = usrdat.design;
    
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
                return;
            end;
            val1 = strmatch(des(val).indvar1, usrdat.factors, 'exact'); if isempty(val1), val1 = 1; end;
            val2 = strmatch(des(val).indvar2, usrdat.factors, 'exact'); if isempty(val2), val2 = 1; end;
            set(findobj(fig, 'tag', 'lbfact0'), 'string', usrdat.factors, 'value', val1);
            set(findobj(fig, 'tag', 'lbfact1'), 'string', usrdat.factors, 'value', val2);
            valfact1 = strmatchmult(des(val).condition, usrdat.factorvals{val1});
            valfact2 = strmatchmult(des(val).group    , usrdat.factorvals{val2});
            set(findobj(fig, 'tag', 'lbval0'), 'string', usrdat.factorvals{val1}, 'value', valfact1);
            set(findobj(fig, 'tag', 'lbval1'), 'string', usrdat.factorvals{val2}, 'value', valfact2);
            valsubj = strmatchmult(des(val).subject    , usrdat.subjects);
            set(findobj(fig, 'tag', 'lbsubj'), 'string', usrdat.subjects, 'value', valsubj);
            if isempty(des(val).includevarlist), str = ''; else str = vararg2str(des(val).includevarlist); end;
            set(findobj(fig, 'tag', 'chk_del'), 'value', des(val).deletepreviousfiles );
            set(findobj(fig, 'tag', 'edit_selectdattrials'), 'string', str );
            set(findobj(fig, 'tag', 'popupselect'), 'value', 1 );
            set(findobj(fig, 'tag', 'lbpair0'), 'value', fastif(isequal(des(val).statvar1,'paired'),1,2));
            set(findobj(fig, 'tag', 'lbpair1'), 'value', fastif(isequal(des(val).statvar2,'paired'),1,2));
            
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
            strs   = get(findobj(fig, 'tag', 'edit_selectdattrials'),  'string');
            valpaired = { 'paired' 'unpaired' };
            
            if ~strcmpi(des(val).indvar1, usrdat.factors{val1})
                des(val).indvar1   = usrdat.factors{val1};
                des(val).condition = usrdat.factorvals{val1}(valf1);
            end;
            if ~strcmpi(des(val).indvar2, usrdat.factors{val2})
                des(val).indvar2   = usrdat.factors{val2};
                des(val).group     = usrdat.factorvals{val2}(valf2);
            end;
            if ~isequal(sort(des(val).condition), sort(usrdat.factorvals{val1}(valf1)))
                des(val).condition = usrdat.factorvals{val1}(valf1);
            end;
            if ~isequal(sort(des(val).group), sort(usrdat.factorvals{val2}(valf2)))
                des(val).group = usrdat.factorvals{val2}(valf2);
            end;
            if ~isequal(sort(des(val).subject), sort(usrdat.subjects(vals)))
                des(val).subject = usrdat.subjects(vals);
            end;
            if ~isequal(des(val).statvar1, valpaired{valp1})
                des(val).statvar1 = valpaired{valp1};
            end;
            if ~isequal(des(val).statvar2, valpaired{valp2})
                des(val).statvar2 = valpaired{valp2};
            end;
            if ~isequal(des(val).deletepreviousfiles, valchk)
                des(val).deletepreviousfiles = 1;
            end;
            if ~isequal(des(val).includevarlist, strs)
                try,
                    des(val).includevarlist = eval( [ '{' strs '}' ]);
                catch,
                    disp('Error while decoding list of parameters');
                    des(val).includevarlist = {};
                    set(findobj(fig, 'tag', 'edit_selectdattrials'),  'string', '');
                end;
            end;
            
        case 'selectfact', % select a specific ind. var. (update value listboxes)
            factval = designind;
            val   = get(findobj(fig, 'tag', 'listboxdesign'), 'value');
            val1  = get(findobj(fig, 'tag', [ 'lbfact' num2str(factval) ]), 'value');
            val2  = get(findobj(fig, 'tag', [ 'lbfact' num2str(~factval) ]), 'value');
            if val1 == val2 && val1 ~= 1
                warndlg2('Cannot select twice the same independent variable');
                val1 = 1;
                set(findobj(fig, 'tag', [ 'lbfact' num2str(factval) ]), 'value', val1);
            end;
            valfact = [1:length(usrdat.factorvals{val1})];
            set(findobj(fig, 'tag', ['lbval' num2str(factval) ]), 'string', usrdat.factorvals{val1}, 'value', valfact);
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
            usrdat.factorvals{val1}{end+1} = usrdat.factorvals{val1}{vals(1)};
            for ind = 2:length(vals)
                usrdat.factorvals{val1}{end} = [ usrdat.factorvals{val1}{end} ' - ' usrdat.factorvals{val1}{vals(ind)} ];
            end;
            set(findobj(fig, 'tag', ['lbval' num2str(factval) ]), 'string', usrdat.factorvals{val1});
            
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
            subjects = unique( { datinfo(indset).subject } );
            
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
                       { 'style' 'listbox' 'string' usrdat.factorvals{1} 'tag' 'lbval2' 'min' 0 'max' 2} };
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
            set(findobj(fig, 'tag', 'lbval2'), 'string', usrdat.factorvals{val1}, 'value', valfact);
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
    [tmp ind] = setdiff(b, a);
    res = setdiff([1:length(b)], ind);
    
