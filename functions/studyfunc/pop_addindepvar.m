% Support function for pop_studydesign

% Copyright (C) 2012- Arnaud Delorme
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

function [var values cat] = pop_addindepvar(varlist, fig, var, values)

if nargin < 3
    var = [];
    values = [];
    cat = 0;
end
if isstruct(varlist)
    
    indVar = 1;
    indCat = 1; % categorical
    indVal = 1;
    if ~isempty(var)
        indVar = strmatch(var, varlist.factors, 'exact');
    else
        indVar = 1;
    end
    
    if isempty(values)
        if varlist.numerical(indVar)
            indCat = 2; % continuous
            indVal = [1:length(varlist.factorvals{indVar})];
        end
    else
        for iVal = 1:length(values)
            indVal(iVal) = std_indvarmatch(values{iVal}, varlist.factorvals{indVar});
        end
    end
    
    curValues = encodevals(varlist.factorvals{indVar});
    cb_selectfact  = 'pop_addindepvar(''selectfact'', gcbf);';
    cb_vartype     = 'pop_addindepvar(''vartype''   , gcbf);';
    cb_combinevals = 'pop_addindepvar(''combinevals'', gcbf);';
    uilist = {     { 'style' 'text'       'string' 'Select independent variable' 'fontweight' 'bold' } ...
        { 'style' 'listbox'    'string' varlist.factors  'value' indVar 'tag' 'lbfact0' 'callback' cb_selectfact } ...
        { 'style' 'popupmenu'  'string' 'This is a categorical var.|This is a continuous var.' 'value' indCat 'tag' 'vartype' 'callback' cb_vartype } ...
        { 'style' 'text'       'string' 'Select variable values ' 'tag' 'text1' 'enable' fastif(indCat == 2, 'off', 'on') } ...
        { 'style' 'listbox'    'string' curValues 'value' indVal 'tag' 'lbval0' 'min' 0 'max' 2 'callback' '' 'enable' fastif(indCat == 2, 'off', 'on') } ...
        { 'style' 'pushbutton' 'string' 'Combine selected values' 'tag' 'combine1' 'callback' cb_combinevals 'enable' fastif(indCat == 2, 'off', 'on') } };
    
    w = 1;
    h = 10;
    geometry = { ...
        {w h [1 1] [1 1] } ...
        {w h [1 2] [0.98 3] } ...
        {w h [1 5] [1 1] } ...
        {w h [1 6] [1 1] } ...
        {w h [1 7] [0.98 3] } ...
        {w h [1 10] [1 1] } ...
        };
    
    for i = 1:length(geometry), geometry{i}{3} = geometry{i}{3}-1; end
    streval = [ 'pop_studydesign2(''selectdesign'', gcf);' ];
    [tmp usrdat tmp2 result] = inputgui('uilist', uilist, 'title', 'Add variable', 'geom', geometry, 'userdata', varlist);
    if isempty(tmp), var = []; cat = 0; return, end
        
    var    = usrdat.factors{ result.lbfact0 };
    values = usrdat.factorvals{ result.lbfact0 }(result.lbval0 );
    cat    = result.vartype;
    if cat == 2, cat = 0; end; % continuous
    if ~cat, values = []; end
    
elseif ischar(varlist)
    com = varlist;
    usrdat = get(fig, 'userdata');
    
    switch com
            
        case 'vartype', % select a specific ind. var. (update value listboxes)
            val1  = get(findobj(fig, 'tag', 'vartype'), 'value');
            if val1 == 1, enableFlag = 'on'; else enableFlag = 'off'; end
            set(findobj(fig, 'tag', 'text1'  ),  'enable', enableFlag);  
            set(findobj(fig, 'tag', 'lbval0') ,  'enable', enableFlag);
            set(findobj(fig, 'tag', 'combine1'), 'enable', enableFlag);
            return;
            
        case 'selectfact', % select a specific ind. var. (update value listboxes)
            val1  = get(findobj(fig, 'tag', 'lbfact0'), 'value');
            valfact = [1:length(usrdat.factorvals{val1})];
            set(findobj(fig, 'tag', 'lbval0'), 'string', encodevals(usrdat.factorvals{val1}), 'value', valfact, 'listboxtop', 1);
            set(findobj(fig, 'tag', 'vartype'), 'ListboxTop', 1);
            if usrdat.numerical(val1), 
                set(findobj(fig, 'tag', 'vartype'), 'value', 2, 'string', 'Categorical variable|Continuous variable', 'ListboxTop', length(usrdat.factorvals{val1})); 
            else
                set(findobj(fig, 'tag', 'vartype'), 'value', 1, 'string', 'Categorical variable', 'ListboxTop', length(usrdat.factorvals{val1})); 
            end
            pop_addindepvar('vartype', fig);
                
            return;
            
        case 'combinevals', %  combine values in value listboxes
            val1    = get(findobj(fig, 'tag', 'lbfact0'), 'value');
            vals    = get(findobj(fig, 'tag', 'lbval0' ), 'value');
            strs    = get(findobj(fig, 'tag', 'lbval0' ), 'string');
            if length(vals) == 1
                warndlg2('You need to select several values to combine them');
                return;
            end
            if ~iscell(usrdat.factorvals{val1})
                warndlg2('Cannot combine values from numerical variables');
                return;
            end
            % combine values for string and integers
            if ischar(usrdat.factorvals{val1}{1}) || iscell(usrdat.factorvals{val1}{1})
                tmpcell = {};
                for indCell = vals(:)'
                    if iscell(usrdat.factorvals{val1}{indCell})
                        tmpcell = { tmpcell{:} usrdat.factorvals{val1}{indCell}{:} };
                    else
                        tmpcell = { tmpcell{:} usrdat.factorvals{val1}{indCell} };
                    end
                end
                usrdat.factorvals{val1}{end+1} = unique_bc(tmpcell);
            else
                usrdat.factorvals{val1}{end+1} = unique_bc([ usrdat.factorvals{val1}{vals} ]);
            end
            set(findobj(fig, 'tag', 'lbval0'), 'string', encodevals(usrdat.factorvals{val1}));
            set(fig, 'userdata', usrdat);
            pop_addindepvar('vartype', fig);
                
            return;            
    end
end

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

function [cellout inds ] = mysetdiff(cell1, cell2);
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
                    cellout{index} = [ cellout{index} ' && ' cellin{index}{indcell} ];
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
                    cellout{index} = [ cellout{index} ' && ' num2str(cellin{index}(indcell)) ];
                end
            end
        end
    end
