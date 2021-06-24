% pop_listfactors() - list independent variables factors for a given design
%
% Usage:
%  >> factors = pop_listfactors(STUDY);
%  >> factors = pop_listfactors(des);
%
% Inputs:
%   STUDY   - existing study structure.
%   des     - existing design
%
% Optional Inputs:
%   'gui'         - ['on'|'off'] pop-up gui to show the list (default:off)
%   'splitreg'   - ['on'|'off'] split regression for different categorical
%                  factors. Default is 'off'.
%   'interaction' - ['on'|'off'] compute interaction when using different
%                  categorical variables. This allows computing interactions
%                  between these variables at the second level. Default
%                  is 'off'.
%   'level'       - ['one'|'two'|'both'] get only first level or second
%                  level factors. Default is 'both'.
%
% Author: Arnaud Delorme, UCSD, 2018

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, arno@salk.edu
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

function allFactorsStruct = pop_listfactors(des, varargin)

if nargin < 1
    help pop_listfactors;
end

g = finputcheck(varargin, { 'gui'         'string' { 'on' 'off' } 'on';
    'splitreg'    'string' { 'on','off' } 'off';
    'contrast'    'string' { 'on','off' } 'off';
    'level'       'string' { 'one','two','both'} 'both';
    'interaction' 'string' { 'on','off' } 'off' });
if ischar(g)
    error(g);
end

if isfield(des,'design')
    des = des.design;
end

allFactors = {}; % (strings) still used to find unique values
allFactorsStruct = [];
count = 1;
for iDes = 1:length(des)
    for iVar = 1:length(des(iDes).variable)
        if strcmpi(des(iDes).variable(iVar).vartype, 'continuous')
            allFactors{count} = sprintf('%s - continuous variable', des(iDes).variable(iVar).label);
            allFactorsStruct(count).vartype = 'continuous';
            allFactorsStruct(count).label = des(iDes).variable(iVar).label;
            allFactorsStruct(count).level = des(iDes).variable(iVar).level;
            count = count+1;
        else
            strVals = getstrval(des(iDes).variable(iVar).value);
            for iSubVal = 1:length(strVals)
                allFactorsStruct(count).vartype = 'categorical';
                allFactorsStruct(count).label = des(iDes).variable(iVar).label;
                allFactorsStruct(count).level = des(iDes).variable(iVar).level;
                allFactorsStruct(count).value = strVals{iSubVal};
                if isnumeric(strVals{iSubVal})
                    allFactors{count} = sprintf('%s - %d', des(iDes).variable(iVar).label, strVals{iSubVal});
                else
                    allFactors{count} = sprintf('%s - %s', des(iDes).variable(iVar).label, strVals{iSubVal});
                end
                count = count+1;
            end
        end
    end
end

% remove duplicates
if length(allFactors) ~= length(unique(allFactors))
    [~, inds ] = unique(allFactors);
    inds = sort(inds);
    allFactors = allFactors(inds);
    allFactorsStruct = allFactorsStruct(inds);
end

% filter first or second level
if ~strcmpi(g.level, 'both')
    allFactorsStruct = filterlevel(allFactorsStruct, g.level);
end

% redorders factors so that all variables are grouped
if strcmpi(g.gui, 'on')
    
    if strcmpi(g.level, 'both')
        allFactorsStruct1 = filterlevel(allFactorsStruct, 'one');
        allFactorsStruct2 = filterlevel(allFactorsStruct, 'two');
        [~,~,des1] = std_limodesign(allFactorsStruct1,[], 'desconly', 'on', 'splitreg', g.splitreg, 'interaction', g.interaction);
        [~,~,des2] = std_limodesign(allFactorsStruct2,[], 'desconly', 'on', 'splitreg', g.splitreg, 'interaction', g.interaction);
    elseif strcmpi(g.level, 'one')
        [~,~,des1] = std_limodesign(allFactorsStruct,[], 'desconly', 'on', 'splitreg', g.splitreg, 'interaction', g.interaction);
        des2 = {};
    elseif strcmpi(g.level, 'two')
        [~,~,des2] = std_limodesign(allFactorsStruct,[], 'desconly', 'on', 'splitreg', g.splitreg, 'interaction', g.interaction);
        des1 = {};
    end
    
    % generate categorical labels
    allLabels1 = getlabels(des1);
    allLabels2 = getlabels(des2);
    
    listui = {};
    if ~isempty(des2)
        listui{2,1} = { 'Style', 'text', 'string' '2nd-level variables' 'fontweight' 'bold' };
        listui{2,2} = { 'Style', 'text', 'string' '(inter participant)' };
        for index = 1:length(allLabels2)
            listui{2,index+2} = { 'Style', 'text', 'string' allLabels2{index} };
        end
    end
    if ~isempty(des1)
        listui{1,1} = { 'Style', 'text', 'string' '1st-level variables' 'fontweight' 'bold' };
        listui{1,2} = { 'Style', 'text', 'string' '(intra participant)' };
        for index = 1:length(allLabels1)
            listui{1,index+2} = { 'Style', 'text', 'string' allLabels1{index} };
        end
    else
        listui(1,:) = [];
    end
    
    if strcmpi(g.contrast, 'on')
        listui{2,1} = { };
        listui{2,2} = { 'Style', 'text', 'string' 'Weight' };
        for index = 1:max(length(allLabels2), length(allLabels1))
            listui{2,index+2} = { 'Style', 'edit', 'string' '' };
        end
        geometry = cell(1,size(listui,2));
        geometry(:) = { [ 1 0.3 ] };
        listui = listui(:)';
        warningmsg = [ 'warndlg2([ ''A contrast is linear combination of variables (model parameters).'' 10 ' ...
            '''Two types of contrasts are valid:'' 10 ' ...
            ''' '' 10 ' ...
            '''- sum/averages: add together variables, e.g., weights of 0.5 (1st var),'' 10 ' ...
            '''and 0.5 (2nd var) creates the average of variables 1 and 2.'' 10 ' ...
            ''' '' 10 ' ...
            '''- differences: weight variables with coefficients adding up to zero, '' 10 ' ...
            '''e.g. 0.5, 0.5, -0.5 and -0.5 creates the difference between the averages '' 10 ' ...
            '''of variables 1&2 vs. 3&4. '' 10 ' ...
            ''' '' 10 ' ...
            '''Make sure you create meaningful contrasts that are easy to interpret. '' 10 ' ...
            '''Plotting and checking the values of an ERP/Spectra/ERSP contrast '' 10 ' ...
            '''[0.5 0.5 -0.5 -0.5] is looking at the difference between two means. '' 10 ' ...
            '''Plotting and checking the values of an ERP/Spectra/ERSP contrast '' 10 ' ...
            '''[1 1 -1 -1] is looking at the difference between two sums. These two '' 10 ' ...
            '''contrasts will give the same statistical result (T/p values) but the '' 10 ' ...
            '''former is easier to interpret than the latter. Contrast validity is '' 10 ' ...
            '''automatically checked.'' ]);' ];
        
        
        [~, ~, allobj] = inputgui( 'geometry', geometry, 'uilist', listui, 'helpcom',  warningmsg);
    else
        geometry = cell(1,size(listui,2));
        geometry(:) = {ones(1,size(listui,1))};
        listui = listui(:);
        listui{end+1} = {};
        geometry{end+1} = [1];
        geometry = { geometry{:} 1 };
        listui = {listui{:} { 'width',80,'align','center','Style', 'pushbutton', 'string', 'OK', 'callback', ['set(gcbf, ''userdata'', ''OK'');'] }  };
        
        fig = figure('visible', 'off');
        [~, ~, allobj] = supergui( 'fig', fig, 'geomhoriz', geometry, 'uilist', listui, ...
            'borders', [0.05 0.015 0.08 0.06], 'spacing', [0 0], 'horizontalalignment', 'left', 'adjustbuttonwidth', 'off' );
        
        waitfor( fig, 'userdata');
    end
    
    try,
        result = get(fig, 'userdata');
        close(fig);
        drawnow;
    end
    
end

% convert nested values to linear sequence
function res = getstrval(vals)

res = {};
if iscell(vals) && length(vals) == 1
    vals = vals{1};
end
if iscell(vals) || (isnumeric(vals) && length(vals) > 1)
    for iVal = 1:length(vals)
        restmp = getstrval(vals(iVal));
        res = { res{:} restmp{:} };
    end
else
    res =  { vals };
end

% format string
function str = formatcond(cellVal)

for iItem = 1:2:length(cellVal)
    if iItem+1>length(cellVal) || isempty(cellVal{iItem+1}) % continuous var
        tmpFactor = sprintf('%s (continuous)', cellVal{iItem});
    elseif isnumeric(cellVal{iItem+1})
        tmpFactor = sprintf('%s = %d', cellVal{iItem}, cellVal{iItem+1});
    else tmpFactor = sprintf('%s = %s', cellVal{iItem}, cellVal{iItem+1});
    end
    if iItem == 1
        str = tmpFactor;
    else str = sprintf('%s & %s', str, tmpFactor);
    end
end

% select variables at a specific level
function allFactorsStruct = filterlevel(allFactorsStruct, level)
rmInd = [];
for ind = 1:length(allFactorsStruct)
    if ~strcmpi(allFactorsStruct(ind).level,level)
        rmInd = [ rmInd ind ];
    end
end
allFactorsStruct(rmInd) = [];

% get labels for design
function allLabels = getlabels(des)
allLabels = {};
count     = 1;
if isempty(des), return; end
for iCat = 1:length(des.categorical)
    for iVal = 1:length(des.categorical{iCat})
        allLabels{count} = [ int2str(count) '. ' formatcond(des.categorical{iCat}{iVal}) ];
        count = count+1;
    end
end
for iCont = 1:length(des.continuous)
    if iscell(des.continuous{iCont})
        allLabels{count} = [ int2str(count) '. ' formatcond(des.continuous{iCont}) ];
    else
        allLabels{count} = [ int2str(count) '. ' formatcond({ des.continuous{iCont} } ) ];
    end
    count = count+1;
end
allLabels{count} = [ int2str(count) '. Constant' ];
