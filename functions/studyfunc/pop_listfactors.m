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
    rmInd = [];
    for ind = 1:length(allFactorsStruct)
        if ~strcmpi(allFactorsStruct(ind).level,g.level)
            rmInd = [ rmInd ind ];
        end
    end
    allFactorsStruct(rmInd) = [];
end

% redorders factors so that all variables are grouped
if strcmpi(g.gui, 'on')
    [~,~,des] = std_limodesign(allFactorsStruct,[], 'desconly', 'on', 'splitreg', g.splitreg, 'interaction', g.interaction);
    if ~isfield(des, 'categorical'), des.categorical = {}; end
    if ~isfield(des, 'continuous'),  des.continuous  = {}; end
    
    % generate categorical labels
    allLabels = {};
    count     = 1;
    for iCat = 1:length(des.categorical)
        for iVal = 1:length(des.categorical{iCat})
            allLabels{count} = [ int2str(count) '. ' formatcond(des.categorical{iCat}{iVal}) ];
            count = count+1;
        end
    end
    for iCont = 1:length(des.continuous)
        allLabels{count} = [ int2str(count) '. ' formatcond(des.continuous{iCont}) ];
        count = count+1;
    end
    
    % add constant (for GUI)
    allLabels{count} = [ int2str(count) '. Constant' ];
    
    warndlg2(strvcat(allLabels), 'List of explanatory variables');
    
%     Prompt = strvcat(allLabels);
%     listui = {};
%     geometry = {};
%     for index = 1:size(Prompt,1)
%         geometry{index} = [1];
%         listui{index} = { 'Style', 'text', 'string' Prompt(index,:) };
%     end
%     listui{end+1} = {};
%     
%     geometry = { geometry{:} 1 ones(1,length(varargin)-1) };
%     for index = 1:length(varargin)-1 % ignoring default val
%         listui = {listui{:} { 'width',80,'align','center','Style', 'pushbutton', 'string', varargin{index}, 'callback', ['set(gcbf, ''userdata'', ''' varargin{index} ''');'] }  };
%         if strcmp(varargin{index}, varargin{end})
%             listui{end}{end+1} = 'fontweight';
%             listui{end}{end+1} = 'bold';
%         end
%     end
%     
%     fig = figure('visible', 'off');
%     [~, ~, allobj] = supergui( 'fig', fig, 'geomhoriz', geometry, 'uilist', listui, ...
%             'borders', [0.05 0.015 0.08 0.06], 'spacing', [0 0], 'horizontalalignment', 'left', 'adjustbuttonwidth', 'off' );
%     
%     waitfor( fig, 'userdata');
%     
%     try,
%         result = get(fig, 'userdata');
%         close(fig);
%         drawnow;
%     end
    
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
    if isempty(cellVal{iItem+1}) % continuous var
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
