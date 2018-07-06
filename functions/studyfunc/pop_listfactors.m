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
%                  is 'on'.
% 
% Author: Arnaud Delorme, UCSD, 2018

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, arno@salk.edu
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

function allFactorsStruct = pop_listfactors(des, varargin)

if nargin < 1
    help pop_listfactors;
end 

g = finputcheck(varargin, { 'gui'         'string' { 'on' 'off' } 'on';
                            'splitreg'    'string' { 'on','off' } 'off';
                            'interaction' 'string' { 'on','off' } 'on' });
if isstr(g)
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

% redorders factors so that all variables are grouped

if strcmpi(g.gui, 'on')
    [~,~,des] = std_limodesign(allFactorsStruct,[], 'desconly', 'on', 'splitreg', g.splitreg, 'interaction', g.interaction);
    
    % generate categorical labels
    allLabels = {};
    count     = 1;
    for iCat = 1:length(des.categorival)
        for iVal = 1:length(des.categorival{iCat})
            allLabels{count} = formatcond(des.categorival{iCat}{iVal});
            count = count+1;
        end
    end
    for iCont = 1:length(des.continuous)
        allLabels{count} = formatcond(des.continuous{iCont});
        count = count+1;
    end
    
    % add constant (for GUI)
    allLabels{count} = 'constant';
    
    warndlg2(strvcat(allLabels), 'List of factors');
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
