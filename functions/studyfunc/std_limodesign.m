% std_limodesign() - create LIMO design matrices for categorical and 
%                    continuous variables.
%
% Usage:
%   >> [catMat,contMat,description] = std_limosavedesignfiles(factors, trialinfo)
%   >> [catMat,contMat,description] = std_limosavedesignfiles(factors, trialinfo, 'key', val)
%
% Inputs:
%   factors      - [struct] list of factors (structure containing fields
%                  'label' [string], 'value' [string or int] and 'vartype'
%                  ('categorical' or 'continuous'). Structure returned by
%                  function pop_listfactors().
%   trialinfo    - [struct] trial information structure as returned by
%                  std_combtrialinfo().
%
% Optional inputs:
%   filepath     - [string] file path. If not empty, categorical_variables.txt
%                  and continuous_variables.txt files are saved.
%   'splitreg'   - ['on'|'off'] split regression for different categorical
%                  factors. Default is 'off'.
%   'interaction' - ['on'|'off'] compute interaction when using different
%                  categorical variables. This allows computing interactions
%                  between these variables at the second level. Default 
%                  is 'on'.
%   'desconly'    - ['on'|'off'] only output description
%
%  Outputs:
%    catMat      - [array] design matrix for categorical variables
%    contMat     - [array] design matrix for continuous variables
%    description - [cell] what each columns contains
%
% Author: Arnaud Delorme, SCCN, UCSD, 2018-
%
% See also: std_limo()

% Copyright (C) Arnaud Delorme
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

function [catMat,contMat,limodesign] = std_limodesign(factors, trialinfo, varargin)

if nargin < 3
    help std_limosavedesignfiles;
    return; 
end

% decode input parameters
% -----------------------
opt = finputcheck( varargin, ...
    { 'splitreg'       'string'  { 'on','off' }   'off';
      'interaction'    'string'  { 'on','off' }   'on';
      'desconly'       'string'  { 'on','off' }   'off';
      'filepath'       'string'  {            }   ''    }, 'std_limosavedesignfiles');
if ischar(opt), error(opt); end

% number of categorical and continuous var
% ----------------------------------------
catMat  = [];
contMat = [];
catVar  = sort(find(cellfun(@(x)strcmpi(x, 'categorical' ), { factors.vartype })));
contVar = find(cellfun(@(x)strcmpi(x, 'continuous' ), { factors.vartype }));
catVarLabel = unique({ factors(catVar).label });

% find all values for each cat. indep. var.
% -----------------------------------------
alloptions = {};
for iVar = 1:length(catVarLabel)
    indVals = find(cellfun(@(x)strcmpi(x, catVarLabel{iVar} ), { factors.label }));
    values = { factors(indVals).value };
    alloptions{iVar} = cellfun(@(x){catVarLabel{iVar} x}, values, 'uniformoutput', false);
end

% compute interactions if necessary
% ---------------------------------
if length(alloptions) > 1
    alloptionsinter = inter(alloptions);
    alloptionsinter = { alloptionsinter }; % 1 condition only
else
    alloptionsinter = alloptions;
end
limodesign.continuous = {};
if strcmpi(opt.interaction, 'on')
    limodesign.categorical = alloptionsinter;
else
    limodesign.categorical = alloptions;
end

% build design matrix for categorical var.
% ----------------------------------------
if strcmpi(opt.desconly, 'off')
    col = 1;
    catMat = zeros(length(trialinfo),1);
    for iVar = 1:length(limodesign.categorical)
        for iVal = 1:length(limodesign.categorical{iVar})
            trialindsx = std_gettrialsind(trialinfo, limodesign.categorical{iVar}{iVal}{:});
            catMat(trialindsx, iVar) = iVal;
        end
    end
    catMat(catMat == 0) = NaN;
    if ~isempty(opt.filepath)
        save(fullfile(opt.filepath, 'categorical_variables.txt'),'catMat','-ascii');
    end
end

% continuous file/matrix
% ----------------------
if ~isempty(contVar)
    contMat = zeros(length(trialinfo),1);
    
    for iVar = 1:length(contVar)
        limodesign.continuous{iVar} = { factors(contVar(iVar)).label '' };
    end
    
    % split variables if necessary
    if strcmpi(opt.splitreg, 'on')
        limodesign.continuous = inter( { limodesign.continuous, alloptionsinter{:} } );
        for iOpt = 1:length(limodesign.continuous)
            limodesign.continuous{iOpt} = [ limodesign.continuous{iOpt}{:} ];
        end
    end
    
    if strcmpi(opt.desconly, 'off')
        % build continuous design matrix
        for iVar = 1:length(limodesign.continuous)
            [trialindsx, eventvals] = std_gettrialsind(trialinfo, limodesign.continuous{iVar}{:});
            contMat(trialindsx,iVar) = eventvals; % non selected values are assigned 0
        end

        % z-score columns if they are split - same as limo_split_continuous(catMat, contMat)
        if strcmpi(opt.splitreg, 'on')
            for iVar = 1:size(contMat,2)
                nonNan = ~isnan(contMat(:,iVar));
                contMat(nonNan,iVar) = zscore(contMat(nonNan,iVar));
                contMat(~nonNan,iVar) = 0;
            end
        end
        if ~isempty(opt.filepath)
            save(fullfile(opt.filepath, 'continuous_variables.txt'),'contMat','-ascii');
        end
    end
end

% compute interactions in a recursive maner
% ex: res = inter( { { {1 2} {3 4} } { { 5 6 } { 7 8 } { 9 10 } } { { 11 12 } {13 14} } })
% -----------------------------------------
function b = inter(a)

if length(a) == 1
    b = a{1};
else
    c = inter(a(2:end));
    for iVal1 = 1:length(a{1})
        for iVal2 = 1:length(c)
            b{(iVal1-1)*length(c)+iVal2} = { a{1}{iVal1}{:} c{iVal2}{:} };
        end
    end
end

