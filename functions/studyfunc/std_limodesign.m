% std_limodesign() - create LIMO design matrices for categorical and 
%                    continuous variables.
%
% Usage:
%   >> std_limosavedesignfiles(factors, trialinfo)
%   >> std_limosavedesignfiles(factors, trialinfo, 'key', val)
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
      'filepath'       'string'  {            }   ''    }, 'std_limosavedesignfiles');
if isstr(opt), error(opt); end;

% number of categorical and continuous var
% ----------------------------------------
catVar  = sort(find(cellfun(@(x)strcmpi(x, 'categorical' ), { factors.vartype })));
contVar = find(cellfun(@(x)strcmpi(x, 'continuous' ), { factors.vartype }));
catVarLabel = unique({ factors(catVar).label });

% categorical file/matrix
% -----------------------
col = 1;
catMat = zeros(length(trialinfo),1);
for iVar = 1:length(catVarLabel)
    % find all values for this label
    indVals = find(cellfun(@(x)strcmpi(x, catVarLabel{iVar} ), { factors.label }));
    values = { factors(indVals).value };
    limodesign.categorival(iVar).label = catVarLabel{iVar};
    limodesign.categorival(iVar).value = values;
    
    % build design matrix
    for iVal = 1:length(values)
        [trialindsx, eventvals] = std_gettrialsind(trialinfo, catVarLabel{iVar}, values{iVal});
        catMat(trialindsx, iVar) = iVal;
    end
end
catMat(catMat == 0) = NaN;
if ~isempty(opt.filepath)
    save(fullfile(opt.filepath, 'categorical_variables.txt'),'catMat','-ascii');
end

% continuous file/matrix
% ----------------------
contMat = [];
if ~isempty(contVar)
    contMat = zeros(length(trialinfo),1);
    
    % build design matrix
    for iVar = 1:length(contVar)
        limodesign.continuous(iVar).label = factors(contVar(iVar)).label;
        [trialindsx, eventvals] = std_gettrialsind(trialinfo, factors(contVar(iVar)).label);
        contMat(:         ,iVar) = NaN;
        contMat(trialindsx,iVar) = eventvals;
    end
    if strcmpi(opt.splitreg, 'on')
        % splits the continuous design based on categorical var. values
        contMat = limo_split_continuous(contMat, catMat);
    end
    if ~isempty(opt.filepath)
        save(fullfile(filepath, 'continuous_variables.txt'),'contMat','-ascii');
    end
end
