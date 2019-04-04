% std_selectdataset() - select datasets and trials for a given independent
%                       variable with a given set of values.
%
% Usage:
%   >> [STUDY] = std_selectdataset(STUDY, ALLEEG, indvar, indvarvals);
%
% Inputs:
%   STUDY       - EELAB STUDY structure
%   ALLEEG      - EELAB dataset structure
%   indvar      - [string] independent variable name
%   indvarvals  - [cell] cell array of string for selected values for the 
%   verboseflag - ['verbose'|'silent'] print info flag
%
%                choosen independent variable
% Output:
%   datind       - [integer array] indices of selected dataset
%   dattrialsind - [cell] trial indices for each dataset (not only the
%                  datasets selected above).
%
% Author: Arnaud Delorme, CERCO, 2010-

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function [datind, dattrialselect] = std_selectdataset(STUDY, ALLEEG, indvar, indvarvals, verboseFlag);

if nargin < 3
    help std_selectdataset;
    return;
end
if nargin < 5
    verboseFlag = 'verbose';
end

% check for multiple condition selection
if ~iscell(indvarvals), 
    pos = strfind(' - ', indvarvals);
    if ~isempty(pos)
        tmpindvar = indvarvals;
        indvarvals = { indvarvals(1:pos(1)-1) };
        pos(end+1) = length(tmpindvar)+1;
        for ind = 1:length(pos)-1
            indvarvals{end+1} = tmpindvar(pos(ind)+3:pos(ind+1)-1);
        end
    else
        indvarvals = { indvarvals }; 
    end
end

% default dattrialselect = all trials
% -----------------------------------
if isfield(STUDY.datasetinfo, 'trialinfo')
     dattrialselect = cellfun(@(x)([1:length(x)]), { STUDY.datasetinfo.trialinfo }, 'uniformoutput', false);
else for i=1:length(ALLEEG), dattrialselect{i} = [1:ALLEEG(i).trials]; end
end
    
if isempty(indvar)
    datind = [1:length(STUDY.datasetinfo)];
elseif isfield(STUDY.datasetinfo, indvar) && ~isempty(getfield(STUDY.datasetinfo(1), indvar))
    % regular selection of dataset in datasetinfo
    % -------------------------------------------
    if strcmpi(verboseFlag, 'verbose'), fprintf('   Selecting datasets with field ''%s'' equal to %s\n', indvar, vararg2str(indvarvals)); end
    eval( [ 'myfieldvals = { STUDY.datasetinfo.' indvar '};' ] );
    datind = [];
    for dat = 1:length(indvarvals)
        datind = union_bc(datind, std_indvarmatch(indvarvals{dat}, myfieldvals));
    end
else
    % selection of trials within datasets
    % -----------------------------------
    if strcmpi(verboseFlag, 'verbose'), fprintf('   Selecting trials with field ''%s'' equal to %s\n', indvar, vararg2str(indvarvals)); end
    dattrials      = cellfun(@(x)(eval(['{ x.' indvar '}'])),  { STUDY.datasetinfo.trialinfo }, 'uniformoutput', false);
    dattrials      = cellfun(@(x)(eval(['{ x.' indvar '}'])),  { STUDY.datasetinfo.trialinfo }, 'uniformoutput', false); % do not remove duplicate line (or Matlab crashes)
    dattrialselect = cell(1,length(STUDY.datasetinfo));
    for dat = 1:length(indvarvals)
        for tmpi = 1:length(dattrials)
            dattrialselect{tmpi} = union_bc(dattrialselect{tmpi}, std_indvarmatch(indvarvals{dat}, dattrials{tmpi}));
        end
    end
    datind = find(~cellfun(@isempty, dattrialselect));
end
