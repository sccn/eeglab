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

function [datind, dattrialselect] = std_selectdataset(STUDY, ALLEEG, indvar, indvarvals, verboseFlag);

if nargin < 3
    help std_selectdataset;
    return;
end;
if nargin < 5
    verboseFlag = 'verbose';
end;

% check for multiple condition selection
if ~iscell(indvarvals), 
    pos = strfind(' - ', indvarvals);
    if ~isempty(pos)
        tmpindvar = indvarvals;
        indvarvals = { indvarvals(1:pos(1)-1) };
        pos(end+1) = length(tmpindvar)+1;
        for ind = 1:length(pos)-1
            indvarvals{end+1} = tmpindvar(pos(ind)+3:pos(ind+1)-1);
        end;
    else
        indvarvals = { indvarvals }; 
    end;
end;

% default dattrialselect = all trials
% -----------------------------------
if isfield(STUDY.datasetinfo, 'trialinfo')
     dattrialselect = cellfun(@(x)([1:length(x)]), { STUDY.datasetinfo.trialinfo }, 'uniformoutput', false);
else for i=1:length(ALLEEG), dattrialselect{i} = [1:ALLEEG(i).trials]; end;
end;
    
if isempty(indvar)
    datind = [1:length(STUDY.datasetinfo)];
elseif isfield(STUDY.datasetinfo, indvar) && ~isempty(getfield(STUDY.datasetinfo(1), indvar))
    % regular selection of dataset in datasetinfo
    % -------------------------------------------
    if strcmpi(verboseFlag, 'verbose'), fprintf('   Selecting datasets with field ''%s'' equal to %s\n', indvar, vararg2str(indvarvals)); end;
    eval( [ 'myfieldvals = { STUDY.datasetinfo.' indvar '};' ] );
    datind = [];
    for dat = 1:length(indvarvals)
        datind = union_bc(datind, std_indvarmatch(indvarvals{dat}, myfieldvals));
    end;
else
    % selection of trials within datasets
    % -----------------------------------
    if strcmpi(verboseFlag, 'verbose'), fprintf('   Selecting trials with field ''%s'' equal to %s\n', indvar, vararg2str(indvarvals)); end;
    dattrials      = cellfun(@(x)(eval(['{ x.' indvar '}'])),  { STUDY.datasetinfo.trialinfo }, 'uniformoutput', false);
    dattrials      = cellfun(@(x)(eval(['{ x.' indvar '}'])),  { STUDY.datasetinfo.trialinfo }, 'uniformoutput', false); % do not remove duplicate line (or Matlab crashes)
    dattrialselect = cell(1,length(STUDY.datasetinfo));
    for dat = 1:length(indvarvals)
        for tmpi = 1:length(dattrials)
            dattrialselect{tmpi} = union_bc(dattrialselect{tmpi}, std_indvarmatch(indvarvals{dat}, dattrials{tmpi}));
        end;
    end;
    datind = find(~cellfun(@isempty, dattrialselect));
end;
