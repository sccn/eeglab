% std_checkconsist() - Create channel groups for plotting.
%
% Usage:    
%                >> boolval = std_checkconsist(STUDY, 'uniform', 'condition');   
%                >> boolval = std_checkconsist(STUDY, 'uniform', 'group');   
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   'uniform'  - ['condition'|'group'] check if there is one group
%                condition per subject 
% Outputs:
%   boolval    - [0|1] 1 if uniform
%
% Authors: Arnaud Delorme, CERCO, 2009

% Copyright (C) Arnaud Delorme, CERCO, arno@salk.edu
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

function [boolval npersubj] = std_checkconsist(STUDY, varargin);

if nargin < 3
    help std_checkconsist;
    return;
end;

opt = struct(varargin{:});

if strcmpi(opt.uniform, 'condition')
    allvals = { STUDY.datasetinfo.condition };
    vallist = STUDY.condition;
elseif strcmpi(opt.uniform, 'group')
    allvals = { STUDY.datasetinfo.group };
    vallist = STUDY.group;
elseif strcmpi(opt.uniform, 'session')
    allvals = { STUDY.datasetinfo.session };
    allvals = cellfun(@num2str, allvals, 'uniformoutput', false);
    vallist = STUDY.session;
    if isempty(vallist), boolval = 1; return; end;
    vallist = cellfun(@num2str, mattocell(vallist), 'uniformoutput', false);
else
    error('unknown option');
end;

if isempty(vallist), boolval = 1; return; end;

for index = 1:length(vallist)
    tmplist = strmatch( vallist{index}, allvals, 'exact');
    vallen(index) = length(unique( { STUDY.datasetinfo(tmplist).subject } ));
end;
if length(unique(vallen)) == 1
     boolval = 1;
else boolval = 0;
end;

