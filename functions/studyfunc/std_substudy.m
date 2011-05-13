% std_substudy()  - create a sub-STUDY set by removing datasets, conditions, groups, or
%                   subjects.
% Usage:  
%  >> [ STUDY ALLEEG ] = std_substudy(STUDY, ALLEEG, 'key', 'val', ...);
%
% Optional Inputs:
%   STUDY                - existing study structure. 
%   ALLEEG               - vector of EEG dataset structures to be included in the STUDY. 
%
% Optional Inputs:
%   'dataset'   - [integer array] indices of dataset to include in sub-STUDY
%                 Default is all datasets.
%   'subject'   - [cell array] name of subjects to include in sub-STUDY.
%                 Default is all subjects.%
%   'condition' - [cell array] name of conditions to include in sub-STUDY
%                 Default is all conditions.
%   'group'     - [cell array] name of gourps to include in sub-STUDY
%                 Default is all groups.
%   'rmdat'     - ['on'|'off'] actually remove datasets 'on', or simply 
%                 remove all references to these datasets for channels and
%                 clusters ('off').
%
% Example:
%    % create a sub-STUDY using only the first 3 subjects
%    % WARNING: make sure your STUDY is saved before creating a sub-STUDY
%    [STUDY ALLEEG] = std_substudy(STUDY, ALLEEG, 'subject', STUDY.subjects(1:3));
%
% Authors: Arnaud Delorme, CERCO/CNSR & SCCN, INC, UCSD, 2009-

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, arno@sccn.ucsd.edu
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

function [ STUDY ALLEEG ] = std_substudy(STUDY, ALLEEG, varargin);

if nargin < 3
    help std_substudy;
    return;
end

opt = finputcheck(varargin, { 'condition' 'cell' {}      {};
                              'dataset'   'integer' {}   [];
                              'group'     'cell' {}      {};
                              'rmdat'     'string' { 'on','off' }      'on';
                              'subject'   'cell' {}      {} }, 'std_substudy');
if isstr(opt), return; end;

% find datasets to remove
% -----------------------
tagdel = [];
if ~isempty(opt.subject)
    for index = 1:length(STUDY.datasetinfo)
        if ~strmatch(STUDY.datasetinfo.subject, opt.subject, 'exact')
            tagdel = [ tagdel index ];
        end;
    end;
end;
if ~isempty(opt.condition)
    for index = 1:length(STUDY.datasetinfo)
        if ~strmatch(STUDY.datasetinfo.condition, opt.condition, 'exact')
            tagdel = [ tagdel index ];
        end;
    end;
end;
if ~isempty(opt.group)
    for index = 1:length(STUDY.datasetinfo)
        if ~strmatch(STUDY.datasetinfo.group, opt.group, 'exact')
            tagdel = [ tagdel index ];
        end;
    end;
end;
if ~isempty(opt.dataset)
    tagdel = [ tagdel setdiff([1:length(ALLEEG)], opt.dataset) ];
end;
tagdel = unique(tagdel);

% find new dataset indices
% ------------------------
alldats = [1:length(ALLEEG)];
if strcmpi(opt.rmdat, 'on')
    alldats(tagdel) = [];
    for index = 1:length(ALLEEG)
        tmp = find(alldats == index);
        if isempty(tmp), tmp = NaN; end;
        datcoresp(index) = tmp;
    end;
    ALLEEG(tagdel) = [];
    STUDY.datasetinfo(tagdel) = [];
    for index = 1:length(STUDY.datasetinfo)
        STUDY.datasetinfo(index).index = index;
    end;
else
    alldats(tagdel) = [];
    for index = 1:length(ALLEEG)
        tmp = find(alldats == index);
        if isempty(tmp), tmp = NaN; else tmp = index; end;
        datcoresp(index) = tmp;
    end;
end;

% check channel consistency
% -------------------------
for i = 1:length(STUDY.changrp)
    for c = 1:size(STUDY.changrp(i).setinds,1)
       for g = 1:size(STUDY.changrp(i).setinds,2)
           newinds = datcoresp(STUDY.changrp(i).setinds{c,g});
           nonnans = find(~isnan(newinds));
           STUDY.changrp(i).setinds{c,g} = newinds(nonnans);
           STUDY.changrp(i).allinds{c,g} = STUDY.changrp(i).allinds{c,g}(nonnans);
       end;
    end;
end;

% check cluster consistency
% -------------------------
for index = 1:length(STUDY.cluster)
    STUDY.cluster(index).sets(:) = datcoresp(STUDY.cluster(index).sets(:));
    for i = size(STUDY.cluster(index).sets,2):-1:1
        if all(isnan(STUDY.cluster(index).sets(:,i)))
            STUDY.cluster(index).sets(:,i) = [];
            STUDY.cluster(index).comps(:,i) = [];
        end;
    end;
    [tmp STUDY.cluster(index).setinds STUDY.cluster(index).allinds] = std_setcomps2cell(STUDY, STUDY.cluster(index).sets, STUDY.cluster(index).comps);
end;

STUDY = std_reset(STUDY);
STUDY = std_checkset(STUDY, ALLEEG);
