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
if ischar(opt), return; end

% find datasets to remove
% -----------------------
tagdel = [];
if ~isempty(opt.subject)
    for index = 1:length(STUDY.datasetinfo)
        if ~strmatch(STUDY.datasetinfo.subject, opt.subject, 'exact')
            tagdel = [ tagdel index ];
        end
    end
end
if ~isempty(opt.condition)
    for index = 1:length(STUDY.datasetinfo)
        if ~strmatch(STUDY.datasetinfo.condition, opt.condition, 'exact')
            tagdel = [ tagdel index ];
        end
    end
end
if ~isempty(opt.group)
    for index = 1:length(STUDY.datasetinfo)
        if ~strmatch(STUDY.datasetinfo.group, opt.group, 'exact')
            tagdel = [ tagdel index ];
        end
    end
end
if ~isempty(opt.dataset)
    tagdel = [ tagdel setdiff([1:length(ALLEEG)], opt.dataset) ];
end
tagdel = unique_bc(tagdel);

% find new dataset indices
% ------------------------
alldats = [1:length(ALLEEG)];
if strcmpi(opt.rmdat, 'on')
    alldats(tagdel) = [];
    for index = 1:length(ALLEEG)
        tmp = find(alldats == index);
        if isempty(tmp), tmp = NaN; end
        datcoresp(index) = tmp;
    end
    ALLEEG(tagdel) = [];
    STUDY.datasetinfo(tagdel) = [];
    for index = 1:length(STUDY.datasetinfo)
        STUDY.datasetinfo(index).index = index;
    end
else
    alldats(tagdel) = [];
    for index = 1:length(ALLEEG)
        tmp = find(alldats == index);
        if isempty(tmp), tmp = NaN; else tmp = index; end
        datcoresp(index) = tmp;
    end
end

% check channel consistency
% -------------------------
for i = 1:length(STUDY.changrp)
    for c = 1:size(STUDY.changrp(i).setinds,1)
       for g = 1:size(STUDY.changrp(i).setinds,2)
           newinds = datcoresp(STUDY.changrp(i).setinds{c,g});
           nonnans = find(~isnan(newinds));
           STUDY.changrp(i).setinds{c,g} = newinds(nonnans);
           STUDY.changrp(i).allinds{c,g} = STUDY.changrp(i).allinds{c,g}(nonnans);
       end
    end
end

% check cluster consistency
% -------------------------
for index = 1:length(STUDY.cluster)
    STUDY.cluster(index).sets(:) = datcoresp(STUDY.cluster(index).sets(:));
    for i = size(STUDY.cluster(index).sets,2):-1:1
        if all(isnan(STUDY.cluster(index).sets(:,i)))
            STUDY.cluster(index).sets(:,i) = [];
            STUDY.cluster(index).comps(:,i) = [];
        end
    end
    [tmp STUDY.cluster(index).setinds STUDY.cluster(index).allinds] = std_setcomps2cell(STUDY, STUDY.cluster(index).sets, STUDY.cluster(index).comps);
end

STUDY = std_reset(STUDY);
STUDY = std_checkset(STUDY, ALLEEG);
