function [ STUDY ALLEEG ] = std_substudy(STUDY, ALLEEG, varargin);

if nargin < 3
    help std_substudy;
    return;
end

opt = finputcheck(varargin, { 'condition' 'string' {}   '';
                              'group'     'string' {}   '';
                              'subject'   'string' {}   '' }, 'std_substudy');
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

% find new dataset indices
% ------------------------
alldats = [1:length(ALLEEG)];
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

% scan clusters to remove datasets
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

% check cluster consistency
% -------------------------
for index = 1:length(STUDY.cluster)
    STUDY.cluster(index).sets(:) = datcoresp(STUDY.cluster(index).sets(:));
    for i = size(STUDY.cluster(index).sets,2):-1:1
        if all(isnan(STUDY.cluster(index).sets(:,i)))
            STUDY.cluster(index).sets(:,i) = [];
        end;
    end;
end;
STUDY = std_checkset(STUDY, ALLEEG);
