% checkstudy() - check study consistency
%
% Usage: >> [STUDY, ALLEEG] = checkstudy(STUDY, ALLEEG);  
%
% Input:
%   STUDY      - STUDY set
%   ALLEEG     - EEGLAB vector of EEG sets included in the STUDY structure 
%
% Output:
%   STUDY      - a new STUDY set containing some or all of the datasets in ALLEEG, 
%                plus additional information from the optional inputs above. 
%   ALLEEG     - an EEGLAB vector of EEG sets included in the STUDY structure 
%
% Authors:  Arnaud Delorme & Hilit Serby, SCCN, INC, UCSD, November 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

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

function [STUDY, ALLEEG] = checkstudy(STUDY, ALLEEG);

if ~isfield(STUDY, 'name'),  STUDY.name  = ''; end;
if ~isfield(STUDY, 'task'),  STUDY.task  = ''; end;
if ~isfield(STUDY, 'notes'), STUDY.notes = ''; end;
if ~isfield(STUDY, 'filename'),  STUDY.filename  = ''; end;
if ~isfield(STUDY, 'filepath'),  STUDY.filepath  = ''; end;
if ~isfield(STUDY, 'history'),   STUDY.history   = ''; end;
if ~isfield(STUDY, 'etc'),       STUDY.etc       = []; end;

% all summary fields
% ------------------
try, STUDY.subject = unique({ STUDY.datasetinfo.subject });
catch, 
     STUDY.subject = ''; 
     disp('Important warning: not all dataset contain subject code, some functions may crash');
end;
try, STUDY.group = unique({ STUDY.datasetinfo.group });
catch, 
     STUDY.group = ''; 
     disp('Important warning: not all dataset contain group info, some functions may crash');
end;
try, STUDY.condition = unique({ STUDY.datasetinfo.condition });
catch, 
     STUDY.condition = ''; 
     disp('Important warning: not all dataset contain condition info, some functions may crash');
end;
try, STUDY.session = unique([STUDY.datasetinfo.session]);
catch, 
     STUDY.session = ''; 
     disp('Important warning: not all dataset contain integer session info, some functions may crash');
end;

% recompute setind matrix
% -----------------------
if ~isempty(STUDY.condition)
    if ~isempty(STUDY.session)
        STUDY.setind = zeros(length(STUDY.condition), length(STUDY.subject) *length(STUDY.session));
    else
        STUDY.setind = zeros(length(STUDY.condition), length(STUDY.subject) );
    end
else
    if ~isempty(STUDY.session)        
        STUDY.setind = zeros(1, length(STUDY.subject) *length(STUDY.session));
    else
        STUDY.setind = zeros(1, length(STUDY.subject) );
    end
end
for k = 1:length(STUDY.datasetinfo)
    setcond = find(strcmp(STUDY.datasetinfo(k).condition, STUDY.condition));
    setsubj = find(strcmp(STUDY.datasetinfo(k).subject,   STUDY.subject));
    setsess = find(strcmp(STUDY.datasetinfo(k).session,   STUDY.session));
    ncomps  = [];
    if ~isempty(setcond)
        if ~isempty(setsess)
            STUDY.setind(setcond, setsubj * setsess) = k; %A 2D matrix of size [conditions (subjects x sessions)]
        else
            STUDY.setind(setcond, setsubj ) = k; 
        end
    else
        if ~isempty(setsess)
            STUDY.setind(1, setsubj * setsess) = k; %A 2D matrix of size [conditions (subjects x sessions)]
        else
            STUDY.setind(1, setsubj) = k; 
        end
    end    
    STUDY.datasetinfo(k).index = k; %The dataset index in the current ALLEEG structure
    STUDY.setind( find(STUDY.setind(:) == 0) ) = NaN;
end

% number of component per dataset
% -------------------------------
for k = 1:length(STUDY.datasetinfo)
   STUDY.datasetinfo(k).ncomps = size(ALLEEG(STUDY.datasetinfo(k).index).icawinv,2);
end

% set cluster array if empty
% --------------------------
if ~isfield(STUDY, 'cluster'), STUDY.cluster = []; end;
if isempty(STUDY.cluster)
    [STUDY] = cls_createclust(STUDY, ALLEEG, 'ParentCluster');
    STUDY.cluster(1).parent = []; 
    for k = 1:size(STUDY.setind,2)
        ind_nonnan = find(~isnan(STUDY.setind(:,k)));
        STUDY.cluster(1).sets =  [STUDY.cluster(1).sets k*ones(1,STUDY.datasetinfo(ind_nonnan(1)).ncomps)];
        STUDY.cluster(1).comps = [STUDY.cluster(1).comps      [1:STUDY.datasetinfo(ind_nonnan(1)).ncomps]];
    end
    if length(STUDY.condition) > 1
        tmp = ones(length(STUDY.condition), length(STUDY.cluster(1).sets));
        for l = 1:length(STUDY.condition)
            tmp(l,:) = STUDY.cluster(1).sets + (l-1)*size(STUDY.setind,2);
        end
        STUDY.cluster(1).sets = tmp;
        clear tmp
    end
end;
