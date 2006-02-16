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

if nargin < 2
    help checkstudy;
    return;
end;
    
modif = 0;
if ~isfield(STUDY, 'name'),  STUDY.name  = ''; modif = 1; end;
if ~isfield(STUDY, 'task'),  STUDY.task  = ''; modif = 1; end;
if ~isfield(STUDY, 'notes'), STUDY.notes = ''; modif = 1; end;
if ~isfield(STUDY, 'filename'),  STUDY.filename  = ''; modif = 1; end;
if ~isfield(STUDY, 'filepath'),  STUDY.filepath  = ''; modif = 1; end;
if ~isfield(STUDY, 'history'),   STUDY.history   = ''; modif = 1; end;
if ~isfield(STUDY, 'subject'),   STUDY.subject   = {}; modif = 1; end;
if ~isfield(STUDY, 'group'),     STUDY.group     = {}; modif = 1; end;
if ~isfield(STUDY, 'session'),   STUDY.session   = {}; modif = 1; end;
if ~isfield(STUDY, 'condition'), STUDY.condition = {}; modif = 1; end;
if ~isfield(STUDY, 'setind'),    STUDY.setind    = []; modif = 1; end;
if ~isfield(STUDY, 'etc'),       STUDY.etc       = []; modif = 1; end;

% all summary fields
% ------------------
try, subject = unique({ STUDY.datasetinfo.subject });
catch, 
     subject = ''; 
     disp('Important warning: not all dataset contain subject code, some functions may crash');
end;
try, group = unique({ STUDY.datasetinfo.group });
catch, 
     group = ''; 
     disp('Important warning: not all dataset contain group info, some functions may crash');
end;
try, condition = unique({ STUDY.datasetinfo.condition });
catch, 
     condition = ''; 
     disp('Important warning: not all dataset contain condition info, some functions may crash');
end;
try, session = unique([STUDY.datasetinfo.session]);
catch, 
     session = ''; 
     disp('Important warning: not all dataset contain integer session info, some functions may crash');
end;
if ~isequal(STUDY.subject,   subject  ), STUDY.subject   = subject;   modif = 1; end;  
if ~isequal(STUDY.group,     group    ), STUDY.group     = group;     modif = 1; end;  
if ~isequal(STUDY.condition, condition), STUDY.condition = condition; modif = 1; end;  
if ~isequal(STUDY.session,   session  ), STUDY.session   = session;   modif = 1; end;  

% recompute setind matrix
% -----------------------
if ~isempty(STUDY.condition)
    if ~isempty(STUDY.session)
        setind = zeros(length(STUDY.condition), length(STUDY.subject) *length(STUDY.session));
    else
        setind = zeros(length(STUDY.condition), length(STUDY.subject) );
    end
else
    if ~isempty(STUDY.session)        
        setind = zeros(1, length(STUDY.subject) *length(STUDY.session));
    else
        setind = zeros(1, length(STUDY.subject) );
    end
end
for k = 1:length(STUDY.datasetinfo)
    setcond = find(strcmp(STUDY.datasetinfo(k).condition, STUDY.condition));
    setsubj = find(strcmp(STUDY.datasetinfo(k).subject,   STUDY.subject));
    setsess = find(strcmp(STUDY.datasetinfo(k).session,   STUDY.session));
    ncomps  = [];
    if ~isempty(setcond)
        if ~isempty(setsess)
            setind(setcond, setsubj * setsess) = k; %A 2D matrix of size [conditions (subjects x sessions)]
        else
            setind(setcond, setsubj ) = k; 
        end
    else
        if ~isempty(setsess)
            setind(1, setsubj * setsess) = k; %A 2D matrix of size [conditions (subjects x sessions)]
        else
            setind(1, setsubj) = k; 
        end
    end
    
    if isfield(STUDY.datasetinfo, 'index')
        if STUDY.datasetinfo(k).index ~= k | isempty(STUDY.datasetinfo(k).index)
            STUDY.datasetinfo(k).index = k; modif = 1; %The dataset index in the current ALLEEG structure
        end;
    else
        STUDY.datasetinfo(k).index = k; modif = 1; %The dataset index in the current ALLEEG structure
    end;
end
% set to NaN empty indices
% ------------------------    
setind( find(setind(:) == 0) ) = NaN;
if ~isequal(setind, STUDY.setind)
    STUDY.setind = setind; modif = 1; 
end;

% set cluster array if empty
% --------------------------
if ~isfield(STUDY, 'cluster'), STUDY.cluster = []; modif = 1; end;
if isempty(STUDY.cluster)
    modif = 1; 
    [STUDY] = cls_createclust(STUDY, ALLEEG, 'ParentCluster');
    STUDY.cluster(1).parent = []; 
    for k = 1:size(STUDY.setind,2)
        ind_nonnan = find(~isnan(STUDY.setind(:,k)));
        ncomps = size(ALLEEG(STUDY.datasetinfo(ind_nonnan(1)).index).icaweights,1);
        STUDY.cluster(1).sets =  [STUDY.cluster(1).sets k*ones(1,ncomps)];
        STUDY.cluster(1).comps = [STUDY.cluster(1).comps      [1:ncomps]];
    end
    if length(STUDY.condition) > 1
        tmp = ones(length(STUDY.condition), length(STUDY.cluster(1).sets));
        for l = 1:length(STUDY.condition)
            tmp(l,:) = STUDY.cluster(1).sets + (l-1)*size(STUDY.setind,2);
        end
        STUDY.cluster(1).sets = tmp;
        clear tmp
    end
else
    for index = 1:length(STUDY.cluster)
        if ~isempty(STUDY.cluster(index).centroid) & ~isstruct(STUDY.cluster(index).centroid)
            STUDY.cluster(index).centroid = [];
            modif = 1; 
        end;
    end;    
end;

% determine if there has been any change
% --------------------------------------
if modif;
    STUDY.saved = 'no';
    eegh('[STUDY, ALLEEG] = checkstudy(STUDYIN, ALLEEG);');
    STUDY.history =  sprintf('%s\n%s',  STUDY.history, '[STUDY, ALLEEG] = checkstudy(STUDYIN, ALLEEG);');
end;
