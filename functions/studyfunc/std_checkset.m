% std_checkset() - check STUDY set consistency. Make changes to STUDY and issue
%                  commandline warnings as needed.
%
% Usage: >> [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG);  
%
% Inputs:
%   STUDY   - EEGLAB STUDY set
%   ALLEEG  - vector of EEG datasets included in the STUDY structure 
%
% Outputs:
%   STUDY   - a possibly modified STUDY set 
%   ALLEEG  - the unmodified ALLEEG input
%
% Authors:  Arnaud Delorme & Hilit Serby, SCCN, INC, UCSD, November 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme & Scott Makeig, SCCN/INC/UCSD, smakeig@ucsd.edu
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

function [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG);

if nargin < 2
    help std_checkset;
    return;
end;
modified = 0; % 0 = NO modifications made; 1 = modified
    
% check existence of all STUDY fields
% ------------------------------------
if ~isfield(STUDY, 'name'),  STUDY.name  = ''; modified = 1; end;
if ~isfield(STUDY, 'task'),  STUDY.task  = ''; modified = 1; end;
if ~isfield(STUDY, 'notes'), STUDY.notes = ''; modified = 1; end;
if ~isfield(STUDY, 'filename'),  STUDY.filename  = ''; modified = 1; end;
if ~isfield(STUDY, 'filepath'),  STUDY.filepath  = ''; modified = 1; end;
if ~isfield(STUDY, 'history'),   STUDY.history   = ''; modified = 1; end;
if ~isfield(STUDY, 'subject'),   STUDY.subject   = {}; modified = 1; end;
if ~isfield(STUDY, 'group'),     STUDY.group     = {}; modified = 1; end;
if ~isfield(STUDY, 'session'),   STUDY.session   = {}; modified = 1; end;
if ~isfield(STUDY, 'condition'), STUDY.condition = {}; modified = 1; end;
if ~isfield(STUDY, 'setind'),    STUDY.setind    = []; modified = 1; end;
if ~isfield(STUDY, 'etc'),       STUDY.etc       = []; modified = 1; end;
if ~isfield(STUDY, 'datasetinfo'), STUDY.datasetinfo = []; modified = 1; end;
if ~isfield(STUDY.datasetinfo, 'comps') & ~isempty(STUDY.datasetinfo), STUDY.datasetinfo(1).comps = []; modified = 1; end;
if ~isfield(STUDY.datasetinfo, 'index') & ~isempty(STUDY.datasetinfo), STUDY.datasetinfo(1).index = []; modified = 1; end;

% check all summary fields
% -------------------------
try, subject = unique({ STUDY.datasetinfo.subject });
catch, 
     subject = ''; 
     disp('Important warning: Some datasets do not have subject codes; some functions may crash!');
end;
try, group = unique({ STUDY.datasetinfo.group });
catch, 
     group = ''; 
     % disp('Important warning: Some datasets do not have group codes; some functions may crash!');
end;
try, condition = unique({ STUDY.datasetinfo.condition });
catch, 
     condition = ''; 
     disp('Important warning: Some datasets do not have condition codes; some functions may crash!');
end;
try, session = unique([STUDY.datasetinfo.session]);
catch, 
     session = ''; 
     % disp('Important warning: Some datasets do not have session numbers; some functions may crash!');
end;
if ~isequal(STUDY.subject,   subject  ), STUDY.subject   = subject;   modified = 1; end;  
if ~isequal(STUDY.group,     group    ), STUDY.group     = group;     modified = 1; end;  
if ~isequal(STUDY.condition, condition), STUDY.condition = condition; modified = 1; end;  
if ~isequal(STUDY.session,   session  ), STUDY.session   = session;   modified = 1; end;  

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
    setsess = find(STUDY.datasetinfo(k).session == STUDY.session);
    ncomps  = [];
    if ~isempty(setcond)
        if ~isempty(setsess)
            setind(setcond, setsubj * length(STUDY.session)+setsess-1) = k; 
            %A 2D matrix of size [conditions (subjects x sessions)]
        else
            setind(setcond, setsubj ) = k; 
        end
    else
        if ~isempty(setsess)
            setind(1, setsubj * length(STUDY.session)+setsess-1) = k; 
        else
            setind(1, setsubj) = k; 
        end
    end
    
    if ~isequal(STUDY.datasetinfo(k).index, k)
        STUDY.datasetinfo(k).index = k; modified = 1; %The dataset index in the current ALLEEG structure
    end;
end

% set to NaN empty indices and remove NaN columns
% -----------------------------------------------
setind( find(setind(:) == 0) ) = NaN;
rmind = [];
for k = 1:size(setind,2)
    ind_nonnan = find(~isnan(setind(:,k)));
    if isempty(ind_nonnan), rmind = [ rmind k ]; end;
end;
setind(:,rmind) = [];
if ~isequal(setind, STUDY.setind)
    STUDY.setind = setind; modified = 1;
end;

% set cluster array if empty
% --------------------------
if ~isfield(STUDY, 'cluster'), STUDY.cluster = []; modified = 1; end;
if isempty(STUDY.cluster)
    modified = 1; 
    [STUDY] = std_createclust(STUDY, ALLEEG, 'ParentCluster');
    STUDY.cluster(1).parent = []; 
    for k = 1:size(STUDY.setind,2)
        
        ind_nonnan = find(~isnan(STUDY.setind(:,k)));
        ind_nonnan = STUDY.setind(ind_nonnan(1),k);
        comps = STUDY.datasetinfo(ind_nonnan).comps;
        if isempty(comps)
            comps = 1:size(ALLEEG(STUDY.datasetinfo(ind_nonnan).index).icaweights,1);
        end;
        STUDY.cluster(1).sets =  [STUDY.cluster(1).sets       STUDY.setind(:,k)*ones(1,length(comps))];
        STUDY.cluster(1).comps = [STUDY.cluster(1).comps      comps];
    end
else
    for index = 1:length(STUDY.cluster)
        if ~isempty(STUDY.cluster(index).centroid) & ~isstruct(STUDY.cluster(index).centroid)
            STUDY.cluster(index).centroid = [];
            modified = 1; 
        end;
    end;    
end;

% determine if there has been any change
% --------------------------------------
if modified;
    STUDY.saved = 'no'; % mark as "changed without saving to disk"
    eegh('[STUDY, ALLEEG] = std_checkset(STUDYIN, ALLEEG);'); % add to EEGLAB session history
    STUDY.history =  sprintf('%s\n%s',  STUDY.history, '[STUDY, ALLEEG] = std_checkset(STUDYIN, ALLEEG);');
end;
