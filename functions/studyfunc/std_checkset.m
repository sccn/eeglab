% std_checkset() - check STUDY set consistency
%
% Usage: >> [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG);  
%
% Input:
%   STUDY      - EEGLAB STUDY set
%   ALLEEG     - vector of EEG datasets included in the STUDY structure 
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

function [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG, option);

if nargin < 2
    help std_checkset;
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
if ~isfield(STUDY, 'preclust'),  STUDY.preclust  = []; modif = 1; end;
if ~isfield(STUDY, 'datasetinfo'), STUDY.datasetinfo = []; modif = 1; end;
if ~isfield(STUDY.preclust, 'erpclusttimes' ),  STUDY.preclust.erpclusttimes = []; modif = 1; end;
if ~isfield(STUDY.preclust, 'specclustfreqs' ), STUDY.preclust.specclustfreqs = []; modif = 1; end;
if ~isfield(STUDY.preclust, 'erspclustfreqs' ), STUDY.preclust.erspclustfreqs = []; modif = 1; end;
if ~isfield(STUDY.preclust, 'erspclusttimes' ), STUDY.preclust.erspclusttimes = []; modif = 1; end;
if ~isfield(STUDY.datasetinfo, 'comps') & ~isempty(STUDY.datasetinfo), STUDY.datasetinfo(1).comps = []; modif = 1; end;
if ~isfield(STUDY.datasetinfo, 'index') & ~isempty(STUDY.datasetinfo), STUDY.datasetinfo(1).index = []; modif = 1; end;

% all summary fields
% ------------------
try, subject = unique({ STUDY.datasetinfo.subject });
catch, 
     subject = ''; 
     disp('Important warning: some datasets do not have subject codes; some functions may crash!');
end;
try, group = unique({ STUDY.datasetinfo.group });
catch, 
     group = {}; 
     % disp('Important warning: some datasets do not have group codes; some functions may crash!');
end;
try, condition = unique({ STUDY.datasetinfo.condition });
catch, 
     condition = {}; 
     disp('Important warning: some datasets do not have condition codes; some functions may crash!');
end;
try, session = unique([STUDY.datasetinfo.session]);
catch, 
     session = ''; 
     % disp('Important warning: some datasets do not have session numbers; some functions may crash!');
end;
if ~isequal(STUDY.subject,   subject  ), STUDY.subject   = subject;   modif = 1; end;  
if ~isequal(STUDY.group,     group    ), STUDY.group     = group;     modif = 1; end;  
if ~isequal(STUDY.condition, condition), STUDY.condition = condition; modif = 1; end;  
if ~isequal(STUDY.session,   session  ), STUDY.session   = session;   modif = 1; end;  

% recompute setind matrix
% -----------------------
notsameica = [];
for is = 1:length(STUDY.subject)
    alldats = strmatch(STUDY.subject{is}, { STUDY.datasetinfo.subject });

    for ig = 1:length(STUDY.group)
        tmpind  = strmatch(STUDY.group{ig}, { STUDY.datasetinfo(alldats).group });
        tmpdats = alldats(tmpind);

        nc = size(ALLEEG(STUDY.datasetinfo(tmpdats(1)).index).icaweights,1);
        for ir = 2:length(tmpdats)
            if nc ~= size(ALLEEG(STUDY.datasetinfo(tmpdats(ir)).index).icaweights,1)
                notsameica = [ 1 tmpdats(1) tmpdats(ir) ];
            end;
        end;
    end;
end;
if ~isempty(notsameica)
    %disp('Different ICA decompositions have been found for the same')
    %disp('subject in two conditions (if the data was recorded at the')
    %disp('time, it is best to have the run ICA on both datasets
    %simultanously.')
    setind = [1:length(STUDY.datasetinfo)];
    if ~isequal(STUDY.setind, setind)
        STUDY.setind = setind; modif = 1;
    end;
else
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
            STUDY.datasetinfo(k).index = k; modif = 1; %The dataset index in the current ALLEEG structure
        end;
    end
end;

% check dataset info consistency 
% ------------------------------
for k = 1:length(STUDY.datasetinfo)
    if ~strcmpi(STUDY.datasetinfo(k).filename, ALLEEG(k).filename)
        STUDY.datasetinfo(k).filename = ALLEEG(k).filename; modif = 1;
        fprintf('Warning: file name has changed for dataset %d and the study has been updated\n', k);
        fprintf('         to discard this change in the study, reload it from disk\n');
    end;
end;

% set to NaN empty indices and remove nan columns
% -----------------------------------------------
setind( find(setind(:) == 0) ) = NaN;
rmind = [];
for k = 1:size(setind,2)
    ind_nonnan = find(~isnan(setind(:,k)));
    if isempty(ind_nonnan), rmind = [ rmind k ]; end;
end;
setind(:,rmind) = [];
if ~isequal(setind, STUDY.setind)
    STUDY.setind = setind; modif = 1;
end;

% set cluster array if empty
% --------------------------
if ~isfield(STUDY, 'cluster'), STUDY.cluster = []; modif = 1; end;
if isempty(STUDY.cluster)
    modif = 1; 
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
            modif = 1; 
        end;
    end;    
end;

% determine if there has been any change
% --------------------------------------
if modif;
    STUDY.saved = 'no';
    eegh('[STUDY, ALLEEG] = std_checkset(STUDYIN, ALLEEG);');
    STUDY.history =  sprintf('%s\n%s',  STUDY.history, '[STUDY, ALLEEG] = std_checkset(STUDYIN, ALLEEG);');
end;
