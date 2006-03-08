%  std_createclust()  - Commandline function, to create a new empty cluster.
%                            After the empty cluster is created, components can be
%                            reassigned to it using the commandline std_movecomp().
% Usage:
%                   >> [STUDY] = std_createclust(STUDY, ALLEEG, name);
% Inputs:
%   STUDY         - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG        - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY.
%                       ALLEEG for a STUDY set is typically created using load_ALLEEG().
% Optional inputs:
%   name          - a name string for the new cluster. {default: 'Cls #', where '#' is
%                     the next available cluster number}.
% Outputs:
%   STUDY    - the input STUDY set structure modified with the new cluster.
%
%   Example:
%                           >> name = 'eye_movements';
%                         >> [STUDY] = std_createclust(STUDY, ALLEEG, name);
%                    A new cluster is created named 'Ceye_movements'.
%
%  See also  pop_clustedit, std_movecomp
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, June, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, June 07, 2005, hilit@sccn.ucsd.edu
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

function [STUDY, clusters] = std_createclust(STUDY, ALLEEG, varargin)

    if length(varargin) > 1
        [STUDY, clusters] = std_createclust2(STUDY, ALLEEG, varargin{:});
        return;
    end;
    
 if isempty(varargin) | strcmpi(varargin,'')
    name = 'Cls';
else
    name = varargin{1};
end
% Find out the highst cluster id number (in cluster name), to find
% next available cluster index
max_id = 0;
for k = 1:length(STUDY.cluster)
    ti = strfind(STUDY.cluster(k).name, ' ');
    clus_id = STUDY.cluster(k).name(ti(end) + 1:end);
    max_id = max(max_id, str2num(clus_id));
end
max_id = max_id + 1;
name = sprintf('%s %d', name, max_id);
STUDY.cluster(end+1).name = name;
% Initialize the new cluster fields.
STUDY.cluster(end).parent{1} = 'manual'; % update parent cluster if exists.
STUDY.cluster(end).child = [];
STUDY.cluster(end).comps = [];
STUDY.cluster(end).sets = [];
STUDY.cluster(end).algorithm = [];
STUDY.cluster(end).centroid = [];
STUDY.cluster(end).preclust.preclustparams = [];
STUDY.cluster(end).preclust.preclustdata = [];


function [STUDY, clusters] = std_createclust2(STUDY,IDX,C, algorithm)

clusters = [];
sets = [];
comp = [];
Ncond = length(STUDY.condition);
if Ncond == 0
    Ncond = 1;
end
nsets = length(STUDY.etc.preclust.preclustcomps); 
for k = 1: nsets
    sets = [sets k*ones(1,length(STUDY.etc.preclust.preclustcomps{k}))];
    comp = [comp STUDY.etc.preclust.preclustcomps{k}];
end
if Ncond > 1
    tmp = ones(Ncond, length(sets));
    for l = 1:Ncond
        tmp(l,:) = sets + (l-1)*nsets;
    end
    sets = tmp;
    clear tmp
end

cls = size(C,1);
% Find the next available cluster index
nc = 0; 
for k =  1:length(STUDY.cluster)
    ti = strfind(STUDY.cluster(k).name, ' ');
    tmp = STUDY.cluster(k).name(ti(end) + 1:end);
    nc = max(nc,str2num(tmp));
    % check if there is a cluster of Notclust components
    if strcmp(STUDY.cluster(k).parent,STUDY.cluster(STUDY.etc.preclust.clustlevel).name) 
        STUDY.cluster(k).preclust.preclustparams = STUDY.etc.preclust.preclustparams;
        clusters = [clusters k];
    end
end
len = length(STUDY.cluster);

if ~isempty(find(IDX==0)) %outliers exist
    l = 0;
    nc = nc + 1;
    len = len +1;
else
    l = 1;
end

for k = l:cls 
    tmp = find(IDX==k);
    STUDY.cluster(k+len).sets = sets(:,tmp);
    STUDY.cluster(k+len).comps = comp(tmp);
    if k == 0
        STUDY.cluster(len).name = ['outlier ' num2str(k+nc)];
    else
        STUDY.cluster(k+len).name = [ 'Cls ' num2str(k+nc)];
    end
    STUDY.cluster(k+len).algorithm = algorithm;
    if STUDY.etc.preclust.clustlevel == 0
        STUDY.cluster(k+len).parent = [];
    else
        STUDY.cluster(k+len).parent{end+1} = STUDY.cluster(STUDY.etc.preclust.clustlevel).name;
        %update parents clusters with cluster child indices
        STUDY.cluster(STUDY.etc.preclust.clustlevel).child{end+1} = STUDY.cluster(k+nc).name;
    end
    STUDY.cluster(k+len).child = [];
    STUDY.cluster(k+len).preclust.preclustdata = STUDY.etc.preclust.preclustdata(tmp,:);
    STUDY.cluster(k+len).preclust.preclustparams = STUDY.etc.preclust.preclustparams;
    STUDY.cluster(k+len).preclust.preclustcomps = STUDY.etc.preclust.preclustcomps;
end

clusters = [ clusters l+len:cls+len];%the new created clusters indices.
