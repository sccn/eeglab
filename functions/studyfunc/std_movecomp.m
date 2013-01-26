% std_movecomp()  - Move ICA component(s) from one cluster to another. 
%
% Usage:    
%       >> [STUDY] = std_movecomp(STUDY, ALLEEG, from_cluster, to_cluster, comps);   
% Inputs:
%   STUDY        - STUDY structure comprising all or some of the EEG datasets in ALLEEG.
%   ALLEEG       - vector of EEG structures in the STUDY, typically created using
%                  load_ALLEEG().  
%   from_cluster - index of the cluster components are to be moved from.  
%   to_cluster   - index of the cluster components are to be moved to.  
%   comps        - [int vector] indices of from_cluster components to move.  
%
% Outputs:
%   STUDY        - input STUDY structure with modified component reassignments.
%
%   Example:
%           >> [STUDY] = std_movecomp(STUDY, ALLEEG, 10, 7, [2 7]);
%           % Move components 2 and 7 of Cluster 10 to Cluster 7. 
%
%  See also:  pop_clustedit         
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, June, 2005

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

function STUDY = std_movecomp(STUDY, ALLEEG, old_clus, new_clus, comps)
icadefs;

% Cannot move components if clusters have children clusters
if ~isempty(STUDY.cluster(old_clus).child)  | ~isempty(STUDY.cluster(new_clus).child)  
    warndlg2('Cannot move components if clusters have children clusters!' , 'Aborting move components');
    return;
end
if isempty(STUDY.cluster(old_clus).parent)  | isempty(STUDY.cluster(new_clus).parent) % The Parent cluster
    warndlg2('Cannot move components to or from the Parent cluster - off all components in STUDY!' , 'Aborting move components');
    return;
end
% Cannot move components if clusters have different parent
% clusters (didn'y come from the same level of clustering), 
% unless the cluster, components are moved to, is an empty new cluster. 
if (length(STUDY.cluster(old_clus).parent) ~= length(STUDY.cluster(new_clus).parent)) & ~strcmp(STUDY.cluster(new_clus).parent, 'manual')
    warndlg2(strvcat('Cannot move components if clusters have different parent clusters!', ...
                     'This limitation will be fixed in the future'), 'Aborting move components');
    return;
end
% Check if all parents are the same
if (~strcmp(STUDY.cluster(new_clus).parent, 'manual')) 
    if ~(sum(strcmp(STUDY.cluster(old_clus).parent, STUDY.cluster(new_clus).parent)) == length(STUDY.cluster(new_clus).parent))% different parent
        warndlg2(strvcat('Cannot move components if clusters have different parent clusters!', ...
                     'This limitation will be fixed in the future') , 'Aborting move components');
        return;
    end
end
for  ci = 1:length(comps)
    comp = STUDY.cluster(old_clus).comps(comps(ci));
    sets = STUDY.cluster(old_clus).sets(:,comps(ci));
    fprintf('Moving component %d from cluster %d to cluster %d, centroids will be recomputed\n',comp, old_clus, new_clus);
    
    %update new cluster
    indcomp = length(STUDY.cluster(new_clus).comps)+1;
    STUDY.cluster(new_clus).comps(indcomp) = comp;%with comp index
    STUDY.cluster(new_clus).sets(:,indcomp) = sets; %with set indices 
    if strcmpi(STUDY.cluster(new_clus).parent, 'manual')
        STUDY.cluster(new_clus).preclust.preclustparams = STUDY.cluster(old_clus).preclust.preclustparams;
        STUDY.cluster(new_clus).parent = STUDY.cluster(old_clus).parent;
        STUDY.cluster(find(strcmp({STUDY.cluster.name},STUDY.cluster(new_clus).parent))).child{end+1} = STUDY.cluster(new_clus).name;
    end
    
    % update preclustering array
    % --------------------------
    if strncmpi('Notclust',STUDY.cluster(old_clus).name,8)
        STUDY.cluster(new_clus).preclust.preclustparams = [];
        STUDY.cluster(new_clus).preclust.preclustdata   = [];
        STUDY.cluster(new_clus).preclust.preclustcomp   = [];
        disp('Important warning: pre-clustering information removed for target cluster');
        disp('(this is because the component moved had no pre-clustering data associated to it)');
    elseif ~strncmpi('Notclust',STUDY.cluster(new_clus).name,8)
        STUDY.cluster(new_clus).preclust.preclustdata(indcomp,:) = STUDY.cluster(old_clus).preclust.preclustdata(comps(ci),:); %with preclustdata
    end;

    % sort by sets
    % ------------
    [tmp,sind] = sort(STUDY.cluster(new_clus).sets(1,:));
    STUDY.cluster(new_clus).sets = STUDY.cluster(new_clus).sets(:,sind);
    STUDY.cluster(new_clus).comps = STUDY.cluster(new_clus).comps(sind);
    if ~isempty(STUDY.cluster(new_clus).preclust.preclustdata)
        STUDY.cluster(new_clus).preclust.preclustdata(sind,:) = STUDY.cluster(new_clus).preclust.preclustdata(:,:);
    end
end
%STUDY.cluster(new_clus).centroid = []; % remove centroid
STUDY = rm_centroid(STUDY, new_clus);
STUDY = rm_centroid(STUDY, old_clus);

% Remove data from old cluster
% left_comps - are all the components of the cluster after the
% components that were moved to the new cluster were removed. 
left_comps = find(~ismember([1:length(STUDY.cluster(old_clus).comps)],comps));
STUDY.cluster(old_clus).comps = STUDY.cluster(old_clus).comps(left_comps);
STUDY.cluster(old_clus).sets = STUDY.cluster(old_clus).sets(:,left_comps);
if ~isempty(STUDY.cluster(old_clus).preclust.preclustdata)
    try,
        STUDY.cluster(old_clus).preclust.preclustdata = STUDY.cluster(old_clus).preclust.preclustdata(left_comps,:);
    catch, % this generates an unknown error but I was not able to reproduce it - AD Sept. 26, 2009        
    end;
end;

% update the component indices
% ----------------------------
STUDY = std_selectdesign(STUDY, ALLEEG, STUDY.currentdesign);
disp('Done.');

% remove cluster information
% --------------------------
function STUDY = rm_centroid(STUDY, clsindex)
    
    keepfields = { 'name' 'parent' 'child' 'comps' 'sets' 'algorithm' 'preclust' };
    allfields  = fieldnames(STUDY.cluster);
    for index = 1:length(allfields)
        if isempty(strmatch(allfields{index}, keepfields))
            STUDY.cluster = setfield( STUDY.cluster, { clsindex }, allfields{index}, []);
        end;
    end;
