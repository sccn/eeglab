%  cls_movecomp()  - Commandline function, to move component(s) from one cluster to another. 
%                            Reassign specified components of one cluster to another cluster. 
% Usage:    
%                   >> [STUDY] = cls_movecomp(STUDY, ALLEEG, from_cluster, to_cluster, comps);   
% Inputs:
%   STUDY         - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG        - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY. 
%                       ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%   from_cluster - cluster number, the cluster components are moved from.  
%   to_cluster     - cluster number, the cluster components are moved to.  
%   comps          - [numeric vector] component indices in the from_cluster to move.  
%
% Outputs:
%   STUDY    - the input STUDY set structure modified with the components reassignment. 
%
%   Example:
%                         >> from_cluster = 10; to_cluster = 7; comps = [2 7];   
%                         >> [STUDY] = cls_movecomp(STUDY,ALLEEG, from_cluster, to_cluster, comps);
%                    Components 2 and 7 of cluster 10 are moved to cluster 7. 
%
%  See also  pop_clustedit         
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

function STUDY = cls_movecomp(STUDY, ALLEEG, old_clus, new_clus, comps)
icadefs;
% Cannot move components if clusters are 'Notclust' clusters
if strncmpi('Notclust',STUDY.cluster(old_clus).name,8) | strncmpi('Notclust',STUDY.cluster(new_clus).name,8)
    warndlg2('cls_movecomp: cannot move components to or from Notclust cluster');
    return;
end

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
    warndlg2('Cannot move components if clusters have different parent clusters!' , 'Aborting move components');
    return;
end
% Check if all parents are the same
if (~strcmp(STUDY.cluster(new_clus).parent, 'manual')) 
    if ~(sum(strcmp(STUDY.cluster(old_clus).parent, STUDY.cluster(new_clus).parent)) == length(STUDY.cluster(new_clus).parent))% different parent
        warndlg2('Cannot move components if clusters have different parent clusters!' , 'Aborting move components');
        return;
    end
end
try 
    h_wait = waitbar(0,['Moving components ...'], 'Color', BACKEEGLABCOLOR,'position', [300, 200, 300, 48]);
catch % for Matlab 5.3
    h_wait = waitbar(0,['Moving components ...'],'position', [300, 200, 300, 48]);
end        
for  ci = 1:length(comps)
    comp =  STUDY.cluster(old_clus).comps(comps(ci));
    sets = STUDY.cluster(old_clus).sets(:,comps(ci));
    %update new cluster
    STUDY.cluster(new_clus).comps(end+1) = comp;%with comp index
    STUDY.cluster(new_clus).sets(:,end+1) = sets; %with set indices 
    STUDY.cluster(new_clus).preclust.preclustdata(end+1,:) = STUDY.cluster(old_clus).preclust.preclustdata(comps(ci),:); %with preclustdata
    if strcmpi(STUDY.cluster(new_clus).parent, 'manual')
        STUDY.cluster(new_clus).preclust.preclustparams = STUDY.cluster(old_clus).preclust.preclustparams;
        STUDY.cluster(new_clus).parent = STUDY.cluster(old_clus).parent;
        STUDY.cluster(find(strcmp({STUDY.cluster.name},STUDY.cluster(new_clus).parent))).child{end+1} = STUDY.cluster(new_clus).name;
    end
    %sort by sets
    [tmp,sind] = sort(STUDY.cluster(new_clus).sets(1,:));
    waitbar(ci/((length(comps)+1)*3),h_wait)
    STUDY.cluster(new_clus).sets = STUDY.cluster(new_clus).sets(:,sind);
    STUDY.cluster(new_clus).comps = STUDY.cluster(new_clus).comps(sind);
    STUDY.cluster(new_clus).preclust.preclustdata = STUDY.cluster(new_clus).preclust.preclustdata(sind,:);
end
STUDY = cls_centroid(STUDY, ALLEEG, new_clus);%new centroid
waitbar((3*(ci-1)+2)/((length(comps)+1)*3),h_wait)
% Remove data from old cluster
% left_comps - are all the components of the cluster after the
% components that were moved to the new cluster were removed. 
left_comps = find(~ismember([1:length(STUDY.cluster(old_clus).comps)],comps));
STUDY.cluster(old_clus).comps = STUDY.cluster(old_clus).comps(left_comps);
STUDY.cluster(old_clus).sets = STUDY.cluster(old_clus).sets(:,left_comps);
STUDY.cluster(old_clus).preclust.preclustdata = STUDY.cluster(old_clus).preclust.preclustdata(left_comps,:);
STUDY = cls_centroid(STUDY, ALLEEG, old_clus);%new centroid
waitbar((3*(ci-1)+3)/((length(comps)+1)*3),h_wait)
waitbar(1,h_wait)
delete(h_wait)
