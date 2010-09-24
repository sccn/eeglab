% std_moveoutlier()  - Commandline function, to reassign specified outlier component(s) 
%                            from a cluster to its outlier cluster. 
% Usage:    
%                   >> STUDY = std_moveoutlier(STUDY, ALLEEG, from_cluster, comps);   
% Inputs:
%   STUDY         - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG        - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY. 
%                       ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%   from_cluster - cluster number, the cluster outlier components are moved from.  
%   comps          - [numeric vector] component indices in the from_cluster to move.  
%
% Outputs:
%   STUDY    - the input STUDY set structure modified with the components reassignment. 
%
%   Example:
%                         >> from_cluster = 10; comps = [2 7];   
%                         >> STUDY = std_movecomp(STUDY,ALLEEG, from_cluster, to_cluster, comps);
%                    Components 2 and 7 of cluster 10 are moved to the its outlier cluster ('Outliers Cls 10'). 
%
%  See also  pop_clustedit         
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

function STUDY = std_moveoutlier(STUDY, ALLEEG, old_clus, comps)

% Cannot move components if the cluster is either a 'Notclust' or
% 'Outliers' cluster
if strncmpi('Notclust',STUDY.cluster(old_clus).name,8) | strncmpi('Outliers',STUDY.cluster(old_clus).name,8)
    warndlg2('std_moveoutlier: cannot move components from Notclust or Outliers cluster');
    return;
end

% Cannot move components if clusters have children clusters
if ~isempty(STUDY.cluster(old_clus).child)   
    warndlg2('Cannot move components if cluster has children clusters!' , 'Aborting remove outliers');
    return;
end

outlier_clust = std_findoutlierclust(STUDY,old_clus); %find the outlier cluster for this cluster
if outlier_clust == 0 %no such cluster exist
    STUDY = std_createclust(STUDY, ALLEEG, ['Outliers ' STUDY.cluster(old_clus).name]); %create an outlier cluster
    outlier_clust = length(STUDY.cluster);
end

%move the compnents to the outliers cluster
STUDY = std_movecomp(STUDY, ALLEEG, old_clus, outlier_clust, comps);   
