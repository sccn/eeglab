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

function STUDY = std_moveoutlier(STUDY, ALLEEG, old_clus, comps)

% Cannot move components if the cluster is either a 'Notclust' or
% 'Outliers' cluster
if strncmpi('Notclust',STUDY.cluster(old_clus).name,8) || strncmpi('Outliers',STUDY.cluster(old_clus).name,8)
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
