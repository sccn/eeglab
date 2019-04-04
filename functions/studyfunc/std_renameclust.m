% std_renameclust()  - Commandline function, to rename clusters using specified (mnemonic) names. 
% Usage:    
%                   >> [STUDY] = std_renameclust(STUDY, ALLEEG, cluster, new_name);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY. 
%                     ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%   cluster     - single cluster number.  
%   new_name  - [string] mnemonic cluster name.
%
% Outputs:
%   STUDY    - the input STUDY set structure modified according to specified new cluster name. 
%
%   Example:
%                         >> cluster = 7; new_name = 'artifacts';  
%                         >> [STUDY] = std_renameclust(STUDY,ALLEEG, cluster, new_name);
%                    Cluster 7 name (i.e.: STUDY.cluster(7).name) will change to 'artifacts 7'. 
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

function STUDY = std_renameclust(STUDY, ALLEEG, cls, new_name)

if ~exist('cls')
    error('std_renameclust: you must provide a cluster number to rename.');
end
if isempty(cls)
   error('std_renameclust: you must provide a cluster number to rename.');
end
if ~exist('new_name')
    error('std_renameclust: you must provide a new cluster name.');
end
if strncmpi('Notclust',STUDY.cluster(cls).name,8)  % Don't rename Notclust 'clusters'
    warndlg2('std_renameclust: Notclust cannot be renamed');
    return;
end

ti = strfind(STUDY.cluster(cls).name, ' ');
clus_id = STUDY.cluster(cls).name(ti(end) + 1:end);
new_name = sprintf('%s %s', new_name, clus_id);
% If the cluster have children cluster update their parent cluster name to the
% new cluster.
if ~isempty(STUDY.cluster(cls).child)
    for k = 1:length(STUDY.cluster(cls).child)
        child_cls = STUDY.cluster(cls).child{k};
        child_id = find(strcmp({STUDY.cluster.name},child_cls));
        parent_id = find(strcmp(STUDY.cluster(child_id).parent,STUDY.cluster(cls).name));
        STUDY.cluster(child_id).parent{parent_id} = new_name;
    end
end
% If the cluster has parent clusters, update the parent clusters with the 
% new cluster name of child cluster.
if ~isempty(STUDY.cluster(cls).parent)
    for k = 1:length(STUDY.cluster(cls).parent)
        parent_cls = STUDY.cluster(cls).parent{k};
        parent_id = find(strcmp({STUDY.cluster.name},parent_cls));
        STUDY.cluster(parent_id).child{find(strcmp(STUDY.cluster(parent_id).child,STUDY.cluster(cls).name))} = new_name;
    end
end
% If the cluster have an Outlier cluster, update the Outlier cluster name.  
outlier_clust = std_findoutlierclust(STUDY,cls); %find the outlier cluster for this cluster
if outlier_clust ~= 0
    ti = strfind(STUDY.cluster(outlier_clust).name, ' ');
    clus_id = STUDY.cluster(outlier_clust).name(ti(end) + 1:end);    
    % If the outlier has parent clusters, update the parent clusters with the 
    % new cluster name of child cluster.
    for k = 1:length(STUDY.cluster(outlier_clust).parent)
        parent_cls = STUDY.cluster(outlier_clust).parent{k};
        parent_id = find(strcmp({STUDY.cluster.name},parent_cls));
        STUDY.cluster(parent_id).child{find(strcmp(STUDY.cluster(parent_id).child,STUDY.cluster(outlier_clust).name))} = sprintf('Outliers %s %s', new_name, clus_id);
    end
    
    STUDY.cluster(outlier_clust).name = sprintf('Outliers %s %s', new_name, clus_id);
end

% Rename cluster
STUDY.cluster(cls).name = new_name;
