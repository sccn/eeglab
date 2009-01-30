%   getmiclusts() - Perform clustering based on a given similarity (or 
%                 dissimilarity) matrix using multi-dimensional scaling 
%                 and K-means clustering. May be applied to a pariwise 
%                 mutual information matrix returned by mi_pairs();
%   Usage:
%                 >>  clusts = getmiclusts(simdata);
%                 >>  clusts = getmiclusts(simdata,N,components, clusterMethod);
%   Inputs:
%
%   simdata       - similarity (or dissimilarity) matrix 
%                   (for example, a mutual information matrix).
%   N             - number of clusters to use {default: 5}
%   components    - [vector] component indices (= rows of simdata matrix) 
%                   to cluster on. For example [1:10] uses only the 
%                   first 10 rows of simdata {default: cluster all rows}
%
%   clusterMethod - which clustering method to use, options are 'linkage'
%                   and 'kmeans' {defgault = 'linkage'}
%
%   Output:   
%
%   clusts        - structure containing information about clusters plus 
%                   2-D, 3-D, and higher dimensional projections of the 
%                   components with fields:
%                    stres2d,stress3d,stress_manyd -> multi-dimensional stresses 
%                         for the best multidimensional scaling (MDS) fits.
%                         The lower the stress, the better the representation.
%                    cluster -> cell array containing information about
%                      each cluster:
%                       components -> list of components in the cluster
%                       coord2d    -> 2-D projection of cluster components
%                       stres2d    -> multi-dimensional stresses for 2-D
%                                     projection of cluster components
%                    coord2d,coord3d,coord_manyd -> coordinates representing 
%                         MDS component projections in different dimensions  
%                    clusterIDs -> cluster ID for each component
%                    silh -> silhouette value for each component
%                    allcomponents -> list of all input components clustered
%   Example:
%
%   % Create four clusters of existing ICA components using mutual information.
%   >>  clusts = getmiclusts(EEG.mutual_info,4,1:30);   
%
%   See also: eeg_miclust(), showmiclusts(), mi_pairs(), mdscale(), pdist(), kmeans()
%
% Author: Nima Bigdely Shamlo, SCCN/INC/UCSD, 2006
   
% Copyright (C) 2006 Nima Bigdely Shamlo, SCCN/INC/UCSD, nima@sccn.ucsd.edu
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

function miclusters = getmiclusts(mi,number_of_clusters,components, clusteringMethod);

if ~exist('mdscale')
    error('Function mdscale not found - requires the Statistics Matlab toolbox.');
end
if ~exist('pdist')
    error('Function pdist not found - requires the Statistics Matlab toolbox.');
end
if ~exist('kmeans')
    error('Function kmeans not found - requires the Statistics Matlab toolbox.');
end
if ~exist('kmeans')
    error('Function kmeans not found - requires the Statistics Matlab toolbox.');
end
HIGH_DIM = 20; % default: compute high-dimensional projection for clustering
K_REPS = 500;  % K-means 'replicates' value

rand('state',0); % so that every run produces the same results

if nargin<2
    number_of_clusters = 5;
end;

if nargin<3
    components = 1:size(mi,1);
end;

if nargin<4
   clusteringMethod = 'linkage' ;
else
    if strcmp(clusteringMethod, 'kmeans') && strcmp(clusteringMethod, 'linkage')
        fprintf('Clustering method not found.');
        return;
    end;
end

if sum(abs( (1:size(mi,1)) - components)) == 0 % if components have indices different from 1:size(mi,1) then use all mi rows, for cases when components start from a number other than 1, like 2:10
    mi = mi(components,:); % choose only input components 
    mi= mi(:,components);   
end;

if trace(abs(mi)) ~=0
    t = mi/2;  % Use mutual information as similarity matrix, else use it as dissimilarity matrix
    for i=1:size(t,1) t(i,i) = 1; end;
else
    t = mi; % use as distance matrix
end

%
% Perform multi-dimensional scaling
%
opts = statset('MaxIter',1500);
[t2d miclusters.stres2d]= mdscale(t,2,'options',opts);
fprintf('Using multidimensional scaling (stress should be in [0 1]):\n');
fprintf('   2-D stress = %6.2f\n',miclusters.stres2d);

[t3d miclusters.stres3d]= mdscale(t,3,'options',opts);
fprintf('   3-D stress = %6.2f\n',miclusters.stres3d);
fprintf('Calculating coordinates in high dimensions, please wait...\n');
[t_manyd miclusters.stresmanyd]= mdscale(t,HIGH_DIM,'options',opts);
fprintf('   High-Dim stress = %6.2f\n',miclusters.stresmanyd);

%
% Perform clustering
%

if strcmp(clusteringMethod , 'kmeans')
    fprintf('clustering coordinates by kmeans.\n');
    if number_of_clusters>1
        [IDX,C,sumd] = kmeans(t_manyd, number_of_clusters,'replicates',K_REPS);
    else
        IDX = ones(size(t_manyd,1),1);
    end;
end;

if strcmp(clusteringMethod , 'linkage')
    fprintf('clustering coordinates by linkage.\n');
    z = linkage (pdist(t_manyd), 'ward');
    IDX = cluster(z,'MaxClust',number_of_clusters);
    figure; 
    dendrogram(z);
    xlabel('Component');
    ylabel('Distance');
end;

[silh] = silhouette(t_manyd,IDX);

for i=1:number_of_clusters
    miclusters.cluster{i}.components = components(find(IDX==i));
    tclust = t(find(IDX==i),:);
    tclust = tclust(:,find(IDX==i));    
    try
      [miclusters.cluster{i}.coord2d miclusters.cluster{i}.stres2d] ...
                        = mdscale(tclust,2,'options',opts);
        p = pdist(miclusters.cluster{i}.coord2d);
        if min(p)/max(p)<.05
            fprintf(...
['Could not generate 2-D coords. dedicated to cluster ' num2str(i) ', using general coords.!\n']);
            miclusters.cluster{i}.coord2d = t2d(IDX==i,:);
            miclusters.cluster{i}.stres2d = miclusters.stres2d;
        end;
    catch
        fprintf(...
['Could not generate 2-D coords. dedicated to cluster ' num2str(i) ', using general coords.!\n']);
        miclusters.cluster{i}.coord2d = t2d(IDX==i,:);
        miclusters.cluster{i}.stres2d = miclusters.stres2d;
    end;
end;

%
% Save cluster structure
%
miclusters.coord2d = t2d;
miclusters.coord3d = t3d;
miclusters.coord_manyd = t_manyd;
miclusters.clusterIDs = IDX;
miclusters.number_of_clusters = number_of_clusters;
miclusters.silh = silh;
miclusters.allcomponents = components;
