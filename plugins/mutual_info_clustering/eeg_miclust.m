% eeg_miclust() - Calculate mutual information matrix between independent component
%                 activations and cluster the components using it. Opens several figures:
%                -  A figure with sorted silhouette values of each cluster
%                -  A figure showing clusters in 3-D, with brightness of
%                     each point proportional to its silhouette value
%                -  Multiple 2-D plots, each for one cluster. Background 
%                     color shows interpolated silhouette value. A bright 
%                     background means the component fits well in the cluster;
%                     a dark background means it belongs near equally 
%                     to another cluster or to no cluster. 
%                -  A dendrogram of componets if 'linkage' method is used
% Usage:
%              >> EEG = eeg_miclust(EEG);
%              >> EEG = eeg_miclust(EEG,N,components, clusterMethod);
% Inputs:
%
%   EEG        - EEG data structure
%   N          - number of clusters to produce {default: 5}
%   components - a vector containing component indices (rows of sim matrix) 
%                to cluster. For example [1:10] uses only first 10 components 
%                {Default: [1:120] or all} 
%
% Optional Inputs:
%
%   clusterMethod - which clustering method to use, options are 'linkage'
%                   and 'kmeans' {defgault = 'linkage'}
% Output:
%
%   EEG        - input EEG structure containing mutual information
%                between specified components in field EEG.etc.miclust.mutual_info and indices of 
%                these components in field EEG.etc.miclust.allcomponents Cluster information 
%                is placed in EEG.etc.miclust field.
% Example:
%
%   % Cluster components 1 to 30 into 4 clusters using mutual information.
%   >> EEG = eeg_miclust(EEG,4,1:30);
%
% See also: %   getmiclusts(), showmiclusts(), mi_pairs()
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

% Edit history:
%
% 5/26/06 added default for N, test for <120 comps, extended help msg, formatted -sm
% 5/30/06 use only first 120 comps -nb

function EEGout = eeg_miclust(EEG,n,comps, clusterMethod)


if nargin<2
   n = 5; % default to 5 clusters
end
if nargin<3
    if size(EEG.icaact,1) >= 120
      comps = 120;
    else
      comps = 1:size(EEG.icaact,1);
    end;
end;

if nargin<4
   clusterMethod = 'linkage' ;
else
    if ~strcmp(clusterMethod, 'kmeans') && ~strcmp(clusterMethod, 'linkage')
        fprintf('Clustering method not found.');
        return;
    end;
end

if ~(isfield(EEG,'etc') && isfield(EEG.etc,'miclust') && isfield(EEG.etc.miclust,'mutual_info')  ) || isempty(EEG.etc.miclust.allcomponents) || ~isfield(EEG.etc.miclust,'allcomponents') || length(EEG.etc.miclust.allcomponents) ~= length(comps) || sum(EEG.etc.miclust.allcomponents ~= comps) > 0
    fprintf('Calculating mutual information matrix (time consuming)...');
    mutual_info = mi_pairs(EEG.icaact(comps,:),'full',false);
else
    mutual_info = EEG.etc.miclust.mutual_info;
end

EEG.etc.miclust = getmiclusts(mutual_info,n, comps, clusterMethod);

EEG.etc.miclust.mutual_info = mutual_info;


plotTopo = @(x) topoplot(EEG.icawinv(:,x), EEG.chanlocs,'electrodes','off');

showmiclusts(EEG.etc.miclust,plotTopo);

if nargout>0
    EEGout = EEG;
end;
