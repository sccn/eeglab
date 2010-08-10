% pop_miclust() - Calculate mutual information matrix between independent component
%                 activations and cluster the components using it. If number of input
%                 arguments is less than 3, pop up an interactive query window.
%                 Calls eeg_miclust(). Opens several figures:
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
%              >> EEG = pop_miclust(EEG);
%              >> EEG = pop_miclust(EEG,N,components, clusterMethod);
% Inputs:
%
%   EEG        - EEG data structure
%   N          - number of clusters to produce
%   components - a vector containing component indices (rows of sim matrix)
%                to cluster. For example [1:10] uses only first 10 components
%
%
% Optional Inputs:
%
%   clusterMethod - which clustering method to use, options are 'linkage'
%                   and 'kmeans' {defgault = 'linkage'}
% Output:
%   EEG        - input EEG structure containing mutual information
%                between specified components in field EEG.etc.miclust.mutual_info and indices of
%                these components in field EEG.etc.miclust.allcomponents Cluster information
%                is placed in EEG.etc.miclust field.
%
% Example:
%
%   % Cluster components 1 to 30 into 4 clusters using mutual information.
%   >> EEG = pop_miclust(EEG,4,1:30);
%
% See also: %   eeg_miclust(), getmiclusts(), showmiclusts(), mi_pairs()
%
% Author: Nima Bigdely Shamlo, SCCN/INC/UCSD, 2007


% Copyright (C) 2007 Nima Bigdely Shamlo, SCCN/INC/UCSD, nima@sccn.ucsd.edu
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

function [EEG, com] = pop_miclust( EEG, n , comps, clusterMethod);

% the command output is a hidden output that does not have to
% be described in the header

com = ''; % this initialization ensure that the function will return something
% if the user press the cancel button

% display help if not enough arguments
% ------------------------------------
if nargin < 1
    help pop_miclust;
    return;
end;

if nargin<4
    clusterMethod = 'linkage' ;
else
    if ~strcmp(clusterMethod, 'kmeans') && ~strcmp(clusterMethod, 'linkage')
        fprintf('Clustering method not found.');
        help pop_miclust;
        return;
    end;
end


% pop up window
% -------------
if nargin < 3
    promptstr    = { 'Components to cluster:' , 'Number of clusters:', 'Clustering method:'};
    inistr       = { '1:30','5', 'linkage' };
    result       = inputdlg2( promptstr, 'Cluster Dataset ICs -- pop_miclust', 1,  inistr, 'pop_miclust');
    if length( result ) == 0 return; end;
    
    comps   	 = eval( [ '[' result{1} ']' ] ); % the brackets allow to process matlab arrays
    n       	 = eval( [ '[' result{2} ']' ] ); % the brackets allow to process matlab arrays
    clusterMethod = result{3};
end;


% check whether ica data is availible
% ---------------------------------------------------
EEG.icaact = eeg_getica(EEG);
if ~isempty( EEG.icaact )
    EEG = eeg_miclust(EEG,n,comps, clusterMethod);
else
    error('You must run ICA first');
end;

% return the string command
% -------------------------
com = sprintf('pop_miclust( %s, %s, [%s], ''%s'' );', inputname(1), result{2}, result{1}, result{3});

return;
