% std_clustmaxelec() - function to find the electrode with maximum absolute projection
%                  for each component of a cluster
% Usage:    
%     >> [STUDY, ALLEEG] = std_clustmaxelec(STUDY, ALLEEG, clustind); 
%
% Inputs:
%   STUDY        - STUDY set structure containing (loaded) EEG dataset structures
%   ALLEEG       - ALLEEG vector of EEG structures, else a single EEG dataset.
%   clustind     - (single) cluster index
%
% Outputs:
%   eleclist     - [cell] electrode list
%   setlist      - [integer] set indices for the cluster
%   complist     - [integer] component indices for the cluster
%
% Authors: Claire Braboszcz & Arnaud Delorme , CERCO, UPS/CRNS, 2011

% Copyright (C) Claire Braboszcz, CERCO, UPS/CRNS, 2011
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

function std_clustmaxelec(STUDY, ALLEEG, clusterind);

if nargin < 1
    help findicmaxelec;
    return;
end;

fprintf('Finding electrodes with max weight for cluster %d\n', clusterind);
fprintf('-------------------------------------------------\n');

for index = 1:length(STUDY.cluster(clusterind).comps)
    set  = STUDY.cluster(clusterind).sets(1,index);
    comp = STUDY.cluster(clusterind).comps( index);
    
    [tmp maxelec] = max( abs(ALLEEG(set).icawinv(:, comp)) );
    indelec = ALLEEG(set).icachansind(maxelec);
    maxallelec{index} = ALLEEG(set).chanlocs(indelec).labels;
    allelec = unique_bc(maxallelec);

    fprintf('The electrode with the max weight for component %d of dataset %d is "%s"\n', comp, set, maxallelec{index});
end;
    
 
for indelec=1:length(allelec)  
    nbelec{indelec} = length(find(strcmp(allelec{indelec}, maxallelec) == 1));
    
    fprintf('Number of occurrence of electrode %s: %d\n', allelec{indelec}, nbelec{indelec});  
end;
     

