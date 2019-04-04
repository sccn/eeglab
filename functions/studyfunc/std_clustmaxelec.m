% std_clustmaxelec() - function to find the electrode with maximum absolute projection
%                  for each component of a cluster
% Usage:    
%     >> [eleclist, setlist, complist] = std_clustmaxelec(STUDY, ALLEEG, clustind); 
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

function std_clustmaxelec(STUDY, ALLEEG, clusterind);

if nargin < 1
    help findicmaxelec;
    return;
end

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
end
    
 
for indelec=1:length(allelec)  
    nbelec{indelec} = length(find(strcmp(allelec{indelec}, maxallelec) == 1));
    
    fprintf('Number of occurrence of electrode %s: %d\n', allelec{indelec}, nbelec{indelec});  
end
     

