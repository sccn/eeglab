% std_readtopoclust() - Compute and return cluster component scalp maps. 
%                  Automatically inverts the polarity of component scalp maps 
%                  to best match the polarity of the cluster mean scalp map.
% Usage:    
%              >> [STUDY clsstruct] = std_readtopoclust(STUDY, ALLEEG, clusters);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the dataset(s) included in 
%                the STUDY. 
%   clusters   - cluster numbers to read.
%
% Outputs:
%   STUDY      - the input STUDY set structure with the computed mean cluster scalp
%                map added (unless cluster scalp map means already exist in the STUDY) 
%                to allow quick replotting. 
%   clsstruct  - STUDY.cluster structure array for the modified clusters.
%
% See also  std_topoplot(), pop_clustedit()
%
% Authors:  Arnaud Delorme, SCCN, INC, UCSD, 2007

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, June 07, 2007, arno@sccn.ucsd.edu
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

function [STUDY, centroid] = std_readtopoclust(STUDY,ALLEEG, clsind);

if nargin < 3
    help readtopoclust;
    return;
end
    
if isempty(clsind)
    for k = 2: length(STUDY.cluster) % don't include the ParentCluster
         if ~strncmpi('Notclust',STUDY.cluster(k).name,8) 
             % don't include 'Notclust' clusters
             clsind = [clsind k];
         end
    end
end

 Ncond = length(STUDY.condition);
if Ncond == 0
    Ncond = 1;
end 
centroid = cell(length(clsind),1);
fprintf('Computing the requested mean cluster scalp maps (only done once)\n');
if ~isfield( STUDY.cluster, 'topo' ), STUDY.cluster(1).topo = []; end
cond = 1;

for clust = 1:length(clsind) %go over all requested clusters
    
    if isempty( STUDY.cluster(clsind(clust)).topo )

        numitems = length(STUDY.cluster(clsind(clust)).comps);
        
        for k = 1:numitems % go through all components
            comp  = STUDY.cluster(clsind(clust)).comps(k);
            abset = STUDY.cluster(clsind(clust)).sets(cond,k);
            if ~isnan(comp) && ~isnan(abset)
                [grid yi xi] = std_readtopo(ALLEEG, abset, comp);
                if ~isfield(centroid{clust}, 'topotmp') || isempty(centroid{clust}.topotmp)
                    centroid{clust}.topotmp = zeros([ size(grid(1:4:end),2) numitems ]);
                end
                centroid{clust}.topotmp(:,k) = grid(1:4:end); % for inversion
                centroid{clust}.topo{k} = grid;
                centroid{clust}.topox = xi;
                centroid{clust}.topoy = yi;
            end
        end
        fprintf('\n');

        %update STUDY
        tmpinds = find(isnan(centroid{clust}.topotmp(:,1)));
        %centroid{clust}.topotmp(tmpinds,:) = [];
        %for clust =  1:length(clsind) %go over all requested clusters
            for cond  = 1
                if clsind(1) > 0
                    ncomp = length(STUDY.cluster(clsind(clust)).comps);
                end
                [ tmp pol ] = std_comppol(centroid{clust}.topotmp);
                fprintf('%d/%d polarities inverted while reading component scalp maps\n', ...
                        length(find(pol == -1)), length(pol));
                nitems = length(centroid{clust}.topo);
                for k = 1:nitems
                    centroid{clust}.topo{k} = pol(k)*centroid{clust}.topo{k};
                    if k == 1, allscalp = centroid{clust}.topo{k}/nitems;
                    else       allscalp = centroid{clust}.topo{k}/nitems + allscalp;
                    end
                end
                STUDY.cluster(clsind(clust)).topox   = centroid{clust}.topox;
                STUDY.cluster(clsind(clust)).topoy   = centroid{clust}.topoy;
                STUDY.cluster(clsind(clust)).topoall = centroid{clust}.topo;
                STUDY.cluster(clsind(clust)).topo    = allscalp;
                STUDY.cluster(clsind(clust)).topopol = pol;
            end
        %end
    else
        
        centroid{clust}.topox = STUDY.cluster(clsind(clust)).topox;
        centroid{clust}.topoy = STUDY.cluster(clsind(clust)).topoy;
        centroid{clust}.topo  = STUDY.cluster(clsind(clust)).topoall;
        
    end
    
end

fprintf('\n');
