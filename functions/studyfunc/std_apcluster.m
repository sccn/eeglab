% std_apcluster() - Affinity propagation cluster for eeglab STUDY
%                   Wrapper for function provided by (Frey/Dueck, Science 2007)
%
% Usage:
%   >>  [IDX,C,sumd] = std_apcluster(clustdata);
%
% Inputs:
% clustdata -Preclustering matrix
%
% Optional inputs:
%
%   'maxits'     -maximum number of iterations (default: 1000)
%   'convits'    -if the estimated exemplars stay fixed for convits
%                 iterations, APCLUSTER terminates early (default: 100)
%   'dampfact'   -update equation damping level in [0.5, 1).  Higher
%                 values correspond to heavy damping, which may be needed
%                 if oscillations occur. (default: 0.9)
%   'dist'      - Same as in pdist:
%       'euclidean'   - Euclidean distance (default)
%       'seuclidean'  - Standardized Euclidean distance. Each coordinate
%                       difference between rows in X is scaled by dividing
%                       by the corresponding element of the standard
%                       deviation S=NANSTD(X). To specify another value for
%                       S, use D=PDIST(X,'seuclidean',S).
%       'cityblock'   - City Block distance
%       'minkowski'   - Minkowski distance. The default exponent is 2. To
%                       specify a different exponent, use
%                       D = PDIST(X,'minkowski',P), where the exponent P is
%                       a scalar positive value.
%       'chebychev'   - Chebychev distance (maximum coordinate difference)
%       'mahalanobis' - Mahalanobis distance, using the sample covariance
%                       of X as computed by NANCOV. To compute the distance
%                       with a different covariance, use
%                       D =  PDIST(X,'mahalanobis',C), where the matrix C
%                       is symmetric and positive definite.
%       'cosine'      - One minus the cosine of the included angle
%                       between observations (treated as vectors)
%       'correlation' - One minus the sample linear correlation between
%                       observations (treated as sequences of values).
%       'spearman'    - One minus the sample Spearman's rank correlation
%                       between observations (treated as sequences of values).
%       'hamming'     - Hamming distance, percentage of coordinates
%                       that differ
%       'jaccard'     - One minus the Jaccard coefficient, the
%                       percentage of nonzero coordinates that differ
%
% Outputs:
%
% See also:
%   std_plotinfocluster
%
% Author: Ramon Martinez-Cancino, SCCN, 2014
%
% Copyright (C) 2014  Ramon Martinez-Cancino,INC, SCCN
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

function [idx,c,sumd] = std_apcluster(clustdata,varargin)

%  Input stuffs
try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            opts.(options{i}) = options{i+1};
        end
    else opts= [];
    end
catch
    disp('std_infocluster() error: calling convention {''key'', value, ... } error'); return;
end

try opts.maxits;           catch, opts.maxits       = 200;          end; % Maximum number of iterations
try opts.convits;          catch, opts.convits      = 100;          end; %
try opts.dampfact;         catch, opts.dampfact     = 0.9;          end; %
try opts.dist;             catch, opts.dist         = 'euclidean';  end; % Distance metric

% Getting distance matrix
S = squareform(pdist(clustdata, opts.dist));

% Spell
[idxtmp,netsim,~,expref] = apcluster(S,mean(S(:)),'maxits', opts.maxits, 'convits', opts.convits, 'dampfact',opts.dampfact);

% --- Adjusting output formats --
% Getting centroids
rmpindx = 1:length(idxtmp);
c_indx  = find(rmpindx(:) == idxtmp(:));
c  = clustdata(c_indx,:);

% Reindexing indxtmp
centroids_realindx = unique(idxtmp);
idx = zeros(size(idxtmp));

dist_tmp = [];
for i = 1 : length(centroids_realindx)
    hit_tmp = find(idxtmp == centroids_realindx(i));
    idx(hit_tmp)  = i;
    
    for j = 1:length(hit_tmp)
        dist_tmp(j) = pdist([c(i,:)' clustdata(hit_tmp(j),:)']', opts.dist);
    end
    sumd(i) = sum(dist_tmp);
    dist_tmp = [];
end
