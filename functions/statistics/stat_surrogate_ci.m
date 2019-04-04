% compute empirical p-vals under the null hypothesis that observed samples
% come from a given surrogate distribution. P-values for Type I error in
% rejecting the null hypothesis are obtained by finding the proportion of
% samples in the distribution that
% (a) are larger than the observed sample (one-sided test)
% (b) are larger or smaller than the observed sample (two-sided test).
%
% Inputs:
%   distribution:   [d1 x d2 x ... x dM x N] matrix of surrogate samples. 
%                   distribution(i,j,k,...,:) is a collection of N samples
%                   from a surrogate distribution.
%   alpha:          [float] alpha value. Default
%   tail:           can be 'one' or 'both' indicating a one-tailed or
%                   two-tailed test. Can also be 'lower' or 'upper'.
% Outputs:
%   
%   ci:             [2 x d1 x d2 x ... x dM] matrix of confidence interval
% 
% Author: Tim Mullen and Arnaud Delorme, SCCN/INC/UCSD

% Copyright: Tim Mullen and Arnaud Delorme
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

function ci = stat_surrogate_ci(distribution,alpha,tail)

if nargin<3
    tail = 'both';
end

if nargin<2
    alpha = 0.05;
end

% reshape matrix
% --------------
nd = size(distribution);
if length(nd) == 2 && nd(2) == 1, nd(2) = []; end
ndori = nd;
nd = nd(1:end-1);
if isempty(nd), nd = 1; end
distribution = reshape(distribution, [prod(nd) size(distribution,myndims(distribution))]);

% append observed to last dimension of surrogate distribution
numDims = myndims(distribution);

% number of samples
N = size(distribution, numDims);

% sort along last dimension (replications)
[tmpsort idx] = sort( distribution, numDims,'ascend');

if strcmpi(tail, 'both'), alpha = alpha/2; end

low  = round(alpha*N);
high = N-low;
switch lower(tail)
    case 'upper'
        cilow  = mean(tmpsort, numDims);
        cihigh = tmpsort(:,high);
    case 'lower'
        cilow   = tmpsort(:,low+1);
        cihigh  = mean(tmpsort, numDims);
    case { 'both' 'one' }
        cilow  = tmpsort(:,low+1);
        cihigh = tmpsort(:,high);
    otherwise
        error('Unknown tail option');
end

ci = reshape(cilow, [1 size(cilow)]);
ci(2,:) = cihigh;
ci = reshape(ci, [2 ndori(1:end-1) 1]);

% get the number of dimensions in a matrix
function val = myndims(a)
if ndims(a) > 2
    val = ndims(a);
else
    if size(a,1) == 1,
        val = 2;
    elseif size(a,2) == 1,
        val = 1;
    else
        val = 2;
    end
end
