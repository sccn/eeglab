% Compute empirical p-vals under the null hypothesis that observed samples
% come from a given surrogate distribution. P-values for Type I error in
% rejecting the null hypothesis are obtained by finding the proportion of
% samples in the distribution that
% (a) are larger than the observed sample (one-sided test)
% (b) are larger or smaller than the observed sample (two-sided test).
%
% This function is based on Arnaud Delorme's statcond:compute_pvals()
% function from EEGLAB
% 
% Inputs:
%
%   distribution:   [d1 x d2 x ... x dM x N] matrix of surrogate samples. 
%                   distribution(i,j,k,...,:) is a collection of N samples
%                   from a surrogate distribution.
%   observed:       [d1 x d2 x ... x dM] matrix of observations.
%   tail:           can be 'one' or 'both' indicating a one-tailed or
%                   two-tailed test
% Outputs:
%   
%   pvals:          [d1 x d2 x ... x dM] matrix of p-values specifying
%                   probability of Type I error in rejecting the null 
%                   hypothesis
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

function pvals = stat_surrogate_pvals(distribution,observed,tail)
numDims = ndims(distribution);
if iscolumn(distribution)
	numDims = 1;
end

n = size(distribution, numDims);
% once support for matlab <= R2016b is dropped:
% pvals = sum(distribution >= observed, numDims) / n;
pvals = sum(bsxfun(@ge, distribution, observed), numDims) / n;

if any(strcmpi(tail, {'right', 'one'}))
	% nothing to be done
	return;
end

p_left = 1 - pvals + sum(bsxfun(@eq, distribution, observed), numDims) / n;

if strcmpi(tail, 'both')
	pvals = 2 * min(pvals, p_left);
elseif strcmpi(tail, 'left')
	pvals = p_left;
else
	error('invalid value for tail: "%s", should be left, right, one or both', tail);
end

end