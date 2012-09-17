function pvals = stat_surrogate_pvals(distribution,observed,tail)
% compute empirical p-vals under the null hypothesis that observed samples
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

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
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

numDims = myndims(distribution);

% append observed to last dimension of surrogate distribution
distribution = cat(numDims,distribution,observed);

numDims = myndims(distribution);

% sort along last dimension (replications)
[tmp idx] = sort( distribution, numDims,'ascend');
[tmp mx]  = max( idx,[], numDims);

len = size(distribution,  numDims );
pvals = 1-(mx-0.5)/len;
if strcmpi(tail, 'both')
    pvals = min(pvals, 1-pvals);
    pvals = 2*pvals;
end;


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
    end;
end;