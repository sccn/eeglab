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
if length(nd) == 2 && nd(2) == 1, nd(2) = []; end;
ndori = nd;
nd = nd(1:end-1);
if isempty(nd), nd = 1; end;
distribution = reshape(distribution, [prod(nd) size(distribution,myndims(distribution))]);

% append observed to last dimension of surrogate distribution
numDims = myndims(distribution);

% number of samples
N = size(distribution, numDims);

% sort along last dimension (replications)
[tmpsort idx] = sort( distribution, numDims,'ascend');

if strcmpi(tail, 'both'), alpha = alpha/2; end;

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
    end;
end;
