% fdr() - compute false detection rate mask
%
% Usage:
%   >> [p_fdr, p_masked] = fdr( pvals, alpha);
%
% Inputs:
%   pvals - vector or array of p-values
%   alpha - threshold value (non-corrected). If no alpha is given
%           each p-value is used as its own alpha and FDR corrected
%           array is returned.
%
% Outputs:
%   p_fdr    - pvalue used for threshold (based on independence
%              or positive dependence of measurements)
%   p_masked - p-value thresholded. Same size as pvals.
%
% Author: Arnaud Delorme, SCCN, 2008-
%         Based on a function by Tom Nichols
%
% See also: eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% $Log: not supported by cvs2svn $
% Revision 1.1  2008/05/06 22:32:52  arno
% Initial revision
%

function [pID, p_masked] = fdr(pvals, q);

p = sort(pvals(:));
V = length(p);
I = (1:V)';

cVID = 1;
cVN = sum(1./(1:V));

if nargin < 2
    pID = ones(size(pvals));
    thresholds = exp(linspace(log(0.1),log(0.000001), 100));
    for index = 1:length(thresholds)
        [tmp p_masked] = fdr(pvals, thresholds(index));
        pID(p_masked) = thresholds(index);    
    end;
else
    pID = p(max(find(p<=I/V*q/cVID))); % standard FDR
    %pN = p(max(find(p<=I/V*q/cVN)));  % non-parametric FDR (not used)
end;
if isempty(pID), pID = 0; end;

if nargout > 1
    p_masked = pvals < pID;
end;
