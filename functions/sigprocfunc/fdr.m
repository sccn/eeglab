% fdr() - compute false detection rate mask
%
% Usage:
%   >> [p_masked, p_fdr] = fdr( pvals, alpha);
%
% Inputs:
%   pvals - vector or array of p-values
%   alpha - threshold value (non-corrected)
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

function [pID, p_masked] = fdr(pvals, q);

p = sort(pvals(:));
V = length(p);
I = (1:V)';

cVID = 1;
cVN = sum(1./(1:V));

pID = p(max(find(p<=I/V*q/cVID))); % standard FDR
pN =  p(max(find(p<=I/V*q/cVN)));  % non-parametric FDR (not used)
if isempty(pID), pID = 0; end;

if nargout > 1
    p_masked = pvals < pID;
end;
