% rsget() - get the p-value for a given collection of l-values
%            (Ramberg-Schmeiser distribution)
%
% Usage: p = getfit(l, val)
%
% Input:
%   l   - [l1 l2 l3 l4] l-values for Ramberg-Schmeiser distribution
%   val - value in the distribution to get a p-value estimate at
%
% Output:
%   p   - p-value
%
% Author: Arnaud Delorme, SCCN, 2003
%
% see also: rspfunc()
%
% Reference: Ramberg, J.S., Tadikamalla, P.R., Dudewicz E.J., Mykkytka, E.F.
%            A probability distribution and its uses in fitting data. 
%            Technimetrics, 1979, 21: 201-214.

% Copyright (C) 2003 Arnaud Delorme, SCCN, arno@salk.edu
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

function p = rsget( l, val);
    
    % plot the curve
    % --------------
    %pval = linspace(0,1, 102); pval(1) = []; pval(end) = [];
    %rp   = l(1) + (pval.^l(3) - (1-pval).^l(4))/l(2);
    %fp   = l(2)*1./(l(3).*(pval.^(l(3)-1)) + l(4).*((1-pval).^(l(4)-1)));
    %figure; plot(pval, rp);
    %figure; plot(rp, fp);
    
    % find p value for a given val
    % ----------------------------
    p = fminbnd('rspfunc', 0, 1, optimset('TolX',1e-300), l, val);
