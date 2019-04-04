% rspfunc() - sub-function used by rsget()
%
% Usage: res = rspfunc(pval, l, rval)
%
% Input:
%   pval - p-value to optimize
%   l    - [l1 l2 l3 l4] l-values for Ramberg-Schmeiser distribution
%   rval - expected r-value
%
% Output:
%   res  - residual
%
% Author: Arnaud Delorme, SCCN, 2003
%
% See also: rsget()
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

function rp = rspfunc( pval, l, rval);
    
    % for fiting rp with fminsearch
    % -----------------------------
    rp   = l(1) + (pval.^l(3) - (1-pval).^l(4))/l(2);
    rp   = abs(rval-rp);
