% rspdfsolv() - sub-function used by rsfit() to searc for optimal
%               parameter for Ramberg-Schmeiser distribution
%
% Usage: res = rspdfsolv(l, l3, l4)
%
% Input:
%   l    - [lambda3 lamda4] parameters to optimize
%   skew - expected skewness
%   kurt - expected kurtosis
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

function res = rspdfsolv( l, a3, a4);
    
    A = 1/(1 + l(1)) - 1/(1 + l(2));
    B = 1/(1 + 2*l(1)) + 1/(1 + 2*l(2)) - 2*beta(1+l(1), 1+l(2));
    C = 1/(1 + 3*l(1)) - 1/(1 + 3*l(2)) ...
             - 3*beta(1+2*l(1), 1+l(2)) + 3*beta(1+l(1), 1+2*l(2));
    D = 1/(1 + 4*l(1)) + 1/(1 + 4*l(2)) ...
             - 4*beta(1+3*l(1), 1+l(2)) - 4*beta(1+l(1), 1+3*l(2)) ...
             + 6*beta(1+2*l(1), 1+2*l(2));

    estim_a3 = (C - 3*A*B + 2*A^3)/(B-A^2)^(3/2);    
    estim_a4 = (D - 4*A*C + 6*A^2*B - 3*A^4)/(B-A^2)^2;
    
    res = (estim_a3 - a3)^2 + (estim_a4 - a4)^2;
    
    % the last term try to ensures that l(1) and l(2) are of the same sign
    if sign(l(1)*l(2)) == -1, res = 2*res; end;
    return;
    
%         original equations
% $$$     A = 1(1 + l(3)) - 1/(1 + l(4));
% $$$     B = 1(1 + 2*l(3)) + 1/(1 + 2*l(4)) - 2*beta(1+l(3), 1+l(4));
% $$$     C = 1(1 + 3*l(3)) - 1/(1 + 3*l(4)) ...
% $$$             - 3*beta(1+2*l(3), 1+l(4)) + 3*beta(1+l(3), 1+2*l(4));
% $$$     D = 1(1 + 4*l(3)) + 1/(1 + 4*l(4)) ...
% $$$             - 4*beta(1+3*l(3), 1+l(4)) - 4*beta(1+l(3), 1+3*l(4)) ...
% $$$             + 6*beta(1+2*l(3), 1+2*l(4));
% $$$ 
% $$$     R(1) = l(1) + A/l(2);
% $$$     R(2) = (B-A^2)/l(2)^2;
% $$$     R(3) = (C - 3*A*B + 2*A^3)/l(2)^3;    
% $$$     R(4) = (D - 4*A*C + 6*A^2*B - 3*A^4)/l(2)^4;
    
