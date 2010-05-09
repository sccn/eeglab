% rsadjust() - adjust l-values (Ramberg-Schmeiser distribution) 
%                with respect to signal mean and variance
%
% Usage: p = rsadjust(l3, l4, m, var, skew)
%
% Input:
%   l3   - value lambda3 for Ramberg-Schmeiser distribution
%   l4   - value lambda4 for Ramberg-Schmeiser distribution
%   m    - mean of the signal distribution
%   var  - variance of the signal distribution
%   skew - skewness of the signal distribution (only the sign of
%          this parameter is used).
%
% Output:
%   l1  - value lambda3 for Ramberg-Schmeiser distribution
%   l2  - value lambda4 for Ramberg-Schmeiser distribution
%   l3  - value lambda3 for Ramberg-Schmeiser distribution (copy
%         from input)
%   l4  - value lambda4 for Ramberg-Schmeiser distribution (copy
%         from input)
%
% Author: Arnaud Delorme, SCCN, 2003
%
% See also: rsfit(), rsget()
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

function [l1,l2,l3,l4] = rsadjust( l3, l4, mu, sigma2, m3);
    
    % swap l3 and l4 for negative skewness
    % ------------------------------------
    if m3 < 0
        ltmp = l4;
        l4   = l3;
        l3   = ltmp;
    end;
    
    A = 1/(1 + l3) - 1/(1 + l4);
    B = 1/(1 + 2*l3) + 1/(1 + 2*l4) - 2*beta(1+l3, 1+l4);
    C = 1/(1 + 3*l3) - 1/(1 + 3*l4) ...
             - 3*beta(1+2*l3, 1+l4) + 3*beta(1+l3, 1+2*l4);

    % compute l2 (and its sign)
    % ------------------------
    l2 = sqrt( (B-A^2)/sigma2 );    
    if m3 == 0, m3 = -0.000000000000001; end;
    if (m3*(C - 2*A*B + 2*A^3)) < 0, l2 = -l2; end;
    %l22 = ((C - 2*A*B + 2*A^3)/m3)^(1/3) % also equal to l2
       
    % compute l1
    % ----------
    l1 = mu - A/l2;
    
    return;
    
    % fitting table 1 of 
    % ---------------
    [l1 l2]  = pdffitsolv2(-.0187,-.0388, 0, 1, 1)
    [l1 l2]  = pdffitsolv2(-.1359,-.1359, 0, 1, 0)

    % sign problem
    [l1 l2]  = pdffitsolv2(1.4501,1.4501, 0, 1)
    
    % numerical problem for l1
    [l1 l2]  = pdffitsolv2(-.00000407,-.001076, 0, 1, 2)
    [l1 l2]  = pdffitsolv2(0,-.000580, 0, 1, 2)
