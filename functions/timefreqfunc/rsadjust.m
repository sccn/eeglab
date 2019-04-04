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

function [l1,l2,l3,l4] = rsadjust( l3, l4, mu, sigma2, m3);
    
    % swap l3 and l4 for negative skewness
    % ------------------------------------
    if m3 < 0
        ltmp = l4;
        l4   = l3;
        l3   = ltmp;
    end
    
    A = 1/(1 + l3) - 1/(1 + l4);
    B = 1/(1 + 2*l3) + 1/(1 + 2*l4) - 2*beta(1+l3, 1+l4);
    C = 1/(1 + 3*l3) - 1/(1 + 3*l4) ...
             - 3*beta(1+2*l3, 1+l4) + 3*beta(1+l3, 1+2*l4);

    % compute l2 (and its sign)
    % ------------------------
    l2 = sqrt( (B-A^2)/sigma2 );    
    if m3 == 0, m3 = -0.000000000000001; end
    if (m3*(C - 2*A*B + 2*A^3)) < 0, l2 = -l2; end
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
