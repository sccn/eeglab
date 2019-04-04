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
    if sign(l(1)*l(2)) == -1, res = 2*res; end
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
    
