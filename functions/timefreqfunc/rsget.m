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
