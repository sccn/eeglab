% correctfit() - correct fit using observed p-values. Use this function
%                if for some reason, the distribution of p values is
%                not uniform between 0 and 1
%
% Usage:
%     >> [p phat pci zerofreq] = correctfit(pval, 'key', 'val');
%
% Inputs:
%    pval  - input p value
%
% Optional inputs:
%    'allpval'   - [float array] collection of p values drawn from random
%                  distributions (theoritically uniform).
%    'gamparams' - [phat pci zerofreq] parameter for gamma function fitting.
%                  zerofreq is the frequency of occurrence of p=0.
%    'zeromode'  - ['on'|'off'] enable estimation of frequency of pval=0
%                  (this might lead to high pval). Default is 'on'.
%
% Outputs:
%    p          - corrected p value.
%    phat       - phat gamfit() parameter.
%    pci        - phat gamfit() parameter.
%    zerofreq   - frequency of occurrence of p=0.
%
% Authors: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2003-
%
% See also: bootstat()

% Copyright (C) 7/02/03  Arnaud Delorme, SCCN/INC/UCSD
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

function [pval, PHAT, PCI, zerofreq] = correctfit(pval, varargin)
    
    if nargin < 2
        help correctfit;
        disp('You need to specify one optional input');
        return;
    end
    
    g = finputcheck( varargin, { 'allpval'    'real'    [0 1]          [];
                                 'zeromode'   'string'  {'on','off'}   'on';
                                 'gamparams'  'real'    []             []}, 'correctfit');
    if ischar(g), error(g); end
    
    if ~isempty(g.gamparams)
        PHAT     = g.gamparams(1);
        PCI      = g.gamparams(2);
        zerofreq = g.gamparams(3);
    elseif ~isempty(g.allpval)
        nonzero     = find(g.allpval(:) ~= 0);
        zerofreq    = (length(g.allpval(:))-length(nonzero))/ length(g.allpval(:));
        tmpdat      = -log10( g.allpval(nonzero) ) + 1E-10;   
        [PHAT, PCI] = gamfit( tmpdat );
        PHAT = PHAT(1);
        PCI  = PCI(2);
    end
    
    if pval == 0
        if strcmpi(g.zeromode, 'on')
            pval = zerofreq;
        end
    else
        tmppval = -log10( pval ) + 1E-10;
        pval    = 1-gamcdf( tmppval, PHAT, PCI);
    end
    
    if 1 % plotting
        if exist('tmpdat') == 1
            figure; hist(tmpdat, 100); hold on; 
            mult = ylim;
            tmpdat = linspace(0.00001,10, 300);
            normy = gampdf( tmpdat, PHAT, PCI);
            plot( tmpdat, normy/max(normy)*mult(2)*3, 'r');
        end
    end
