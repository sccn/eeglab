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
%                  zerofreq is the frequency of occurence of p=0.
%    'zeromode'  - ['on'|'off'] enable estimation of frequency of pval=0
%                  (this might lead to hight pval). Default is 'on'.
%
% Outputs:
%    p          - corrected p value.
%    phat       - phat gamfit() parameter.
%    pci        - phat gamfit() parameter.
%    zerofreq   - frequency of occurence of p=0.
%
% Authors: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2003-
%
% See also: bootstat()

% Copyright (C) 7/02/03  Arnaud Delorme, SCCN/INC/UCSD
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

function [pval, PHAT, PCI, zerofreq] = correctfit(pval, varargin)
    
    if nargin < 2
        help correctfit;
        disp('You need to specify one optional input');
        return;
    end;
    
    g = finputcheck( varargin, { 'allpval'    'real'    [0 1]          [];
                                 'zeromode'   'string'  {'on','off'}   'on';
                                 'gamparams'  'real'    []             []}, 'correctfit');
    if isstr(g), error(g); end;
    
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
    end;
    
    if pval == 0
        if strcmpi(g.zeromode, 'on')
            pval = zerofreq;
        end;
    else
        tmppval = -log10( pval ) + 1E-10;
        pval    = 1-gamcdf( tmppval, PHAT, PCI);
    end;
    
    if 1 % plotting
        if exist('tmpdat') == 1
            figure; hist(tmpdat, 100); hold on; 
            mult = ylim;
            tmpdat = linspace(0.00001,10, 300);
            normy = gampdf( tmpdat, PHAT, PCI);
            plot( tmpdat, normy/max(normy)*mult(2)*3, 'r');
        end;
    end;
