% fastregress - perform fast regression and return p-value
%
% Usage:
% [ypred, alpha, rsq, slope, intercept] = fastregress(x, y);
% [ypred, alpha, rsq, slope, intercept] = fastregress(x, y, plotflag, plotleg);
%
% Inputs
%  y - y values
%  x - x values
%  plotflag - [0|1] plot regression and legend. Default 0.
%  plotleg  - [0|1] plot legend. Default same as plotflag.
%
% Outputs
%  ypred - y prediction
%  alpha - significance level
%  R^2   - r square
%  slope - slope of the fit
%  intercept - intercept of the fit
%
% Arnaud Delorme, 25 Feb 2003

% Copyright (C) Arnaud Delorme, 25 Feb 2003
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

function [ypred, alpha, rsq, B, intercept, h] = fastregress(x, y, plotflag, plotleg)
    
    if nargin < 1
        help fastregress; return;
    end
    if nargin < 3
        plotflag = false;
    end
    if nargin < 4
        plotleg = plotflag;
    end
    
    % this part is useless but still works
    %B=polyfit(x, y, 1);         % B is the slope
    %ypred = polyval(B,x);       % predictions
    %dev = y - mean(y);          % deviations - measure of spread
    %SST = sum(dev.^2);          % total variation to be accounted for
    %resid = y - ypred;          % residuals - measure of mismatch
    %SSE = sum(resid.^2);        % variation NOT accounted for
    %rsq = 1 - SSE/SST;          % percent of error explained
    % see the page  http://www.facstaff.bucknell.edu/maneval/help211/fitting.html

    [B,BINT,R,RINT,STATS] = regress(y(:), [ones(length(x),1) x(:)]);
    alpha = STATS(3);
    rsq   = STATS(1);
    
    %note that we also have 
    %ypred =  [ones(size(x,2),1) x(:)]*B;
    ypred =  x*B(2) + B(1);

    % B(1) contain the offset, B(2) the slope
    intercept = B(1);
    B = B(2);

    if plotflag
        hold on;
        [ynew, tmp] = sort(ypred);
        xnew        = x(tmp);
        h = plot(xnew, ynew, 'r');
    end
    if plotleg
        legend(sprintf('R^2=%f', rsq), sprintf('p  =%f', alpha));
    end
