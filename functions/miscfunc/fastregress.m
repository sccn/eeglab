% fastregress - perform fast regression and return p-value
%
% Usage:
% [ypred, alpha, rsq, B] = myregress(x, y, plotflag);
%
% Inputs
%  y - y values
%  x - x values
%  plotflag - [0|1] plot regression
%
% Outputs
%  ypred - y prediction
%  alpha - significance level
%  R^2   - r square
%  slope - slope of the fit
%
% Arnaud Delorme, 25 Feb 2003

function [ypred, alpha, rsq, B, BINT] = fastregress(x, y, plotflag)
    
    if nargin < 1
        help fastregress; return;
    end;
    
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
    B = B(2);

    if nargin > 2
        hold on;
        [ynew tmp] = sort(ypred);
        xnew       = x(tmp);
        plot(xnew, ynew, 'r');
        legend(sprintf('R^2=%f', rsq), sprintf('p  =%f', alpha));
    end;