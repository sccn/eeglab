% rsfit() - find p value for a given value in a given distribution
%             using Ramberg-Schmeiser distribution
%
% Usage: >> p = rsfit(x, val)
%        >> [p c l chi2] = rsfit(x, val, plot)
%
% Input:
%   x    - [float array] accumulation values
%   val  - [float] value to test
%   plot - [0|1|2] plot fit. Using 2, the function avoids creating
%          a new figure. Default: 0.
%
% Output:
%   p    - p value
%   c    - [mean var skewness kurtosis] distribution cumulants
%   l    - [4x float vector] Ramberg-Schmeiser distribution best fit
%          parameters.
%   chi2 - [float] chi2 for goodness of fit (based on 12 bins). 
%          Fit is significantly different from data histogram if 
%          chi2 > 19 (5%) 
%
% Author: Arnaud Delorme, SCCN, 2003
%
% See also: rsadjust(), rsget(), rspdfsolv(), rspfunc()
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

function [p, c, l, res] = rsfit(x, val, plotflag)

    if nargin < 2
        help rsfit;
        return;
    end;
    if nargin < 3
        plotflag  = 0;
    end;
    
    % moments
    % -------
    m1  = mean(x);
    m2  = sum((x-m1).^2)/length(x);
    m3  = sum((x-m1).^3)/length(x);
    m4  = sum((x-m1).^4)/length(x);

    xmean = m1;
    xvar  = m2;
    xskew = m3/(m2^1.5);
    xkurt = m4/(m2^2);
    c     = [ xmean xvar xskew xkurt ];
    
    if xkurt < 0
        disp('rsfit error: Can not fit negative kurtosis');
        save('/home/arno/temp/dattmp.mat', '-mat', 'x');
        disp('data saved to disk in /home/arno/temp/dattmp.mat');        
    end;
    
    % find fit
    % --------
    try, 
        [sol tmp exitcode] = fminsearch('rspdfsolv', [0.1 0.1], optimset('TolX',1e-12, 'MaxFunEvals', 100000000), abs(xskew), xkurt);    
    catch, exitcode = 0; % did not converge
    end;
    if ~exitcode
        try, [sol tmp exitcode] = fminsearch('rspdfsolv', -[0.1 0.1], optimset('TolX',1e-12, 'MaxFunEvals', 100000000), abs(xskew), xkurt);
        catch, exitcode = 0; end;
    end;
    if ~exitcode,           error('No convergence'); end;
    if sol(2)*sol(1) == -1, error('Wrong sign for convergence'); end;
    %fprintf('          l-val:%f\n', sol);
    
    res = rspdfsolv(sol, abs(xskew), xkurt);
    l3 = sol(1);
    l4 = sol(2);

    %load res;
    %[tmp indalpha3] = min( abs(rangealpha3 - xskew) );
    %[tmp indalpha4] = min( abs(rangealpha4 - xkurt) );
    %l3  = res(indalpha3,indalpha4,1);
    %l4  = res(indalpha3,indalpha4,2);    
    %res = res(indalpha3,indalpha4,3);
    
    % adjust fit
    % ----------
    [l1 l2 l3 l4] = rsadjust(l3, l4, xmean, xvar, xskew);
    l = [l1 l2 l3 l4];
    p = rsget(l, val);

    % compute goodness of fit
    % -----------------------
    if nargout > 3 | plotflag

        % histogram of value 12 bins
        % --------------------------
        [N X] = hist(x, 25);
        interval = X(2)-X(1);
        X = [X-interval/2 X(end)+interval/2]; % borders
        
        % regroup bin with less than 5 values
        % -----------------------------------
        indices2rm = [];
        for index = 1:length(N)-1
            if N(index) < 5
                N(index+1) = N(index+1) + N(index);
                indices2rm = [ indices2rm index];
            end;
        end;
        N(indices2rm)   = [];
        X(indices2rm+1) = [];
        indices2rm = [];
        for index = length(N):-1:2
            if N(index) < 5
                N(index-1) = N(index-1) + N(index);
                indices2rm = [ indices2rm index];
            end;
        end;
        N(indices2rm)   = [];
        X(indices2rm)   = [];
        
        % compute expected values
        % -----------------------        
        for index = 1:length(X)-1
            p1 = rsget( l, X(index+1));
            p2 = rsget( l, X(index  ));
            expect(index) = length(x)*(p1-p2); 
        end;
        
        % value of X2
        % -----------
        res = sum(((expect - N).^2)./expect);
        
        % plot fit
        % --------
        if plotflag
            if plotflag ~= 2, figure('paperpositionmode', 'auto'); end;
            hist(x, 10);
            
            % plot fit
            % --------
            xdiff = X(end)-X(1);
            abscisia   = linspace(X(1)-0.2*xdiff, X(end)+0.2*xdiff, 100);
            %abscisia  = (X(1:end-1)+X(2:end))/2;
            expectplot = zeros(1,length(abscisia)-1);
            for index = 2:length(abscisia); 
                p1 = rsget( l, abscisia(index-1));
                p2 = rsget( l, abscisia(index  ));
                expectplot(index-1) = length(x)*(p2-p1); 
                % have to do this subtraction since this a cumulate density distribution
            end;
            abscisia = (abscisia(2:end)+abscisia(1:end-1))/2;
            hold on; plot(abscisia, expectplot, 'r');
        
            % plot PDF
            % ----------
            pval = linspace(0,1, 102); pval(1) = []; pval(end) = [];
            rp   = l(1) + (pval.^l(3) - (1-pval).^l(4))/l(2);
            fp   = l(2)*1./(l(3).*(pval.^(l(3)-1)) + l(4).*((1-pval).^(l(4)-1)));
            [maxval index]  = max(expect);
            [tmp closestind] = min(abs(rp - abscisia(index)));
            fp = fp./fp(closestind)*maxval;
            plot(rp, fp, 'g');
            legend('Chi2 fit (some bins have been grouped)', 'Pdf', 'Data histogram'  );
            xlabel('Bins');
            ylabel('# of data point per bin');            
            title (sprintf('Fit of distribution using Ramberg-Schmeiser distribution (Chi2 = %2.4g)', res));            
        end;
    end;
    return
