% qqdiagram() - Empirical quantile-quantile diagram.
%               The quantiles (percentiles) of the input distribution Y are plotted (Y-axis)
%               against the corresponding quantiles of the input distribution X.
%               If only X is given, the corresponding quantiles are plotted (Y-axis)
%               against the quantiles of a Gaussian distribution ('Normal plot').
%               Two black dots indicate the lower and upper quartiles.
%               If the data in X and Y belong the same distribution the plot will be linear.
%               This will be true also if the data in X and Y belong to two distributions with
%               the same shape, one distribution being rescaled and shifted with respect to the
%               other.
%               If only X is given, a line is plotted to indicate the mean of X, and a segment
%               is plotted to indicate the standard deviation of X.
%
% Usage:
%   >>  ah  =  qqdiagram( x, y, pk );
%
% Inputs:
%   x       - vector of observations
%
% Optional inputs:
%   y       - second vector of observation to compare the first to
%   pk      - the empirical quantiles will be estimated at the values in pk [0..1]
%
% Author: Luca Finelli, CNL / Salk Institute - SCCN, 20 August 2002
%
% Reference: Stahel W., Statistische Datenanalyse, Vieweg, Braunschweig/Wiesbaden, 1995
%
% See also: 
%   quantile(), signalstat(), eeglab() 

% Copyright (C) 2002 Luca Finelli, Salk/SCCN, La Jolla, CA
%
% Reference: Stahel, W. Statistische Datenanalyse, Vieweg, Braunschweig/Wiesbaden 1995
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

% $Log: not supported by cvs2svn $

function qqdiagram( x , y, pk )

if nargin < 1
	help qqdiagram;
	return;
end;	

if (nargin == 3 & (any(pk > 1) | any(pk < 0)))
    error('qqdiagram(): elements in pk must be between 0 and 1');
end

if nargin==1
	y=x;    
	nn=max(1000,10*length(y))+1;
	x=randn(1,nn);
end

if nargin < 3
	nx=sum(~isnan(x));
	ny=sum(~isnan(y));
   	k=min(nx,ny);
    pk=((1:k) - 0.5) ./ k;  % values to estimate the empirical quantiles at 
else 
    k=length(pk);
end

if nx==k
    xx=sort(x(~isnan(x)));
else
    xx=quantile(x(~isnan(x)),pk);
end

if ny==k
    yy=sort(y(~isnan(y)));
else
    yy=quantile(y(~isnan(y)),pk);
end

% Quartiles
xqrt1=quantile(x,0.25); xqrt3=quantile(x,0.75);
yqrt1=quantile(y,0.25); yqrt3=quantile(y,0.75);

plot(xx,yy,'+')
hold on

% Dawing the line
if nargin ==1
    maxy=max(y);
    miny=min(y);
    rangey=maxy-miny;
	ymin=miny-rangey/50;
	ymax=maxy+rangey/50;
	
	plot([(miny-mean(y))/std(y) (maxy-mean(y))/std(y)],[miny maxy],'r-.')
	xlim=get(gca,'XLim');
	plot([1 1],[ymin  (mean(y)+std(y))],'k--')
	plot([1 1],[mean(y)  (mean(y)+std(y))],'k-','LineWidth',2)
	text(1, mean(y)+3*rangey/50,' St.dev. < X > ')
	plot([0 0],[ymin  mean(y)],'k--')
	plot(xlim,[mean(y) mean(y)],'k--')
    text(xlim(1), mean(y)+rangey/50,' Mean < X > ')
	plot([xqrt1  xqrt3],[yqrt1 yqrt3],'k.','MarkerSize',10)
	set(gca,'YLim',[ymin ymax])
	xlabel('Standard Normal Quantiles [Std.Dev.]')
	ylabel('X Quantiles')
else
	% For normally distributed data, the slope of the plot line is equal to the ratio
	% of the standard deviation of the distributions
	sigma=(yqrt3-yqrt1)/(xqrt3-xqrt1); 
	cx=(xqrt1 + xqrt3)/2
	cy=(yqrt1 + yqrt3)/2
	maxy=cy+sigma*(max(x)-cx);
	miny=cy-sigma*(cx-min(x));
	
	plot([xqrt1 xqrt3],[yqrt1 yqrt3],'k-','LineWidth',2); % IQR range
	plot([min(x) max(x)],[miny maxy],'r-.');
    xlabel('X Quantiles');
    ylabel('Y Quantiles');
end
