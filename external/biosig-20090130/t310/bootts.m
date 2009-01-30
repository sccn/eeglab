function [ts_boot, cfl, cfu] = bootts(ts, B, alpha, med)
% Calculates non-parametric t-percentile bootstrap statistics.
%
% Calculates non-parametric t-percentile bootstrap statistics of a given time
% series.
%
% Input parameters:
%   ts    ... Time series (trials x epochs)
%   B     ... Number of resamplings (default: 300)
%   alpha ... Alpha significance of confidence intervals (default: [0.1 0.05 0.01])
%   med   ... 0: mean, 1: median (default: 0)
% 
% Output parameters:
%   ts_boot ... Bootstrapped time series
%   cfl     ... Lower limit of confidence interval
%   cfu     ... Upper limit of confidence interval

% Copyright by Bernhard Graimann, modified by Clemens Brunner
% $Revision: 1.1 $ $Date: 2009-01-30 06:04:44 $
% E-Mail: clemens.brunner@tugraz.at

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


if nargin < 1, help bootts; return; end;

if nargin < 2, B = 300; end;
if nargin < 3, alpha = [0.1 0.05 0.01]; end;
if nargin < 4, med = 0; end;


%==============================================================================
% Initialization
%==============================================================================

[trials,tpts] = size(ts);    % trials and number of samples (time points)

nbootlevels = length(alpha);

cfl = zeros(nbootlevels,tpts);
cfu = zeros(nbootlevels,tpts);

ts_boot = zeros(1,tpts);     % bootstrapped time series
midx = fix(B/2);             % index of median value

ts_mean = mean(ts);          % mean over all trials
ts_std = std(ts);            % standard deviation over all trials


%==============================================================================
% Calculate percentiles
%==============================================================================

for k=1:nbootlevels
  q1(k)=floor(B*alpha(k)/2);
  if q1(k)<1
    q1(k)=1;
    B=ceil(2/alpha(k));
    fprintf('BOOTTS: At least %d bootstrap resamplings are necessary for alpha=%.3f\n',B,alpha(k));
  end;
  q2(k)=B-q1(k)+1;
end;


%==============================================================================
% Bootstrap each single sample
%==============================================================================

for j=1:tpts

  % do the resampling
  idx = ceil(trials*rand(trials,B));
  strials = ts(:,j);
  resamples = strials(idx);               % draw B resamples  (trials x boostraps)
  resamples_stdest = std(resamples);      % standard deviation of resamples (1 x bootstraps)

  % average the resampled values
  if med == 0, resamples = mean(resamples); else, resamples = median(resamples); end;

  bootstat = sort((resamples-ts_mean(j))./resamples_stdest);
  ts_boot(j) = bootstat(midx)*resamples_stdest(midx)+ts_mean(j);

  % for each bootstrap level, i.e. percentile, according to alpha, calc. conf. limits
  for kk = 1:nbootlevels
    cfl(kk,j) = ts_mean(j)-bootstat(q2(kk))*ts_std(j);
    cfu(kk,j) = ts_mean(j)-bootstat(q1(kk))*ts_std(j);
  end;

end;
