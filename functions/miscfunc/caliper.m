% caliper() - Measure a set of spatial components of a given data epoch relative to 
%           a reference epoch and decomposition. 
% Usage: 
%   >> [amp,window]=caliper(newepoch,refepoch,weights,compnums,filtnums,times,'noplot');
%
% Inputs:
%     newepoch = (nchannels,ntimes) new data epoch
%     refepoch = a (nchannels,ntimes) reference data epoch
%     weights  = (nchannels,ncomponents) unmixing matrix (e.g., ICA weights*sphere)
%     compnums = vector of component numbers to return amplitudes for {def|0: all}
%     filtnums = [srate highpass lowpass] filter limits for refepoch {def|0: allpass}
%     times    = [start_ms end_ms] epoch latency limits, else latencies vector {def|0: 0:EEG.pnts-1}
%     'noplot' = produce no plots {default: plots windows for the first <=3 components}
%
% Outputs:
%     amps = (1,length(compnums)) vector of mean signed component rms amplitudes 
%     windows = (length(compnums)),length(times)) matrix of tapering windows used
%
% Notes:
%   Function caliper() works as follows: First the reference epoch is decomposed using 
%   the given weight matrix (may be ICA, PCA, or etc). Next, the time course of the 
%   main lobe of the activation in the reference epoch (from max to 1st min, forward 
%   and backward in time from abs max, optionally after bandpass filtering) is used 
%   to window the new epoch. Then, the windowed new epoch is decomposed by the same 
%   weight matrix, and signed rms amplitude (across the channels) is returned of the 
%   projection of each of the specified component numbers integrated across the windowed 
%   epoch. If not otherwise specified, plots the windows for the first <= 3 components listed.
%
% Example: Given a grand mean response epoch and weight matrix, it can be used 
%   to measure amps of grand mean components in single-subject averages.
%
% Authors: Scott Makeig & Marissa Westerfield, SCCN/INC/UCSD, La Jolla, 11/2000 

% Copyright (C) 11/2000 Scott Makeig & Marissa Westerfield,, SCCN/INC/UCSD 
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

% Edit History:
% 12/05/00 -- added fig showing data, ref activation, and window vector -mw
% 12/19/00 -- adjusted new icaproj() args -sm
% 01-25-02 reformated help & license -ad 

function [amps,windows] = caliper(newepoch,refepoch,weights,compnums,filtnums,times,noplot)

if nargin < 3
  help caliper
  return
end

if nargin<4
  compnums = 0;
end

nchans = size(newepoch,1);
ntimes = size(newepoch,2);

if compnums(1) == 0 | isempty(compnums(1)) 
   compnums = 1:nchans;
end

if min(compnums) < 1 | max(compnums) > size(weights,2)
   help caliper
   return
end

if nargin<5
  filtnums = [];
else
  if isempty(filtnums)
      filtnums = [];
  elseif length(filtnums)==1 & filtnums(1)==0
     filtnums = [];
  elseif length(filtnums) ~= 3
     fprintf('\ncaliper(): filter parameters (filtnums) must have length 3.\n')
     return
  end
end

if nargin< 6  | isempty(times) | (length(times)==1 & times(1)==0)
  times = 0:ntimes-1;
else
  if length(times) ~= ntimes
    if length(times) ~= 2
      fprintf('caliper(): times argument should be [startms endms] or vector.\n')
      return
    else
      times = times(1):(times(2)-times(1))/(ntimes-1):times(2);
      times = times(1:ntimes);
    end
  end
end

if nargin < 7
   noplot = 0;
else
   noplot = 1;
end

refact = weights(compnums,:)*refepoch; % size (length(compnums),ntimes)
newact = weights(compnums,:)*newepoch;

if length(filtnums) == 3
  if ~exist('eegfilt')
    fprintf('caliper(): function eegfilt() not found - not filtering refepoch.\n');
  else
   try
      refact = eegfilt(refact,filtnums(1),filtnums(2),filtnums(3));
   catch
      fprintf('\n%s',lasterr)
      return
   end
  end
end

if size(weights,1) == size(weights,2)
   winv = inv(weights);
else
   winv = pinv(weights);
end
winvrms = sqrt(mean(winv.*winv)); % map rms amplitudes
  
amps = [];
windows = [];
n = 1; % component index
for c=compnums
  if floor(c) ~= c
    help caliper
    fprintf('\ncaliper(): component numbers must be integers.\n')
    return
  end

  [lobemax i] = max(abs(refact(n,:)));
  f = i+1;
  oldact = lobemax;
  while f>0 & f<=ntimes & abs(refact(n,f)) < oldact
    oldact = abs(refact(n,f));
    f = f+1;
  end % f is now one past end of "main lobe"
  lobeend = f-1;

  f = i-1;
  oldact = lobemax;
  while f>0 & f<=ntimes & abs(refact(n,f)) < oldact
    oldact = abs(refact(n,f));
    f = f-1;
  end % f is now one past start of "main lobe"
  lobestart = f+1;

  refact(n,1:lobestart) = zeros(1,length(1:lobestart));
  refact(n,lobeend:end) = zeros(1,length(lobeend:ntimes));
  windows = [windows; refact(n,:)];

  refnorm = sum(refact(n,:));
  if abs(refnorm)<1e-25
    fprintf('caliper(): near-zero activation for component %d - returning NaN amp.\n',c)
    amps = [amps NaN];
  else
    refact(n,:) = refact(n,:)/refnorm; % make reference epoch window sum to 1
    amps = [amps winvrms(c)*sum(refact(n,:).*newact(n,:))];
  end
  n = n+1;
  if ~noplot & n <= 4  %%% only plot out at most the first 3 components
     refproj = icaproj(refepoch,weights,c);
     refproj = env(refproj);
     windproj = winv(:,c)*(refact(n-1,:)*refnorm);
     windproj = env(windproj);
     figure; plot(times,newepoch(2:nchans,:),'g');
     hold on;h=plot(times,newepoch(1,:),'g');
     hold on;h=plot(times,newepoch(1,:),'g',...
                    times,refproj(1,:),'b',...
                    times,windproj(1,:),'r','LineWidth',2);
     hold on;plot(times,refproj(2,:),'b',times,windproj(2,:),'r','LineWidth',2);
     set(h(1),'linewidth',1);
     legend(h,'new data','comp. act.','window');
     title(['Component ',int2str(c),';  rms amplitude = ',num2str(amps(n-1))],...
            'FontSize',14);
     ylabel('Potential (uV)');
  end;             %% endif
end % c

