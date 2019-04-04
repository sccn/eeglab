% phasecoher() - Implements inter-trial amp/coherence using Gaussian wavelets.
%                Returns same data length as input frames.
%                Plots results when nargin>6. Outputs have flat ends 
%                at data indices [1:halfwin] and [frames-halfwin:frames].
% Usage:
%     >> [amps,cohers ] = phasecoher(data,frames,srate,freq,cycles); 
%     >> [amps,cohers,cohsig,ampsig,allamps,allphs] ...
%                 = phasecoher(data,frames,...
%                                srate,freq,cycles,...
%                                  alpha,times,titl,...
%                                    warpframes, events);
% Inputs:
%   data   = input data, (1,frames*trials) or NB: (frames,trials) 
%   frames = frames per trial
%   srate  = sampling rate (in Hz)
%   freq   = frequency to work on (in Hz)
%   cycles = cycles in Gaussian wavelet window (float) {3}
%   alpha  = (0 0.1] significance probability threshold. Requires 
%            >=3 output arguments. alpha=0 -> no signif {default: 0}.
%   times  = vector of latencies (in ms) for plotting {default: no plot}
%   titl   = [string] plot title {default none}
%   warpframes = frame numbers of warped events (below)
%   events = matrix of events in each trial, size (nevents, trials)
%            as frame numbers.
% Outputs:
%   amps    = mean amplitude at each time point
%   cohers  = phase coherence at each time point [0,1]
%   cohsig  = coherence significance threshold (bootstrap, alpha level)
%   ampsig  = amplitude significance thresholds [lo high] (bootstrap, alpha level)
%   allamps = amplitudes at each trial and time point (frames,trials)
%   allphs  = phase (deg) at each trial and time point (frames,trials)
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 5-5-98 
%
% See also: erpimage()

% Copyright (C) 5-5-98 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 5-7-98 added frames, made input data format one-row -sm
% 5-8-98 added MIN_AMP, times, PLOT_IT -sm
% 10-27-98 added cohsig, alpha -sm
% 07-24-99 added allamps -sm
% 02-29-00 added ampsig -sm
% 12-05-00 changed complex abs() to sqrt( .^2+ .^2) to avoid possible i ambiguity -sm & tpj
% 12-31-00 added ...20...40... frames fprinting -sm
% 08-17-01 added allphs option -sm
% 08-18-01 debugged cohsig plotting line (302) -sm
% 01-25-02 reformated help & license, added links -ad 

function [amps,cohers,cohsig,ampsig,allamps,allphs] = ...
      phasecoher(data, frames, srate, freq, cycles, alpha, times, titl, timeStretchRef, timeStretchMarks)

MIN_AMP = 10^-6;
DEFAULT_ALPHA = nan;% no bootstrap computed
COHSIG_REPS = 500;  % default number of bootstrap surrogate coherence values 
                    % to compute 
TITLEFONT= 18;
TEXTFONT = 16;
TICKFONT = 14;
DEFAULT_CYCLES=3;
ampsig = []; % initialize for null output
cohsig = [];

if nargin<4
  help phasecoher
  return
end
if nargin<5
  cycles = DEFAULT_CYCLES;
end

if nargin < 8 || isempty(titl)
   titl = '';
end

if nargin < 7 || isempty(titl) || isempty(times)
  PLOT_IT = 0;
elseif length(times) ~= frames
  fprintf('phasecoher(): times vector length must be same as frames.\n')
  return
else
  PLOT_IT = 1;
end
 
if nargin < 6
    alpha = nan; % no alpha given
end
if nargout > 2 && isnan(alpha) % if still want cohsig
  alpha = DEFAULT_ALPHA; 
elseif nargout > 2
  if alpha < 0 || alpha > 0.1
    help phasecoher
    fprintf('phasecoher(): alpha out of bounds.\n');
    return
  end
  if alpha==0
    alpha = nan; % no bootstrap
  end
  if alpha*COHSIG_REPS < 5
    COHSIG_REPS = ceil(5/alpha);
    fprintf('   Computing %d bootstrap replications.\n',COHSIG_REPS);
  end
elseif ~isnan(alpha)
  alpha = nan; % no cohsig calculation
end

if length(frames)>1
   help phasecoher
   fprintf('phasecoher(): frames should be a single integer.\n');
   return
end
if frames == 0,
   frames = size(data,1);
end
trials = size(data,1)*size(data,2)/frames;
if floor(trials) ~= trials
   fprintf('phasecoher(): data length not divisible by %d frames.\n',frames);
   return
end

fprintf('phasecoher(): Analyzing %d data trials of %d frames ',trials,frames);
if trials < 10
  fprintf(...
  'Low number of trials (%d) may not give accurate coherences.\n',trials)
end

if size(freq,1)*size(freq,2)~=1
   fprintf('\nphasecoher(): only one frequency can be analyzed at a time.\n');
   help phasecoher
   return
end

if size(data,1) == 1
  data = reshape(data,frames,trials); % trials are columns
end
window = gauss(ceil(srate/freq*cycles),2); % gauss(std,+/-stds)
winlength = length(window);
halfwin = floor(winlength/2);
fprintf('\n  with a moving %d-frame analysis window...',winlength);

if frames < winlength
  fprintf(...
  '\nProblem: Epoch length (%d frames) too short for analysis with %g cycles!\n',...
                frames,                                 cycles);
  return
end

%
% Extend the data to minimize edge effects
%
data = [data([halfwin+1:-1:1],:); ...
                data; ...
        data([frames:-1:frames+1-(winlength-halfwin)],:)];
%
% Remove epoch means
%
%data = data - ones(frames+winlength+1,1)* mean(data); % remove column means

angleinc = cycles*2*pi/winlength;
cosx = cos(-cycles*pi:angleinc:cycles*pi); % sinusoids
cosx = cosx(1:winlength);
sinx = sin(-cycles*pi:angleinc:cycles*pi);
sinx = sinx(1:winlength);
coswin = window.*cosx;               % window sinusoids
sinwin = window.*sinx;
coswin = coswin/(coswin*cosx');      % normalize windowed sinusoids
sinwin = sinwin/(sinwin*sinx');

% figure;plot(coswin,'r');hold on; plot(sinwin,'b');
% iang = -cycles*pi:angleinc:cycles*pi;
% iang = iang(1:winlength);
% figure;plot(iang,[sinwin;coswin]);
amps    = zeros(1,frames);
if nargout > 3
  allamps = zeros(frames,trials);
end
if nargout > 5
  allphs = zeros(frames,trials);
end
cohers  = zeros(1,frames);
ix = 0:winlength-1;
% nsums = zeros(1,frames);  % never called

C = [];

for f = 1:frames %%%%%%%%%%%%%%% frames %%%%%%%%%%%%%%%%%%%%
  epoch = data(ix+f,:);
  epoch = epoch - ones(winlength,1)*mean(epoch); % remove epoch means
  if rem(f,50)== 0
    fprintf(' %d',f)
  end
  for t = 1:trials  %%%%%%%%%%%%%%% trials %%%%%%%%%%%%%%%%%%%
    realpart = coswin*epoch(:,t);
    imagpart = sinwin*epoch(:,t);
    C(f,t) = complex(realpart, imagpart);
  end
end

allamps = sqrt(C.*conj(C)); %compute all amplitudes for all frames, all trials
allphs = angle(C); %get the phase

if exist('timeStretchRef') && exist('timeStretchMarks') && ...
    length(timeStretchRef) > 0 & length(timeStretchMarks) > 0 %Added -Jean
  for t=1:trials
    M = timewarp(timeStretchMarks(:,t)', timeStretchRef');
    allamps(:,t) = M*allamps(:,t);
    allphs(:,t) = angtimewarp(timeStretchMarks(:,t)', timeStretchRef', ...
                              allphs(:,t));
  end
end

[amps, cohers, nsums]=getAmpCoh(allamps, allphs, MIN_AMP);

% Old routine, for archeological purposes
% $$$ realcoh = zeros(1,frames);
% $$$ imagcoh = zeros(1,frames);
% $$$ for f = 1:frames %%%%%%%%%%%%%%% frames %%%%%%%%%%%%%%%%%%%%
% $$$   epoch = data(ix+f,:);
% $$$   %epoch = epoch - ones(winlength,1)*mean(epoch); % remove epoch means
% $$$   if rem(f,50)== 0
% $$$     fprintf(' %d',f)
% $$$   end
% $$$   for t = 1:trials  %%%%%%%%%%%%%%% trials %%%%%%%%%%%%%%%%%%%
% $$$     realpart = coswin*epoch(:,t);
% $$$     imagpart = sinwin*epoch(:,t);
% $$$     amp = sqrt(realpart.*realpart+imagpart.*imagpart);
% $$$     if amp >= MIN_AMP
% $$$       amps(f) = amps(f) + amp; % sum of amps
% $$$       realcoh(f) = realcoh(f) + realpart/amp;
% $$$       imagcoh(f) = imagcoh(f) + imagpart/amp;
% $$$       nsums(f) = nsums(f)+1;
% $$$     end
% $$$     if nargout > 3 
% $$$       if amp < MIN_AMP
% $$$         amp = MIN_AMP;
% $$$       end
% $$$       allamps(f,t) = amp;
% $$$     end
% $$$     if nargout > 5
% $$$       allphs(f,t) = 180/pi*angle(realpart+i*imagpart);
% $$$     end
% $$$   end
% $$$   if nsums(f)>0
% $$$     amps(f) = amps(f)/nsums(f);
% $$$     realcoh(f) = realcoh(f)/nsums(f);
% $$$     imagcoh(f) = imagcoh(f)/nsums(f);
% $$$   else
% $$$     amps(f) = 0;
% $$$     realcoh(f) = 0;
% $$$     imagcoh(f) = 0;
% $$$   end
% $$$ end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ cohers = sqrt(realcoh.^2+imagcoh.^2);

fprintf('\n');




cohsig = [];

if ~isnan(alpha)  %%%%%%%%%%%%%% Compute cohsig/ampsig %%%%%%%%%%%%%%
 % ix = 0:winlength-1;  % never called
 bootcoher = zeros(1,COHSIG_REPS);
 bootamp   = zeros(1,COHSIG_REPS);
 
 bootallamps = zeros(COHSIG_REPS, trials); %Added -Jean
 bootallphs = zeros(COHSIG_REPS, trials); %Added -Jean
 
 fprintf('Computing %d bootstrap coherence values... ',COHSIG_REPS); 
 for f = 1:COHSIG_REPS %%%%%%%%%%%%%%% Bootstrap replications %%%%%%%%%%%
  if rem(f,50) == 0
    fprintf('%d ',f);
  end
  randoff = floor(rand(1,trials)*(frames-winlength))+1; % random offsets
  
  %Create randomized time-stretched allamps and allphs arrays (see above)
  for t = 1:trials
    bootallamps(f, t) = allamps(randoff(t), t);
    bootallphs(f, t) = allphs(randoff(t), t);
  end
 end

 [bootamp, bootcoher]=getAmpCoh(bootallamps, bootallphs, MIN_AMP);
 
 fprintf('\n');

 bootcoher = sort(bootcoher); % sort low to high
 cohsig = bootcoher(round(COHSIG_REPS*(1.0-alpha)));

 bootamp = sort(bootamp); % sort low to high
 ampsig = [bootamp(round(COHSIG_REPS*(alpha))) ...
           bootamp(round(COHSIG_REPS*(1.0-alpha)))];
% keyboard

end %%%%%%%%%%%%%%%%%%%%%%%%%%%% end cohsig %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for f=1:halfwin                 % pad amps, cohers to front of input data
  amps(f) = amps(halfwin+1);
  cohers(f) = cohers(halfwin+1);
end
for f=frames:-1:frames-halfwin  % pad amps, cohers to end of input data
  amps(f) = amps(frames-halfwin);
  cohers(f) = cohers(frames-halfwin);
end

if PLOT_IT %%%%%%%%%%%%%% make two-panel plot of results %%%%%%%%

  subplot(2,1,1);plot(times,amps');
  title(titl,'fontsize',TITLEFONT,'fontweight','bold');
     ylabel(['Amplitude (' num2str(freq) ' Hz)'],...
           'fontsize',TEXTFONT,'fontweight','bold');
     ax = axis; 
     hold on; plot([0 0],[0 1],'k');   % vertical line at time 0
     axis([ax(1) ax(2) ax(3) ax(4)*1.25]);
     set(gca,'FontSize',TICKFONT);
     set(gca,'FontWeight','bold');

  subplot(2,1,2);plot(times,cohers','r');
     ylabel(['Intertrial Coherence (' num2str(freq) ' Hz)'],...
           'fontsize',TEXTFONT,'fontweight','bold');
     xlabel('Time (ms)','fontsize',TEXTFONT,'fontweight','bold');
     hold on
  winstframe = floor(frames/7);
  winframes = [winstframe:winstframe+winlength-1];
  wintimes = times(winframes);
  ax = axis; 
  plot(wintimes,0.8+window*0.1,'k');
  plot(wintimes,0.8-window*0.1,'k');
  % ax2 = axis;     % never called
  hold on; plot([0 0],[0 1000],'k');   % vertical line at time 0
  axis([ax(1) ax(2) 0 1]);
     set(gca,'fontSize',TICKFONT);
     set(gca,'FontWeight','bold');
     alpha
  if ~isnan(alpha)                      % plot coher significance
    plot([wintimes(1) wintimes(end)],[cohsig cohsig],'r'); 
       % was [times(1) times(winframes)] !??
  end                                  
end



function [amps, cohers, nsums] = getAmpCoh(allamps, allphs, MIN_AMP)
  minampfilter = allamps >= MIN_AMP;
  nsums = sum(minampfilter,2);

  amps(find(nsums == 0)) = 0; %zero the amplitude if no trial shows
                              %significant power at that frame
  cohers(find(nsums == 0)) = 0; %zero the coherence too if no trial shows
                                %significant power at that frame
                                %Now average out amplitudes over trials
  % allminamps is never used. TF 04/02/2007                              
  %allminamps = allamps;
  % nargout is never greater than 3. Bug 262. TF 04/02/2007
  %if nargout > 3
  %  allminamps(~minampfilter) = MIN_AMP;
  %end
  allzeramps = allamps .* minampfilter;
  allzeramps = allzeramps(find(nsums ~= 0),:);
  amps(find(nsums ~= 0)) = sum(allzeramps,2) ./ nsums(find(nsums ~= 0));

  %Convert angles to complex for summing
  allzerphs = complex(cos(allphs), sin(allphs)) .* minampfilter;
  allzerphs = allzerphs(find(nsums ~= 0), :);
  cohers(find(nsums ~= 0)) = sum(allzerphs,2) ./ nsums(find(nsums ~= 0));
  cohers = sqrt(cohers .* conj(cohers));

  function outvec = gauss(frames,sds)

outvec = [];
if nargin < 2
  help gauss
  return
end
if sds <=0 || frames < 1
  help gauss
  return
end

incr = 2*sds/(frames-1);
outvec = exp(-(-sds:incr:sds).^2);
