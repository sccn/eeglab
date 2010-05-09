% timefrq() - progressive Power Spectral Density estimates on a single 
%             EEG channel using out-of-bounds and muscle activity rejection 
%             tests. Uses Matlab FFT-based psd().
% Usage: 
%   >> [Power,frqs,times,rejections] = timefrq(data,srate,subwindow);
%   >> [Power,frqs,times,rejections] = ...
%                        timefrq(data,subwindow,fftwindow,substep, ...
%                          epochstep,overlap,srate,nfreqs, ...
%                             rejthresh,minmuscle,maxmuscle,musthresh);
%
% Inputs:
%       data        = single-channel (1,frames) EEG data      {none}
%       srate       = data sampling rate (Hz)                 {256 Hz}
%       subwindow   = subepoch data length per psd()          {<=256}
%
%       fftwindow   = subepoch FFT window length after zero-padding 
%                            (determines freq bin width)      {subwindow}
%       substep     = subepoch step interval in frames        {subwindow/4}
%       epochstep   = output epoch step in frames             {subwindow*2}
%       overlap     = overlap between output epochs in frames {subwindow*2}
%                     total epoch length is (overlap+epochstep)
%       nfreqs      = nfreqs to output (2:nfreqs+1), no DC    {fftwindow/4}
%       rejthresh   = abs() rejection threshold for subepochs {off}
%                     If in (0,1) == percentage of data to reject; else 
%                     reject subepochs reaching > the given abs value.
%       minmuscle   = lower bound of muscle band (Hz)         {30 Hz}
%       maxmuscle   = upper bound of muscle band (Hz)         {50 Hz}
%       musthresh   = mean muscle-band power rejection threshold 
%                            (no percentile option)           {off}
%
% Note: frequency of resulting rejections in tty output 
%       ('o'=out-of-bounds rejection;'+' = muscle band rejection)
%
% Outputs: 
%            Power - time-frequency transform of data (nfreqs,data_epochs)
%            frqs  - frequency bin centers (in Hz) [DC bin not returned]
%            times - midpoints of the output analysis epochs (in sec.)
%            rejections - 8-element vector of rejection statistics =
% [rejthresh,musthresh,goodepochs,badepochs,subacc,subrej,oobrej,musrej]
%            *  rejthresh = out-of-bounds abs rejection threshold 
%            *  musthresh = muscle-band power rejection threshold 
%            *  goodepochs,badepochs = numbers of epochs accepted/rejected
%            *  subacc,subrej = numbers of subepochs accepted/rejected
%            *  oobrej,musrej = numbers of subepochs rejected for oob/muscle
%
% Authors: Tzyy-Ping Jung & Scott Makeig, SCCN/INC/UCSD, La Jolla, 10/1/97 
%
% See also: timef()
 
% Copyright (C) 10/1/97 Tzyy-Ping Jung & Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 11-17-97 added rejections vector output -sm
% 01-25-02 reformated help & license, added links -ad 

function [Power,frqs,times,rejections] = timefrq(data,subwindow,fftwindow,substep,epochstep,overlap,srate,nfreqs,rejthresh,minmuscle,maxmuscle,musclethresh);

MINMUSCLE = 30; % Hz
MAXMUSCLE = 50; % Hz
DEFAULT_SRATE = 256; % Hz
MIN_SUBWINDOW = 8;
MAX_SUBWINDOW = 256;
MINSUBEPOCHS = 3;  % min number of non-rejected subepochs to median-average
SURFPLOT = 1;      % flag to make a surf(|) plot of the time*freq results
HISTPLOT = 0;      % flag to plot the EEG histogram
OOB      = 1e22;   % out-of-bounds large number

if nargin<1
  help timefrq
  return
end
[chans,frames] = size(data);
if chans>1,
  fprintf('timefrq(): data must be one-channel.\n');
  help timefrq
  return
end

if nargin < 12
  musclethresh = 0;
end
if musclethresh==0,
  musclethresh = OOB; % DEFAULT
end

if nargin < 11
  maxmuscle = 0;
end
if maxmuscle==0,
  maxmuscle = MAXMUSCLE; % DEFAULT
end

if nargin < 10
  minmuscle = 0;
end
if minmuscle==0,
  minmuscle = MINMUSCLE; % DEFAULT
end

if nargin < 9
  rejthresh = 0;
end
if rejthresh==0,
  rejthresh = OOB; % DEFAULT
end

if nargin < 8
  nfreqs = 0;
end
if nargin < 7
  srate = 0;
end
if srate==0,
  srate = DEFAULT_SRATE;
end

if nargin < 6
  overlap = 0;
end
if nargin < 5,
  epochstep =0;
end
if nargin < 4,
  substep =0;
end
if nargin < 3
  fftwindow =0;
end
if nargin < 2,
  subwindow = 0;
end
if subwindow==0,
  subwindow = 2^round((log(frames/64)/log(2))); % DEFAULT
  if subwindow > MAX_SUBWINDOW,
     fprintf('timefrq() - reducing subwindow length to %d.\n',subwindow);
     subwindow = MAX_SUBWINDOW;
  end
end
if subwindow < MIN_SUBWINDOW
  fprintf('timefrq() - subwindow length (%d) too short.\n',subwindow);
  return
end

if fftwindow==0,
  fftwindow=subwindow; % DEFAULT
end
if fftwindow < subwindow
  fprintf('timefrq() - fftwindow length (%d) too short.\n',fftwindow);
  return
end
  
if substep==0,
  substep=round(subwindow/4); % DEFAULT
end
if substep < 1
  fprintf('timefrq() - substep length (%d) too short.\n',substep);
  return
end
  
if epochstep==0,
  epochstep=subwindow*2; % DEFAULT
end
if epochstep < 1
  fprintf('timefrq() - epochstep length (%d) too short.\n',epochstep);
  return
end
  
if overlap==0,
  overlap=subwindow*2; % DEFAULT
end
if overlap > epochstep
  fprintf('timefrq() - overlap (%d) too large.\n',overlap);
  return
end
  
if nfreqs==0,
  nfreqs=floor(fftwindow/2); % DEFAULT
end
if nfreqs > floor(fftwindow/2)
  fprintf('timefrq() - nfreqs (%d) too large.\n',nfreqs);
  return
end
%
%%%%%%%%%%%%%%% Compute rejection threshold from percentile %%%%%%%%%%%%%%%
%
if rejthresh > 0 & rejthresh < 1.0,
disp yes
    data = data - mean(data); % make data mean-zero
    fprintf('Sorting mean-zeroed data to compute rejection threshold...\n');
    sortdat = sort(abs(data));
    rejpc = 1.0-rejthresh;
    idx = max(1,round(rejpc*frames));
    rejthresh = sortdat(idx);
    titl = [ 'Out-of-bounds rejection threshold is +/-' ...
                   num2str(rejthresh) ...
                      ' (' num2str(100*rejpc) ' %ile)'];
  if HISTPLOT
  %
  %%%%%%%%%%%%%% Plot data histogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
    figure;        % plot figure for reference
    hist(data,50); % histogram with 50 bins
    hold on; 
    ax = axis; % get current axis limits
    axis([ax(1) ax(2) ax(3) frames/50]); % see side bins clearly
    plot([rejthresh rejthresh],[0 1e10],'g');
    plot([-rejthresh -rejthresh],[0 1e10],'g');
    title(titl);
  else
    fprintf('%s\n',titl);
    fprintf('timefrq(): data histogram plotting disabled.\n');
  end
end
if nfreqs > fftwindow/2
   fprintf('fftwindow of %d will output only %d frequencies.\n',...
                  fftwindow,fftwindow/2);
   return
elseif nfreqs<1
   help timefrq
   fprintf('Number of output frequencies must be >=1.\n');
   return
end
Power=zeros(nfreqs+1,floor(frames/epochstep)); 
                                   % don't use a final partial epoch
badepochs=0;   % epochs rejected
totacc = 0;    % subepochs accepted
totrej = 0;    % subepochs rejected
musrej = 0;    % subepochs rejected for muscle artifact
oobrej = 0;    % subepochs rejected for out-of-bounds values
firstpsd = 1;  % logical variable
zcol = zeros(fftwindow/2+1,1); % column of zeros
%
%%%%%%%%%%%%%%%%%%%%%%%% Print header info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
start = 1;
stop=min(epochstep,frames);

fprintf('\nMoving power spectrum will use epochs of %d frames\n',...
             epochstep+overlap)
fprintf('output at intervals of %d frames.\n',epochstep);
fprintf('Each epoch is composed of %d subepochs of %d frames\n',...
             floor((epochstep+overlap+1-subwindow)/substep), subwindow);
fprintf('starting at %d-frame intervals', substep);
if fftwindow>subwindow,
   fprintf(' and zero-padded to %d frames\n',fftwindow);
else
   fprintf('.\n')
end
fprintf(...
'Rejection criteria: "o" out-of-bounds (>%g), "+" muscle activity (>%g)\n',...
             rejthresh,musclethresh);
%
%%%%%%%%%%%%%%%%%%%%% Process data epochs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
times = zeros(1,1+floor((frames-subwindow)/epochstep));
frqs = zeros(1,nfreqs+1); % initialize in case no valid psd returns
ptmpzeros=zeros(fftwindow/2+1,ceil((stop-start+1)/substep));
epoch = 0;       % epoch counter

for I=1:epochstep:frames-subwindow ;  % for each epoch . . .
  epoch=epoch+1;
  start=max(1,I-overlap);
  stop=min(I+epochstep+overlap-1,frames);
  times(epoch) = round((stop+start)/2); % midpoint of data epoch
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%% Process subepochs %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  %
  tmp=data(start:stop); % copy data subepoch
  subepoch=0; % subepoch counter
  rej = 0; % initialize subepoch rejections counter
  ptmp=ptmpzeros; % start with zeros

  for j=1:substep:stop-start+1-subwindow , % for each subepoch . . .
    subepoch=subepoch+1;

    % idx=find(abs(tmp(j:min(j+subwindow-1,stop-start+1))) <= rejthresh);
    % if length(idx) > subwindow/2 | subepoch == 1

    idx=find(abs(tmp(j:min(j+subwindow-1,stop-start+1))) > rejthresh);
    if length(idx) == 0                % If no point in subepoch out of bounds
      datwin = tmp(j:j+subwindow-1)-mean(tmp(j:j+subwindow-1));
                                       % get 1 power est. for subepoch.
      [ptmp(:,subepoch),frqs]=psd(datwin,fftwindow,srate,hanning(subwindow),0);

      if firstpsd > 0                  % On very first subepoch
         muscle = find(frqs>=minmuscle & frqs<=maxmuscle);
         firstpsd = 0;                 % compute muscle band frequency bins.
      end
      if mean(ptmp(muscle,subepoch))>musclethresh 
         rej = rej+1; fprintf('+')     % If muscle-band out of bounds 
         musrej = musrej+1;
         if subepoch > 1                       
            ptmp(:,subepoch)=ptmp(:,subepoch-1); % reject for muscle noise.
         else
            ptmp(:,subepoch)=zcol;     % set back to zeros
         end
      end
    else                               % some data value out of bounds
      rej = rej+1;fprintf('o')         % so reject for out of bounds
      oobrej = oobrej+1;
      if subepoch > 1
         ptmp(:,subepoch)=ptmp(:,subepoch-1);        
      end
    end
  end                 % end of subepochs

  k = 1;
  while ptmp(:,k) == zcol   % while subepoch power values all zeros . . .
   % sum(ptmp(:,k))
     k = k+1;
     if k >  subepoch, 
         break
     end
  end
  if rej>0,
     fprintf('\n');
     fprintf(' epoch %d: %d of %d subepochs rejected ',epoch,rej,subepoch);
     totrej = totrej+rej;
     totacc = totacc+(subepoch-rej);
  else
     fprintf('.');
  end
  if k>subepoch
     fprintf('(no non-zero subepochs)\n',k);
  elseif k>1
     fprintf('(first non-zero subepoch %d)\n',k);
  elseif rej>0
     fprintf('\n')
  end
  if k<= subepoch & subepoch-rej >= MINSUBEPOCHS
     Power(:,epoch)=[median(ptmp(1:nfreqs+1,k:subepoch)')]'; 
                                          % omit initial zero cols
  elseif epoch > 1
     Power(:,epoch) = Power(:,epoch-1);
     fprintf ('')
     badepochs = badepochs+1;
  end
end % end epochs

Power = Power(2:nfreqs+1,:);  % omit DC power bin
frqs = frqs(2:nfreqs+1);    % omit DC power bin
times = times/srate;  % convert to seconds
if length(times)>1
       time_interval = times(2)-times(1);
else
       time_interval = 0;
end
if nargout>3,
   rejections = [rejthresh,musclethresh,epoch-badepochs,badepochs,...
                        totacc,totrej,oobrej,musrej]; 
end
if epoch<1
  fprintf('No epochs processed: too little data.\n');
  return
end
%
%%%%%%%%%%%%%%%% Print trailer info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
fprintf('\nTotal of %d epochs processed (%d rejected).\n',...
                         epoch,badepochs);
fprintf('Output is %d freqs by %d epochs.\n',...
                         nfreqs,length(times));
fprintf('Output epoch length %d frames (%g secs).\n', ...
           epochstep+overlap,(epochstep+overlap)/srate);
fprintf('Output interval %d frames (%g secs).\n', ...
                         epochstep,time_interval);
fprintf('First and last time points: %g and %g secs.\n',...
                         times(1),times(length(times)));
%
%%%%%%%%%%%%%%% Make surf() plot of time-frequency distribution %%%%%%
%
if nfreqs>1 & epoch>1 & SURFPLOT 
  if min(min(Power))>0
    fprintf('Plot shows dB log(Power) - Power output is not log scaled.\n');
    off    = [50 -50 0 0];      % successive figure offset in pixels
    pos = get(gcf,'Position');
    figure('Position',pos+off); % make the 2nd plot offset from the 1st
    
    surf(times,frqs,10*log(Power)/log(10))
    view([0 90]); % top view
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    shading interp
    title('timefrq()');
    c=colorbar;
    t=axes('Position',[0 0 1 1],'Visible','off');
    text(0.85,0.08,'dB','Parent',t);
  else
    fprintf('Some or all output Power estimates were zero - too many rejections?\n');
  end
end
