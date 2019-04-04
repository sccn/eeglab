% crossf() - Returns estimates and plot of event-related coherence (ERC) changes
%            between data from two input channels. The lower panel gives the
%            coherent phase difference between the processes. In this panel, for Ex.
%               -90 degrees (blue) means xdata leads ydata by a quarter cycle.
%                90 degrees (orange) means ydata leads xdata by a quarter cycle.
%            Click on each subplot to view separately and zoom in/out.
%
% Function description:
%            Uses EITHER fixed-window, zero-padded FFTs (faster) OR constant-Q 
%            0-padded DFTs (better sensitivity), both Hanning-tapered. Output 
%            frequency spacing is the lowest frequency (srate/winsize) divided 
%            by the padratio.
%
%            If number of output arguments > 4, then bootstrap statistics are 
%            computed (from a distribution of 200 (NACCU) surrogate baseline
%            data epochs) for the baseline epoch, and non-significant features 
%            of the output plots are zeroed (e.g., plotted in green). Baseline
%            epoch is all windows with center times < 0 (MAX_BASELN)
%
%            If number of output arguments > 5, coherency angles (lags) at
%            significant coherency (time,frequency) points are plotted as well.
%
% Usage: 
%      >> [coh,mcoh,timesout,freqsout,cohboot,cohangles] = crossf(xdata,ydata,...
%                                              frames,tlimits,titl,          ...
%                                              srate,cycles,winsize,timesout,...
%                                              padratio,maxfreq,alpha,verts);
%
% Inputs:
%       xdata       = first single-channel (1,frames*nepochs) data  {none}
%       ydata       = second single-channel (1,frames*nepochs) data {none}
%       frames      = frames per epoch                        {768}
%       tlimits     = epoch time limits (ms) [mintime maxtime]{-1000 2000}
%       titl        = figure title                            {none}
%       srate       = data sampling rate (Hz)                 {256}
%       cycles      = >0 -> number of cycles in each analysis window (slower)
%                     =0 -> use FFT (constant window length)  {0}
%       winsize     = cycles==0: data subwindow length (2^k<frames)
%                     cycles >0: *longest* window length to use; 
%                     determines the lowest output frequency  {~frames/8}
%       timesout    = number of output times (int<frames-winsize){200}
%       padratio    = FFT-length/winsize (2^k)                {2}
%                     Multiplies the number of output frequencies.
%       maxfreq     = maximum frequency to plot (Hz)          {50}
%       alpha       = Two-tailed bootstrap signif. probability {0.02}
%                     Sets n.s. plotted output values to green (0). 
%                     NOTE that it requires at least FIVE output arguments!
%       verts       = times of vertical lines (other than time 0) {none}
%       caxma       = color axis maximum (magnitude) {default: from data}
%
% Outputs: 
%       coh         = between-channel coherency changes (nfreqs,timesout)
%       mcoh        = vector of mean baseline coherence at each frequency
%       timesout    = vector of output times (subwindow centers) in ms.
%       freqsout    = vector of frequency bin centers in Hz.
%       cohboot     = [2,nfreqs] matrix of [lower;upper] coh significance diffs.
%       cohangle    = coherency angles (nfreqs,timesout) 
%
% Note: when cycles==0, nfreqs is total number of FFT frequencies.
%
% Authors: Sigurd Enghoff & Scott Makeig, SCCN/INC/UCSD, La Jolla, 1998 
%
% See also: timef()

% Copyright (C) 8/1/98 Sigurd Enghoff & Scott Makeig, SCCN/INC/UCSD
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

% 11-20-98 defined LINEWIDTH constant -sm
% 04-01-99 made number of frequencies consistent -se
% 06-29-99 fixed constant-Q freq indexing -se
% 08-13-99 added cohangle plotting -sm
% 08-20-99 made bootstrap more efficient -sm
% 08-24-99 allow nan values introduced by possible eventlock() preproc. -sm
% 03-16-00 added lead/lag interpretation to help msg - sm & eric visser
% 03-16-00 added axcopy() feature -sm & tpj
% 04-20-00 fixed Rangle sign for wavelets, added verts array -sm
% 01-22-01 corrected help msg when nargin<2 -sm & arno delorme
% 01-25-02 reformated help & license, added links -ad 

function [R,mbase,times,freqs,Rboot,Rangle,Rsignif] = crossf(X,Y,epoch,timelim,ftitle,Fs,varwin,winsize,nwin,oversmp,maxfreq,alpha,verts,caxmax)

% Constants set here:
MAX_BASELN      = 0;            % Windows with center times < this are in baseline.
NACCU           = 200;			% Number of sub-windows to accumulate
if nargin>13
  COH_CAXIS_LIMIT = caxmax;
else
  COH_CAXIS_LIMIT = 0;          % 0 -> use data limits; else positive value
end
                                % giving symmetric +/- caxis limits.
AXES_FONT       = 10;
LINEWIDTH       = 2;
TITLE_FONT      = 8;
ANGLEUNITS      = 'deg';        % angle plotting units - 'ms' or 'deg'

% Commandline arg defaults:
DEFAULT_EPOCH	= 768;			% Frames per epoch
DEFAULT_TIMELIM = [-1000 2000];	% Time range of epochs (ms)
DEFAULT_FS		= 256;			% Sampling frequency (Hz)
DEFAULT_NWIN	= 200;			% Number of windows = horizontal resolution
DEFAULT_VARWIN	= 0;			% Fixed window length or base on cycles.
								% =0: fix window length to nwin
								% >0: set window length equal varwin cycles
								%     bounded above by winsize, also determines
								%     the min. freq. to be computed.
DEFAULT_OVERSMP	= 2;			% Number of times to oversample = vertical resolution
DEFAULT_MAXFREQ = 50;			% Maximum frequency to display (Hz)
DEFAULT_TITLE	= '';			% Figure title
DEFAULT_ALPHA   = 0.02;			% Default two-sided significance probability threshold
MARGIN          = 0.12;         % width of marginal plots
DEFAULT_VERTS   = [];           % default no vertical lines

if (nargin < 2)
	help crossf
	return
end

if (min(size(X))~=1 || length(X)<2)
	fprintf('crossf(): xdata must be a row or column vector.\n');
    return
elseif (min(size(Y))~=1 || length(Y)<2)
	fprintf('crossf(): ydata must be a row or column vector.\n');
    return
elseif (length(X) ~= length(Y))
	fprintf('crossf(): xdata and ydata must have same length.\n');
    return
end

if (nargin < 3)
	epoch = DEFAULT_EPOCH;
elseif (~isnumeric(epoch) || length(epoch)~=1 || epoch~=round(epoch))
	fprintf('crossf(): Value of frames must be an integer.\n');
    return
elseif (epoch <= 0)
	fprintf('crossf(): Value of frames must be positive.\n');
    return
elseif (rem(length(X),epoch) ~= 0)
	fprintf('crossf(): Length of data vectors must be divisible by frames.\n');
    return
end

if (nargin < 4)
	timelim = DEFAULT_TIMELIM;
elseif (~isnumeric(timelim) || sum(size(timelim))~=3)
	error('crossf(): Value of tlimits must be a vector containing two numbers.');
elseif (timelim(1) >= timelim(2))
	error('crossf(): tlimits interval must be [min,max].');
end

if (nargin < 5)
	ftitle = DEFAULT_TITLE;
elseif (~ischar(ftitle))
	error('crossf(): Plot title argument must be a quoted string.');
end

if (nargin < 6)
	Fs = DEFAULT_FS;
elseif (~isnumeric(Fs) || length(Fs)~=1)
	error('crossf(): Value of srate must be a number.');
elseif (Fs <= 0)
	error('crossf(): Value of srate must be positive.');
end

if (nargin < 7)
	varwin = DEFAULT_VARWIN;
elseif (~isnumeric(varwin) || length(varwin)~=1)
	error('crossf(): Value of cycles must be a number.');
elseif (varwin < MAX_BASELN)
	error('crossf(): Value of cycles must be either zero or positive.');
end

if (nargin < 8)
	winsize = max(pow2(nextpow2(epoch)-3),4);
elseif (~isnumeric(winsize) || length(winsize)~=1 || winsize~=round(winsize))
	error('crossf(): Value of winsize must be an integer number.');
elseif (winsize <= 0)
	error('crossf(): Value of winsize must be positive.');
elseif (varwin == 0 && pow2(nextpow2(winsize)) ~= winsize)
	error('crossf(): Value of winsize must be an integer power of two [1,2,4,8,16,...]');
elseif (winsize > epoch)
	error('crossf(): Value of winsize must be less than epoch length.');
end

if (nargin < 9)
	nwin = DEFAULT_NWIN;
elseif (~isnumeric(nwin) || length(nwin)~=1 || nwin~=round(nwin))
	error('crossf(): Value of nwin must be an integer number.');
elseif (nwin <= 0)
	error('crossf(): Value of nwin must be positive.');
end
if (nwin > epoch-winsize)
	error('crossf(): Value of nwin must be <= epoch-winsize.');
end

if (nargin < 10)
	oversmp = DEFAULT_OVERSMP;
elseif (~isnumeric(oversmp) || length(oversmp)~=1 || oversmp~=round(oversmp))
	error('crossf(): Value of oversmp must be an integer number.');
elseif (oversmp <= 0)
	error('crossf(): Value of oversmp must be positive.');
elseif (pow2(nextpow2(oversmp)) ~= oversmp)
	error('crossf(): Value of oversmp must be an integer power of two [1,2,4,8,16,...]');
end

if (nargin < 11)
	maxfreq = DEFAULT_MAXFREQ;
elseif (~isnumeric(maxfreq) || length(maxfreq)~=1)
	error('crossf(): Value of maxfreq must be a number.');
elseif (maxfreq <= 0)
	error('crossf(): Value of maxfreq must be positive.');
end

if (nargin < 12)
	alpha = DEFAULT_ALPHA;
elseif (~isnumeric(alpha) || length(alpha)~=1)
	error('crossf(): Value of alpha must be a number.');
elseif (round(NACCU*alpha) < 2 || alpha > .5)
	fprintf('crossf(): Value of alpha must be in the range (~0,0.5]');
    return 
else
    if round(NACCU*alpha)<1,
      alpha = 1/NACCU;
	  fprintf(...
        'Using alpha = %0.3f. To decrease, must raise NACCU in source code.\n',...
                       alpha);
    end
end
if (nargin < 13)
   verts = DEFAULT_VERTS;
end

if (varwin == 0) % FFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	freqs = Fs/winsize*[1:2/oversmp:winsize]/2;
	win = hanning(winsize);

	R = zeros(oversmp*winsize/2,nwin);
	RR = zeros(oversmp*winsize/2,nwin);
	Rboot = zeros(oversmp*winsize/2,NACCU);

else % wavelet DFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	freqs = Fs*varwin/winsize*[2:2/oversmp:winsize]/2;
	win = dftfilt(winsize,maxfreq/Fs,varwin,oversmp,1);

	R  = zeros(size(win,2),nwin);
	RR = zeros(size(win,2),nwin);
	Rboot = zeros(size(win,2),NACCU);
end

wintime = 500*winsize/Fs;
times = timelim(1)+wintime:(timelim(2)-timelim(1)-2*wintime)/(nwin-1):timelim(2)-wintime;

baseln = find(times < 0);
dispf = find(freqs <= maxfreq);
stp = (epoch-winsize)/(nwin-1);
trials = length(X)/epoch;

fprintf('\nComputing the Event-Related Cross-Coherence image\n');
fprintf(' based on %d trials of %d frames sampled at %g Hz.\n',...
                       trials,epoch,Fs);
fprintf('Trial timebase is %d ms before to %d ms after the stimulus\n',...
                       timelim(1),timelim(2));
fprintf('The frequency range displayed is %g-%g Hz.\n',min(dispf),maxfreq);
if varwin==0
  fprintf('The data window size is %d samples (%g ms).\n',winsize,2*wintime);
  fprintf('The FFT length is %d samples\n',winsize*oversmp);
else
  fprintf('The window size is %d cycles.\n',varwin);
  fprintf('The maximum window size is %d samples (%g ms).\n',winsize,2*wintime);
end
fprintf('The window is applied %d times\n',nwin);
fprintf(' with an average step size of %g samples (%g ms).\n',...
                       stp,1000*stp/Fs);
fprintf('Results are oversampled %d times.\n',oversmp);
if nargout>4
  fprintf('Bootstrap confidence limits will be computed based on alpha = %g\n',...
              alpha);
else
  fprintf('Bootstrap confidence limits will NOT be computed.\n'); 
end
if nargout>5
  fprintf(['Coherence angles will be imaged in ',ANGLEUNITS,' and saved.\n']);
end

fprintf('\nProcessing trial (of %d):',trials);
firstboot = 1;
Rn=zeros(1,nwin);
X = X(:)'; % make X and Y column vectors
Y = Y(:)';
for t=1:trials,
	if (rem(t,10) == 0)
		fprintf(' %d',t);
	end
    if rem(t,120) == 0
        fprintf('\n');
    end

	for j=1:nwin, % for each time window
		tmpX = X([1:winsize]+floor((j-1)*stp)+(t-1)*epoch);
		tmpY = Y([1:winsize]+floor((j-1)*stp)+(t-1)*epoch);

        if ~any(isnan(tmpX))
		  tmpX = tmpX - mean(tmpX);
		  tmpY = tmpY - mean(tmpY);

		  if varwin == 0 % use FFTs
			tmpX = win .* tmpX(:);
			tmpY = win .* tmpY(:);
			tmpX = fft(tmpX,oversmp*winsize);
			tmpY = fft(tmpY,oversmp*winsize);
			tmpX = tmpX(2:oversmp*winsize/2+1);
			tmpY = tmpY(2:oversmp*winsize/2+1);
		  else 
			tmpX = win' * tmpX(:);
			tmpY = win' * tmpY(:);
		  end

          if nargout > 4
           if firstboot==1
             tmpsX = repmat(nan,length(tmpX),nwin);
             tmpsY = repmat(nan,length(tmpY),nwin);
             firstboot = 0;
           end
           tmpsX(:,j) = tmpX;
           tmpsY(:,j) = tmpY;
          end
		
		  RR(:,j) = tmpX.*conj(tmpY) ./ (abs(tmpX).*abs(tmpY)); % complex coher.
	      R(:,j) = R(:,j) + RR(:,j);
          Rn(j) = Rn(j)+1;
        end % ~any(isnan())
	end % time window
	
	if (nargout > 4) % get NACCU bootstrap estimates for each trial
        j=1;
		while j<=NACCU
           s = ceil(rand([1 2])*nwin); % random ints [1,nwin]
           tmpX = tmpsX(:,s(1));
           tmpY = tmpsY(:,s(2));
           if ~any(isnan(tmpX)) && ~any(isnan(tmpY))
		      RR = tmpX.*conj(tmpY) ./ (abs(tmpX).*abs(tmpY)); % complex coher.
		      Rboot(:,j) = Rboot(:,j) + RR;
              j = j+1;
           end
        end
	end
	
end % t = trial

fprintf('\nNow plotting...\n');

Rangle = angle(R);
if varwin ~= 0
   Rangle = -Rangle; % make lead/lag the same for FFT and wavelet analysis
end
R = abs(R) ./ (ones(size(R,1),1)*Rn);               % coherence magnitude
Rraw = R;						% raw coherence values
mbase = mean(R(:,baseln)');     % mean baseline coherence magnitude
% R = R - repmat(mbase',[1 nwin]);% remove baseline mean

if (nargout > 4) % bootstrap
	i = round(NACCU*alpha);
	Rboot = abs(Rboot) / trials; % normalize bootstrap magnitude to [0,1]
	Rboot = sort(Rboot');  
	Rsignif = mean(Rboot(NACCU-i+1:NACCU,:)); % significance levels for Rraw
%	Rboot = [mean(Rboot(1:i,:))-mbase ; mean(Rboot(NACCU-i+1:NACCU,:))-mbase];
 	Rboot = [mean(Rboot(1:i,:)) ; mean(Rboot(NACCU-i+1:NACCU,:))];
end % NOTE: above, mean ?????

set(gcf,'DefaultAxesFontSize',AXES_FONT)
colormap(jet(256));

pos = get(gca,'position'); % plot relative to current axes
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)];
axis('off')

if nargout>5 % image coherence lag as well as coherence magnitude
  ybase = 0.5;
  subheight = 0.4;
  MARGIN = MARGIN*0.75;
else
  ybase = 0.0;
  subheight = 0.9;
end

%
% Image the coherence [% perturbations]
%
RR = R;
if (nargout > 4) % zero out (and 'green out') nonsignif. R values
	RR(find((RR > repmat(Rboot(1,:)',[1 nwin])) ...
              & (RR < repmat(Rboot(2,:)',[1 nwin])))) = 0;
end
if (nargout > 5) % zero out nonsignif. Rraw values
	Rraw(find(repmat(Rsignif',[1,size(Rraw,2)])>=Rraw))=0;
end

if COH_CAXIS_LIMIT == 0
    coh_caxis = max(max(R(dispf,:)))*[-1 1];
else
    coh_caxis = COH_CAXIS_LIMIT*[-1 1];
end

h(6) = axes('Units','Normalized',...
               'Position',[MARGIN ybase+MARGIN 0.9-MARGIN subheight].*s+q);

map=hsv(300); % install circular color map - green=0, yellow, orng, red, violet = max
              %                                         cyan, blue, violet = min
map = flipud([map(251:end,:);map(1:250,:)]);
map(151,:) = map(151,:)*0.9; % tone down the (0=) green!
colormap(map);

imagesc(times,freqs(dispf),RR(dispf,:),coh_caxis); % plot the coherence image

set(h(6),'Units','Normalized',...
               'Position',[MARGIN ybase+MARGIN 0.9-MARGIN subheight].*s+q);
hold on
plot([0 0],[0 freqs(max(dispf))],'--m','LineWidth',LINEWIDTH)
for i=1:length(verts)
  plot([verts(i) verts(i)],[0 freqs(max(dispf))],'--m','LineWidth',LINEWIDTH);
end
hold off
set(h(6),'YTickLabel',[],'YTick',[])
set(h(6),'XTickLabel',[],'XTick',[])
title('Event-Related Coherence')
h(8) = axes('Position',[.92 ybase+MARGIN .05 subheight].*s+q);
cbar(h(8),151:300,[0 coh_caxis(2)]); % use only positive colors (gyorv) 
%                                      for coherences
%
% Plot delta-mean min and max coherence at each time point on bottom of image
%
h(10) = axes('Units','Normalized','Position',[MARGIN ybase 0.9-MARGIN MARGIN].*s+q); 
Emax = max(R(dispf,:)); % mean coherence at each time point
Emin = min(R(dispf,:)); % mean coherence at each time point
plot(times,Emax,'b');
hold on
plot(times,Emin,'b');
plot([times(1) times(length(times))],[0 0],'LineWidth',0.7);
plot([0 0],[-500 500],'--m','LineWidth',LINEWIDTH);
for i=1:length(verts)
  plot([verts(i) verts(i)],[-500 500],'--m','LineWidth',LINEWIDTH);
end
axis([min(times) max(times) 0 max(Emax)*1.2])
tick = get(h(10),'YTick');
set(h(10),'YTick',[tick(1) ; tick(length(tick))])
set(h(10),'YAxisLocation','right')
midpos = get(h(10),'Position');
if nargout<6
 xlabel('Time (ms)')
end
ylabel('coh.')

%
% Plot mean baseline coherence at each freq on left side of image
%

h(11) = axes('Units','Normalized','Position',[0 ybase+MARGIN MARGIN subheight].*s+q);
E = mbase(dispf); % baseline mean coherence at each frequency
if (nargout > 4) % plot bootstrap significance limits (base mean +/-)
	plot(freqs(dispf),E,'m','LineWidth',LINEWIDTH); % plot mbase
    hold on
	% plot(freqs(dispf),Rboot(:,dispf)+[E;E],'g','LineWidth',LINEWIDTH);
	plot(freqs(dispf),Rboot([1 2],dispf),'g','LineWidth',LINEWIDTH);
	plot(freqs(dispf),Rsignif(dispf),'k:','LineWidth',LINEWIDTH);
	axis([freqs(1) freqs(max(dispf)) 0 max([E Rsignif])*1.2]);
else             % plot marginal mean coherence only
	plot(freqs(dispf),E,'LineWidth',LINEWIDTH);
	% axis([freqs(1) freqs(max(dispf)) min(E)-max(E)/3 max(E)+max(E)/3]);
	if ~isnan(max(E))
	   axis([freqs(1) freqs(max(dispf)) 0 max(E)*1.2]);
	end;   
end

tick = get(h(11),'YTick');
set(h(11),'YTick',[tick(1) ; tick(length(tick))])
set(h(11),'View',[90 90])
xlabel('Freq. (Hz)')
ylabel('coh.')

if (length(ftitle) > 0) % plot title
	axes('Position',pos,'Visible','Off');               
	h(12) = text(-.05,1.01,ftitle);
	set(h(12),'VerticalAlignment','bottom')
	set(h(12),'HorizontalAlignment','left')
	set(h(12),'FontSize',TITLE_FONT)
end
%
% Plot coherence time lags in bottom panel
%
if nargout>5
   h(13) = axes('Units','Normalized','Position',[MARGIN MARGIN 0.9-MARGIN subheight].*s+q);
   if strcmp(ANGLEUNITS,'ms')  % convert to ms
     Rangle = (Rangle/(2*pi)).*repmat(1000./freqs(dispf)',1,length(times)); 
     maxangle = max(max(abs(Rangle)));
   else
     Rangle = Rangle*180/pi; % convert to degrees
     maxangle = 180; % use full-cycle plotting 
   end
   Rangle(find(Rraw==0)) = 0; % set angle at non-signif coher points to 0

   imagesc(times,freqs(dispf),Rangle(dispf,:),[-maxangle maxangle]); % plot the 
   hold on                                             % coherence phase angles
   plot([0 0],[0 freqs(max(dispf))],'--m','LineWidth',LINEWIDTH); % zero-time line
   for i=1:length(verts)
     plot([verts(i) verts(i)],[0 freqs(max(dispf))],'--m','LineWidth',LINEWIDTH);
   end

   pos13 = get(h(13),'Position');
   set(h(13),'Position',[pos13(1) pos13(2) midpos(3) pos13(4)]);
   ylabel('Freq. (Hz)')
   xlabel('Time (ms)')
   h(14)=axes('Position',[.92 MARGIN .05 subheight].*s+q);
   cbar(h(14),0,[-maxangle maxangle]); % two-sided colorbar

   if (length(ftitle) > 0) % plot title
	axes('Position',pos,'Visible','Off');               
	h(13) = text(-.05,1.01,ftitle);
	set(h(13),'VerticalAlignment','bottom')
	set(h(13),'HorizontalAlignment','left')
	set(h(13),'FontSize',TITLE_FONT)
   end
 end
axcopy(gcf);
