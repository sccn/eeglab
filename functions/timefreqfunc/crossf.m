% crossf() - Returns estimates and plots event-related coherence (ERCOH) changes
%            between two input time series (x,y). A lower panel (optionally) shows 
%            the coherence phase difference between the processes. In this panel: 
%               -90 degrees (blue)   means x leads y by a quarter cycle.
%                90 degrees (orange) means y leads x by a quarter cycle.
%            Click on each subplot to view separately and zoom in/out.
%
% Function description:
%            Uses EITHER fixed-window, zero-padded FFTs (faster) OR constant-Q 
%            0-padded wavelet DFTs (better sensitivity), both Hanning-tapered. 
%            Output frequency spacing is the lowest frequency ('srate'/'winsize') 
%            divided by the 'padratio'.
%
%            If an 'alpha' value is given, then bootstrap statistics are 
%            computed (from a distribution of 200 ('naccu') surrogate baseline
%            data epochs) for the baseline epoch, and non-significant features 
%            of the output plots are zeroed (e.g., plotted in green). The baseline
%            epoch is all windows with center times < the 'baseline' value or, 
%            if 'baseboot' is 1, the whole epoch. 
% Usage: 
%        >> [coh,mcoh,timesout,freqsout,cohboot,cohangles] ...
%                       = crossf(x,y,frames,tlimits,titl,          ...
%                                    srate,cycles,winsize,timesout,...
%                                              padratio,maxfreq,alpha,verts);
%
% Required inputs:
%       x           = first single-channel data  (1,frames*nepochs)      {none}
%       y           = second single-channel data (1,frames*nepochs)      {none}
%       frames      = frames per epoch                                   {750}
%       tlimits     = [mintime maxtime] (ms) epoch time limits  {[-1000 2000]}
%       srate       = data sampling rate (Hz)                            {250}
%       cycles      = If >0 -> Number of cycles in each analysis wavelet 
%                     If==0 -> Use FFTs (constant window length 'winsize') {0}
%
%    Optional Coherence Type:
%       'type'      = ['coher'|'phasecoher'] Compute either linear coherence
%                      ('coher') or phase coherence ('phasecoher') also known
%                      as phase coupling factor' { 'phasecoher' }.
%    Optional Detrend:
%       'detret'    = ['on'|'off'], Detrend data within epochs.   {'off'}
%       'detrep'    = ['on'|'off'], Detrend data across trials    {'off'}
%
%    Optional FFT/DFT:
%       'winsize'   = If cycles==0: data subwindow length (fastest, 2^n<frames);
%                      if cycles >0: *longest* window length to use. This
%                      determines the lowest output frequency  {~frames/8}
%       'timesout'  = Number of output times (int<frames-winsize) {200}
%       'padratio'  = FFTlength/winsize (2^k)                     {2}
%                      Multiplies the number of output frequencies by
%                      dividing their spacing. When cycles==0, frequency
%                      spacing is (low_frequency/padratio).
%       'maxfreq'   = Maximum frequency (Hz) to plot (& output if cycles>0) {50}
%                      If cycles==0, all FFT frequencies are output.
%       'baseline'  = Coherence baseline end time (ms).           {0}
%       'powbase'   = Baseline spectrum to log-subtract.          {from data}
%
%    Optional Bootstrap:
%       'alpha'     = If non-0, compute Two-tailed bootstrap significance prob. 
%                      level. Show non-signif output values as green. {0}
%       'naccu'     = Number of bootstrap replications to compute {200}
%       'baseboot'  = Bootstrap baseline subtract (0=same as 'baseline';
%                      1=whole epoch)                              {0}
%       'boottype'  = ['trials'|'times'|'both'] Bootstrap type: Either shuffle
%                      trials but not windows ('trials'), windows but not trials
%                      ('times') or both ('both')                  {'times' }
%    Optional Scalp Map:
%       'topovec'   = Scalp topography (map) to plot              {[]}
%       'elocs'     = Electrode location file for scalp map       {none}
%                      File should be ascii in format of  >> topoplot example   
%
%    Optional Plot Features:
%       'plotamps'  = ['on'|'off'], Plot coherence magnitude      {'on'}
%       'plotphase' = ['on'|'off'], Plot coherence phase angle    {'on'}
%       'title'     = Optional figure title                       {none}
%       'marktimes' = Times to mark with a dotted vertical line   {none}
%       'linewidth' = Line width for marktimes traces (thick=2, thin=1) {2}
%       'cmax'      = Maximum amplitude for color scale  { use data limits }
%       'angleunit' = Phase units: 'ms' for msec or 'deg' for degrees {'deg'}
%       'pboot'     = Bootstrap power limits (e.g., from timef()) {from data}
%       'rboot'     = Bootstrap coherence limits (e.g., from timef()) {from data}
%       'axesfont'  = Axes font size                               {10}
%       'titlefont' = Title font size                              {8}
%
% Outputs: 
%       coh         = Matrix (nfreqs,timesout) of coherence magnitudes 
%       mcoh        = Vector of mean baseline coherence at each frequency
%       timesout    = Vector of output times (window centers) (ms).
%       freqsout    = Vector of frequency bin centers (Hz).
%       cohboot     = Matrix (2,nfreqs) of [lower;upper] coh signif. limits
%       cohangle    = (nfreqs,timesout) matrix of coherence angles 
%
% Note: when cycles==0, nfreqs is total number of FFT frequencies.
%
% Authors: Sigurd Enghoff, Arnaud Delorme & Scott Makeig
%          SCCN/INC/UCSD, La Jolla, 1998-2002 
%
% See also: timef()

% Copyright (C) 8/1/98 Sigurd Enghoff, Arnaud Delorme & Scott Makeig, SCCN/INC/UCSD
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
% Revision 1.2  2002/04/07 02:24:36  scott
% worked on hlpe message, changed some defaults -sm
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 11-20-98 defined g.linewidth constant -sm
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
% 03-09-02 function restructuration -ad
%  add 'key', val arguments (+ external baseboot, baseline, color axis, angleunit...)
%  add detrending (across time and trials) + 'coher' option for amplitude coherence
%  significance only if alpha is given, ploting options in 'plotamp' and 'plotphase'
% 03-16-02 timeout automatically adjusted if too high -ad 
% 04-03-02 added new options for bootstrap -ad 

function [R,mbase,times,freqs,Rboot,Rangle,Rsignif] = crossf(X, Y, frame, tlimits, Fs, varwin, varargin)

%varwin,winsize,nwin,oversmp,maxfreq,alpha,verts,caxmax)

% Commandline arg defaults:
DEFAULT_ANGLEUNITS = 'deg';     % angle plotting units - 'ms' or 'deg'
DEFAULT_EPOCH	= 750;			% Frames per epoch
DEFAULT_TIMELIM = [-1000 2000];	% Time range of epochs (ms)
DEFAULT_FS		= 250;			% Sampling frequency (Hz)
DEFAULT_NWIN	= 200;			% Number of windows = horizontal resolution
DEFAULT_VARWIN	= 0;			% Fixed window length or base on cycles.
								% =0: fix window length to nwin
								% >0: set window length equal varwin cycles
								%     bounded above by winsize, also determines
								%     the min. freq. to be computed.
DEFAULT_OVERSMP	= 2;			% Number of times to oversample = vertical resolution
DEFAULT_MAXFREQ = 50;			% Maximum frequency to display (Hz)
DEFAULT_TITLE	= 'Event-Related Coherence';			% Figure title
DEFAULT_ALPHA   = NaN;			% Default two-sided significance probability threshold
AXES_FONT       = 10;
TITLE_FONT = 14;
           
if (nargin < 2)
	help crossf
	return
end

if (min(size(X))~=1 | length(X)<2)
	fprintf('crossf(): x must be a row or column vector.\n');
    return
elseif (min(size(Y))~=1 | length(Y)<2)
	fprintf('crossf(): y must be a row or column vector.\n');
    return
elseif (length(X) ~= length(Y))
	fprintf('crossf(): x and y must have same length.\n');
    return
end

if (nargin < 3)
	frame = DEFAULT_EPOCH;
elseif (~isnumeric(frame) | length(frame)~=1 | frame~=round(frame))
	fprintf('crossf(): Value of frames must be an integer.\n');
    return
elseif (frame <= 0)
	fprintf('crossf(): Value of frames must be positive.\n');
    return
elseif (rem(length(X),frame) ~= 0)
	fprintf('crossf(): Length of data vectors must be divisible by frames.\n');
    return
end

if (nargin < 4)
	tlimits = DEFAULT_TIMELIM;
elseif (~isnumeric(tlimits) | sum(size(tlimits))~=3)
	error('crossf(): Value of tlimits must be a vector containing two numbers.');
elseif (tlimits(1) >= tlimits(2))
	error('crossf(): tlimits interval must be [min,max].');
end

if (nargin < 5)
	Fs = DEFAULT_FS;
elseif (~isnumeric(Fs) | length(Fs)~=1)
	error('crossf(): Value of srate must be a number.');
elseif (Fs <= 0)
	error('crossf(): Value of srate must be positive.');
end

if (nargin < 6)
	varwin = DEFAULT_VARWIN;
elseif (~isnumeric(varwin) | length(varwin)~=1)
	error('crossf(): Value of cycles must be a number.');
elseif (varwin < 0)
	error('crossf(): Value of cycles must be either zero or positive.');
end

% consider structure for these arguments
% --------------------------------------
if ~isempty(varargin)
    try, g = struct(varargin{:}); 
    catch, error('Argument error in the {''param'', value} sequence'); end; 
end;
g.tlimits = tlimits;
g.frame   = frame;
g.srate   = Fs;
g.cycles  = varwin;

try, g.title;      catch, g.title = DEFAULT_TITLE; end;
try, g.winsize;    catch, g.winsize = max(pow2(nextpow2(g.frame)-3),4); end;
try, g.pad;        catch, g.pad = max(pow2(nextpow2(g.winsize)),4); end;
try, g.timesout;   catch, g.timesout = DEFAULT_NWIN; end;
try, g.padratio;   catch, g.padratio = DEFAULT_OVERSMP; end;
try, g.maxfreq;    catch, g.maxfreq = DEFAULT_MAXFREQ; end;
try, g.topovec;    catch, g.topovec = []; end;
try, g.elocs;      catch, g.elocs = ''; end;
try, g.alpha;      catch, g.alpha = DEFAULT_ALPHA; end;  
try, g.marktimes;  catch, g.marktimes = []; end; % default no vertical lines
try, g.powbase;    catch, g.powbase = nan; end;
try, g.pboot;      catch, g.pboot = nan; end;
try, g.rboot;      catch, g.rboot = nan; end;
try, g.plotamp;    catch, g.plotamp = 'on'; end;
try, g.plotphase;  catch, g.plotphase  = 'on'; end;
try, g.detrep;     catch, g.detrep = 'off'; end;
try, g.detret;     catch, g.detret = 'off'; end;
try, g.baseline;   catch, g.baseline = 0; end;
try, g.baseboot;   catch, g.baseboot = 0; end;
try, g.linewidth;  catch, g.linewidth = 2; end;
try, g.naccu;      catch, g.naccu = 200; end;
try, g.angleunit;  catch, g.angleunit = DEFAULT_ANGLEUNITS; end;
try, g.cmax;       catch, g.cmax = 0; end; % 0=use data limits
try, g.type;       catch, g.type = 'phasecoher'; end; 
try, g.boottype;   catch, g.boottype = 'times'; end; 

% testing arguments consistency
% -----------------------------
if (~ischar(g.title))
	error('Title must be a string.');
end

if (~isnumeric(g.winsize) | length(g.winsize)~=1 | g.winsize~=round(g.winsize))
	error('Value of winsize must be an integer number.');
elseif (g.winsize <= 0)
	error('Value of winsize must be positive.');
elseif (g.cycles == 0 & pow2(nextpow2(g.winsize)) ~= g.winsize)
	error('Value of winsize must be an integer power of two [1,2,4,8,16,...]');
elseif (g.winsize > g.frame)
	error('Value of winsize must be less than frame length.');
end

if (~isnumeric(g.timesout) | length(g.timesout)~=1 | g.timesout~=round(g.timesout))
	error('Value of timesout must be an integer number.');
elseif (g.timesout <= 0)
	error('Value of timesout must be positive.');
end
if (g.timesout > g.frame-g.winsize)
	g.timesout = g.frame-g.winsize;
	disp(['Value of timesout must be <= frame-winsize, timeout adjusted to ' int2str(g.timesout) ]);
end

if (~isnumeric(g.padratio) | length(g.padratio)~=1 | g.padratio~=round(g.padratio))
	error('Value of padratio must be an integer.');
elseif (g.padratio <= 0)
	error('Value of padratio must be positive.');
elseif (pow2(nextpow2(g.padratio)) ~= g.padratio)
	error('Value of padratio must be an integer power of two [1,2,4,8,16,...]');
end

if (~isnumeric(g.maxfreq) | length(g.maxfreq)~=1)
	error('Value of g.maxfreq must be a number.');
elseif (g.maxfreq <= 0)
	error('Value of g.maxfreq must be positive.');
elseif (g.maxfreq > Fs/2)
	fprintf('Warning: value of g.maxfreq greater that Nyquist rate\n\n');
end

if isempty(g.topovec)
	g.topovec = [];
elseif (min(size(g.topovec))~=1)
	error('tvec must be a row or column vector.');
end

if isempty(g.elocs)
	g.elocs = '';
elseif (~ischar(g.elocs))
	error('Channel location file must be a valid text file.');
end

if (~isnumeric(g.alpha) | length(g.alpha)~=1)
	error('timef(): Value of g.alpha must be a number.\n');
elseif (round(g.naccu*g.alpha) < 2)
	fprintf('Value of g.alpha is out of the normal range [%g,0.5]\n',2/g.naccu);
    g.naccu = round(2/g.alpha);
	fprintf('  Increasing the number of bootstrap iterations to %d\n',g.naccu);
end
if g.alpha>0.5 | g.alpha<=0
    error('Value of g.alpha is out of the allowed range (0.00,0.5).');
end
if ~isnan(g.alpha)
   if g.baseboot > 0
     fprintf('Bootstrap analysis will use data in baseline (pre-0) subwindows only.\n')
   else
     fprintf('Bootstrap analysis will use data in all subwindows.\n')
   end
end
switch g.angleunit
    case { 'ms', 'deg' },;
    otherwise error('Angleunit must be either ''deg'' or ''ms''');
end;    
switch g.type
    case { 'coher', 'phasecoher' },;
    otherwise error('Type must be either ''coher'' or ''phasecoher''');
end;    
switch g.boottype
    case { 'trials', 'times', 'both' },;
    otherwise error('Boot type must be either ''trials'', ''times'' or ''both''');
end;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% winsize, nwin, oversmp, maxfreq, alpha, vert =marktimes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (g.cycles == 0) %%%%%%%%%%%%%% constant window-length FFTs %%%%%%%%%%%%%%%%
    freqs = g.srate/g.winsize*[1:2/g.padratio:g.winsize]/2;
    win = hanning(g.winsize);

	R  = zeros(g.padratio*g.winsize/2,g.timesout); % mean coherence
	RR = zeros(g.padratio*g.winsize/2,g.timesout); % (coherence)
	Rboot = zeros(g.padratio*g.winsize/2,g.naccu); % summed bootstrap coher
	switch g.type
	    case 'coher',
           cumulX = zeros(g.padratio*g.winsize/2,g.timesout);
           cumulY = zeros(g.padratio*g.winsize/2,g.timesout);
           cumulXboot = zeros(g.padratio*g.winsize/2,g.naccu);
           cumulYboot = zeros(g.padratio*g.winsize/2,g.naccu);
    end;
    switch g.boottype
        case 'trials'
	       Rboot = zeros(g.padratio*g.winsize/2, g.timesout,g.naccu); % summed bootstrap coher
           cumulXboot = zeros(g.padratio*g.winsize/2, g.timesout, g.naccu);
           cumulYboot = zeros(g.padratio*g.winsize/2, g.timesout, g.naccu);
    end;            
else % %%%%%%%%%%%%%%%%%% Constant-Q (wavelet) DFTs %%%%%%%%%%%%%%%%%%%%%%%%%%%%

   	freqs = g.srate*g.cycles/g.winsize*[2:2/g.padratio:g.winsize]/2;
	win = dftfilt(g.winsize,g.maxfreq/g.srate,g.cycles,g.padratio,.5);

	R = zeros(size(win,2),g.timesout);       % mean coherence
	RR = repmat(nan,size(win,2),g.timesout); % initialize with nans
	Rboot = zeros(size(win,2),g.naccu);  % summed bootstrap coher
	switch g.type
	    case 'coher',
           cumulX = zeros(size(win,2),g.timesout);
           cumulY = zeros(size(win,2),g.timesout);
           cumulXboot = zeros(size(win,2),g.naccu);
           cumulYboot = zeros(size(win,2),g.naccu);
    end;        
    switch g.boottype
        case 'trials'
	       Rboot = zeros(size(win,2), g.timesout,g.naccu); % summed bootstrap coher
           cumulXboot = zeros(size(win,2), g.timesout, g.naccu);
           cumulYboot = zeros(size(win,2), g.timesout, g.naccu);
    end;            
end

wintime = 500*g.winsize/g.srate;
times = [g.tlimits(1)+wintime:(g.tlimits(2)-g.tlimits(1)-2*wintime)/(g.timesout-1):g.tlimits(2)-wintime];

if g.baseboot
   baseln = 1:length(times); % use all times as baseline
else   baseln = find(times < g.baseline); % subtract means of pre-0 (centered) windows
   if isempty(baseln)
       baseln = 1:length(times); % use all times as baseline
       disp('Bootstrap baseline empty, using the whole epoch');
   end;     
end
dispf = find(freqs <= g.maxfreq);
stp = (g.frame-g.winsize)/(g.timesout-1);
trials = length(X)/g.frame;

fprintf('\nComputing the Event-Related Cross-Coherence image\n');
fprintf(' based on %d trials of %d frames sampled at %g Hz.\n',...
                       trials,g.frame,g.srate);
fprintf('Trial timebase is %d ms before to %d ms after the stimulus\n',...
                       g.tlimits(1),g.tlimits(2));
fprintf('The frequency range displayed is %g-%g Hz.\n',min(dispf),g.maxfreq);
if g.cycles==0
  fprintf('The data window size is %d samples (%g ms).\n',g.winsize,2*wintime);
  fprintf('The FFT length is %d samples\n',g.winsize*g.padratio);
else
  fprintf('The window size is %d cycles.\n',g.cycles);
  fprintf('The maximum window size is %d samples (%g ms).\n',g.winsize,2*wintime);
end
fprintf('The window is applied %d times\n',g.timesout);
fprintf(' with an average step size of %g samples (%g ms).\n',...
                       stp,1000*stp/g.srate);
fprintf('Results are oversampled %d times.\n',g.padratio);
if ~isnan(g.alpha)
  fprintf('Bootstrap confidence limits will be computed based on alpha = %g\n',...
              g.alpha);
else
  fprintf('Bootstrap confidence limits will NOT be computed.\n'); 
end
switch g.plotphase
    case 'on', fprintf(['Coherence angles will be imaged in ',g.angleunit,'\n']);
end;

fprintf('\nProcessing trial (of %d):',trials);

% detrend over epochs (trials) if requested
% -----------------------------------------
switch g.detrep
    case 'on'
        X = reshape(X, g.frame, length(X)/g.frame);
        X = X - mean(X,2)*ones(1, length(X(:))/g.frame);
        Y = reshape(Y, g.frame, length(Y)/g.frame);
        Y = Y - mean(Y,2)*ones(1, length(Y(:))/g.frame);
end;        

firstboot = 1;
Rn=zeros(1,g.timesout);
X = X(:)'; % make X and Y column vectors
Y = Y(:)';
for t=1:trials,
	if (rem(t,10) == 0)
		fprintf(' %d',t);
	end
    if rem(t,120) == 0
        fprintf('\n');
    end

	for j=1:g.timesout, % for each time window
		tmpX = X([1:g.winsize]+floor((j-1)*stp)+(t-1)*g.frame);
		tmpY = Y([1:g.winsize]+floor((j-1)*stp)+(t-1)*g.frame);

        if ~any(isnan(tmpX))
		  tmpX = tmpX - mean(tmpX);
		  tmpY = tmpY - mean(tmpY);
          switch g.detret, case 'on', 
              tmpX = detrend(tmpX); 
              tmpY = detrend(tmpY); 
          end;

		  if g.cycles == 0 % use FFTs
			tmpX = win .* tmpX(:);
			tmpY = win .* tmpY(:);
			tmpX = fft(tmpX,g.padratio*g.winsize);
			tmpY = fft(tmpY,g.padratio*g.winsize);
			tmpX = tmpX(2:g.padratio*g.winsize/2+1);
			tmpY = tmpY(2:g.padratio*g.winsize/2+1);
		  else 
			tmpX = win' * tmpX(:);
			tmpY = win' * tmpY(:);
		  end

          if ~isnan(g.alpha)
           if firstboot==1
             tmpsX = repmat(nan,length(tmpX),g.timesout);
             tmpsY = repmat(nan,length(tmpY),g.timesout);
             firstboot = 0;
           end
           tmpsX(:,j) = tmpX;
           tmpsY(:,j) = tmpY;
          end

		  switch g.type
		      case 'coher',
		          R(:,j)      = R(:,j) + tmpX.*conj(tmpY); % complex coher.
                  cumulX(:,j) = cumulX(:,j)+abs(tmpX);
                  cumulY(:,j) = cumulY(:,j)+abs(tmpY);
		      case 'phasecoher',
		          R(:,j) = R(:,j) + tmpX.*conj(tmpY) ./ (abs(tmpX).*abs(tmpY)); % complex coher.
          end;
          Rn(j) = Rn(j)+1;
        end % ~any(isnan())
	end % time window
	
	if ~isnan(g.alpha) 
	   if strcmp(g.boottype, 'times') % get g.naccu bootstrap estimates for each trial
        j=1;
		while j<=g.naccu
           s = ceil(rand([1 2])*g.timesout); % random ints [1,g.timesout]
           tmpX = tmpsX(:,s(1));
           tmpY = tmpsY(:,s(2));
           if ~any(isnan(tmpX)) & ~any(isnan(tmpY))
		      switch g.type
		          case 'coher',
		              Rboot(:,j) = Rboot(:,j) + tmpX.*conj(tmpY); % complex coher.
                      cumulXboot(:,j) = cumulXboot(:,j)+abs(tmpX);
                      cumulYboot(:,j) = cumulYboot(:,j)+abs(tmpY);
		          case 'phasecoher',
		              Rboot(:,j) = Rboot(:,j) + tmpX.*conj(tmpY) ./ (abs(tmpX).*abs(tmpY)); % complex coher.
              end;
              j = j+1;
           end
        end
      else
        alltmpsX{t} = tmpsX;
        alltmpsY{t} = tmpsX;
      end;
	end
	
end % t = trial

% handle specific bootstrap types
% -------------------------------
if ~isnan(g.alpha) & ~strcmp(g.boottype, 'times')
    fprintf('\nProcessing bootstrap (of %d):',trials);
    for allt=1:trials
		if (rem(allt,10) == 0)
			fprintf(' %d',allt);
		end
	    if rem(allt,120) == 0
	        fprintf('\n');
	    end
	    j=1;
	    while j<=g.naccu
	        switch g.boottype
	            case 'trials',
	                t = ceil(rand([1 2])*trials); % random ints [1,g.timesout]
	                tmpsX = alltmpsX{t(1)};
	                tmpsY = alltmpsY{t(2)};
	                if ~any(isnan(tmpsX(:))) & ~any(isnan(tmpsY(:)))
	 	               switch g.type
				          case 'coher',
				              Rboot(:,:,j) = Rboot(:,:,j) + tmpsX.*conj(tmpsY); % complex coher.
		                      cumulXboot(:,:,j) = cumulXboot(:,:,j)+abs(tmpsX);
		                      cumulYboot(:,:,j) = cumulYboot(:,:,j)+abs(tmpsY);
				          case 'phasecoher',
				              Rboot(:,:,j) = Rboot(:,:,j) + tmpsX.*conj(tmpsY) ./ (abs(tmpsX).*abs(tmpsY)); % complex coher.
		               end;
			           j = j+1;
			        end
	            case 'both',
	                t = ceil(rand([1 2])*trials); % random ints [1,g.timesout]
	                tmpsX = alltmpsX{t(1)};
	                tmpsY = alltmpsY{t(2)};
	                tmpX = tmpsX(:,s(1));
	                tmpY = tmpsY(:,s(2));
			        if ~any(isnan(tmpX)) & ~any(isnan(tmpY))
				        switch g.type
					          case 'coher',
					              Rboot(:,j) = Rboot(:,j) + tmpX.*conj(tmpY); % complex coher.
			                      cumulXboot(:,j) = cumulXboot(:,j)+abs(tmpX);
			                      cumulYboot(:,j) = cumulYboot(:,j)+abs(tmpY);
					          case 'phasecoher',
					              Rboot(:,j) = Rboot(:,j) + tmpX.*conj(tmpY) ./ (abs(tmpX).*abs(tmpY)); % complex coher.
			            end;
			            j = j+1;
			        end
	        end;            
	    end
	end;    
end;

% if coherence, perform the division
% ----------------------------------
switch g.type
   case 'coher',
        R = R ./ ( cumulX .* cumulY );
	    if ~isnan(g.alpha)
            Rboot = Rboot ./ ( cumulXboot .* cumulYboot );
	    end;   
end;        

switch lower(g.plotphase)
   case 'on',  
       switch lower(g.plotamp), 
          case 'on', ordinate1 = 0.67; ordinate2 = 0.1; height = 0.33; g.plot = 1;
          case 'off', ordinate2 = 0.1; height = 0.9; g.plot = 1;
       end;     
   case 'off', ordinate1 = 0.1; height = 0.9; 
       switch lower(g.plotamp), 
          case 'on', ordinate1 = 0.1; height = 0.9;  g.plot = 1;
          case 'off', g.plot = 0;
       end;     
end;    

if g.plot
    fprintf('\nNow plotting...\n');
end;

Rangle = angle(R);
if g.cycles ~= 0
   Rangle = -Rangle; % make lead/lag the same for FFT and wavelet analysis
end
R = abs(R) ./ (ones(size(R,1),1)*Rn);               % coherence magnitude
Rraw = R;						% raw coherence values
mbase = mean(R(:,baseln)');     % mean baseline coherence magnitude
% R = R - repmat(mbase',[1 g.timesout]);% remove baseline mean

if ~isnan(g.alpha) % if bootstrap analysis included . . .
    switch g.boottype
	    case 'trials',
			i = round(g.naccu*g.alpha);
			Rboot = abs(Rboot) / trials; % normalize bootstrap magnitude to [0,1]
			Rboot = sort(Rboot,3);  
			Rsignif = mean(Rboot(:,:,g.naccu-i+1:g.naccu),3); % significance levels for Rraw
		%	Rboot = [mean(Rboot(1:i,:))-mbase ; mean(Rboot(g.naccu-i+1:g.naccu,:))-mbase];
		 	Rbootplus = mean(Rboot(:,:,1:i),3);
		 	Rbootminus = mean(Rboot(:,:,g.naccu-i+1:g.naccu),3);
        otherwise
			i = round(g.naccu*g.alpha);
			Rboot = abs(Rboot) / trials; % normalize bootstrap magnitude to [0,1]
			Rboot = sort(Rboot');  
			Rsignif = mean(Rboot(g.naccu-i+1:g.naccu,:)); % significance levels for Rraw
		%	Rboot = [mean(Rboot(1:i,:))-mbase ; mean(Rboot(g.naccu-i+1:g.naccu,:))-mbase];
		 	Rboot = [mean(Rboot(1:i,:)) ; mean(Rboot(g.naccu-i+1:g.naccu,:))];
	end;
end % NOTE: above, mean ?????

set(gcf,'DefaultAxesFontSize',AXES_FONT)
colormap(jet(256));

pos = get(gca,'position'); % plot relative to current axes
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)];
axis('off')

switch lower(g.plotamp)
 case 'on' 
    %
    % Image the coherence [% perturbations] 
    %
	RR = R;
	if ~isnan(g.alpha) % zero out (and 'green out') nonsignif. R values
        switch g.boottype
	       case 'trials',
	          size(RR)
	          size(Rboot)
		      RR(find((RR > Rbootplus) & (RR < Rbootminus))) = 0;
		      Rraw(find(Rsignif >= Rraw))=0;
		      Rboottime = [mean(Rbootplus(dispf,:),1); mean(Rbootminus(dispf,:),1)];
		      Rsigniftime = mean(Rsignif(dispf,:),1);
		      Rboot = [mean(Rbootplus,2) mean(Rbootminus,2)]';
		      Rsignif = mean(Rsignif,2)';
		   otherwise
		      RR(find((RR > repmat(Rboot(1,:)',[1 g.timesout])) ...
	              & (RR < repmat(Rboot(2,:)',[1 g.timesout])))) = 0;
		      Rraw(find(repmat(Rsignif',[1,size(Rraw,2)])>=Rraw))=0;
		   end;   
	end

	if g.cmax == 0
	    coh_caxis = max(max(R(dispf,:)))*[-1 1];
	else
	    coh_caxis = g.cmax*[-1 1];
	end

	h(6) = axes('Units','Normalized', 'Position',[.1 ordinate1 .8 height].*s+q);

	map=hsv(300); % install circular color map - green=0, yellow, orng, red, violet = max
	              %                                         cyan, blue, violet = min
	map = flipud([map(251:end,:);map(1:250,:)]);
	map(151,:) = map(151,:)*0.9; % tone down the (0=) green!
	colormap(map);

	imagesc(times,freqs(dispf),RR(dispf,:),coh_caxis); % plot the coherence image

	hold on
	plot([0 0],[0 freqs(max(dispf))],'--m','LineWidth',g.linewidth)
	for i=1:length(g.marktimes)
	  plot([g.marktimes(i) g.marktimes(i)],[0 freqs(max(dispf))],'--m','LineWidth',g.linewidth);
	end;
	hold off
	set(h(6),'YTickLabel',[],'YTick',[])
	set(h(6),'XTickLabel',[],'XTick',[])
	%title('Event-Related Coherence')

	h(8) = axes('Position',[.95 ordinate1 .05 height].*s+q);
	cbar(h(8),151:300,[0 coh_caxis(2)]); % use only positive colors (gyorv) 

	%
	% Plot delta-mean min and max coherence at each time point on bottom of image
	%
	h(10) = axes('Units','Normalized','Position',[.1 ordinate1-0.1 .8 .1].*s+q); % plot marginal means below
	Emax = max(R(dispf,:)); % mean coherence at each time point
	Emin = min(R(dispf,:)); % mean coherence at each time point
	if ~isnan(g.alpha) & strcmp(g.boottype, 'trials') % plot bootstrap significance limits (base mean +/-)
	    plot(times,Rboottime([1 2],:),'g','LineWidth',g.linewidth); hold on;
	    plot(times,Rsigniftime,'k:','LineWidth',g.linewidth);
		plot(times,Emax,'b');
		plot(times,Emin,'b');
		plot([times(1) times(length(times))],[0 0],'LineWidth',0.7);
		plot([0 0],[-500 500],'--m','LineWidth',g.linewidth);
		for i=1:length(g.marktimes)
		  plot([g.marktimes(i) g.marktimes(i)],[-500 500],'--m','LineWidth',g.linewidth);
		end;
		axis([min(times) max(times) 0 max([Emax(:)' Rsignif(:)'])*1.2])
    else
		plot(times,Emax,'b');
		hold on
		plot(times,Emin,'b');
		plot([times(1) times(length(times))],[0 0],'LineWidth',0.7);
		plot([0 0],[-500 500],'--m','LineWidth',g.linewidth);
		for i=1:length(g.marktimes)
		  plot([g.marktimes(i) g.marktimes(i)],[-500 500],'--m','LineWidth',g.linewidth);
		end;
		axis([min(times) max(times) 0 max(Emax)*1.2])
    end;
	tick = get(h(10),'YTick');
	set(h(10),'YTick',[tick(1) ; tick(length(tick))])
	set(h(10),'YAxisLocation','right')
    xlabel('Time (ms)')
	ylabel('coh.')

	%
	% Plot mean baseline coherence at each freq on left side of image
	%

	h(11) = axes('Units','Normalized','Position',[0 ordinate1 .1 height].*s+q); % plot mean spectrum
	E = mbase(dispf); % baseline mean coherence at each frequency
	if ~isnan(g.alpha) % plot bootstrap significance limits (base mean +/-)
		plot(freqs(dispf),E,'m','LineWidth',g.linewidth); % plot mbase
	    hold on
		% plot(freqs(dispf),Rboot(:,dispf)+[E;E],'g','LineWidth',g.linewidth);
		plot(freqs(dispf),Rboot([1 2],dispf),'g','LineWidth',g.linewidth);
		plot(freqs(dispf),Rsignif(dispf),'k:','LineWidth',g.linewidth);
		axis([freqs(1) freqs(max(dispf)) 0 max([E Rsignif])*1.2]);
	else             % plot marginal mean coherence only
		plot(freqs(dispf),E,'LineWidth',g.linewidth);
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
end;

switch lower(g.plotphase)
  case 'on'
   %
   % Plot coherence phase lags in bottom panel
   %
   h(13) = axes('Units','Normalized','Position',[.1 ordinate2 .8 height].*s+q);
   if strcmp(g.angleunit,'ms')  % convert to ms
     Rangle = (Rangle/(2*pi)).*repmat(1000./freqs(dispf)',1,length(times)); 
     maxangle = max(max(abs(Rangle)));
   else
     Rangle = Rangle*180/pi; % convert to degrees
     maxangle = 180; % use full-cycle plotting 
   end
   Rangle(find(Rraw==0)) = 0; % set angle at non-signif coher points to 0

   imagesc(times,freqs(dispf),Rangle(dispf,:),[-maxangle maxangle]); % plot the 
   hold on                                             % coherence phase angles
   plot([0 0],[0 freqs(max(dispf))],'--m','LineWidth',g.linewidth); % zero-time line
   for i=1:length(g.marktimes)
     plot([g.marktimes(i) g.marktimes(i)],[0 freqs(max(dispf))],'--m','LineWidth',g.linewidth);
   end;

   ylabel('Freq. (Hz)')
   xlabel('Time (ms)')

   h(14)=axes('Position',[.95 ordinate2 .05 height].*s+q);
   cbar(h(14),0,[-maxangle maxangle]); % two-sided colorbar
end

if g.plot
   if (length(g.title) > 0) % plot title
	   axes('Position',pos,'Visible','Off');               
	   h(13) = text(-.05,1.01,g.title);
	   set(h(13),'VerticalAlignment','bottom')
	   set(h(13),'HorizontalAlignment','left')
	   set(h(13),'FontSize',TITLE_FONT)
   end
   axcopy(gcf);
end;
