% crossf() - Returns estimates and plots event-related coherence (ERCOH) 
%        between two input data time series (X,Y). A lower panel (optionally) 
%        shows the coherence phase difference between the processes. 
%        In this panel, output by   > crossf(X,Y,...);
%            90 degrees (orange) means X leads Y by a quarter cycle.
%           -90 degrees (blue)   means Y leads X by a quarter cycle.
%        Coherence phase units may be radians, degrees, or msec.
%        Click on any subplot to view separately and zoom in/out.
%
% Function description:
%        Uses EITHER fixed-window, zero-padded FFTs (fastest) OR constant-Q 
%        0-padded wavelet DFTs (more even sensitivity across frequencies), 
%        both Hanning-tapered.  Output frequency spacing is the lowest 
%        frequency ('srate'/'winsize') divided by the 'padratio'.
%
%        If an 'alpha' value is given, then bootstrap statistics are 
%        computed (from a distribution of 'naccu' {200} surrogate baseline
%        data epochs) for the baseline epoch, and non-significant features 
%        of the output plots are zeroed (and shown in green). The baseline
%        epoch is all windows with center latencies < the given 'baseline' 
%        value, or if 'baseboot' is 1, the whole epoch. 
% Usage: 
%        >> [coh,mcoh,timesout,freqsout,cohboot,cohangles] ...
%                       = crossf(X,Y,frames,tlimits,srate,cycles, ...
%                                        'key1', 'val1', 'key2', val2' ...);
% Required inputs:
%       X       = first single-channel data set (1,frames*nepochs)      
%       Y       = second single-channel data set (1,frames*nepochs)     
%       frames  = frames per epoch                                 {default: 750}
%       tlimits = [mintime maxtime] (ms) epoch latency limits {def: [-1000 2000]}
%       srate   = data sampling rate (Hz)                          {default: 250}
%       cycles  = 0  -> Use FFTs (with constant window length) 
%               = >0 -> Number of cycles in each analysis wavelet 
%               = [cycles expfactor] -> if 0 < expfactor < 1,  the number 
%                 of wavelet cycles expands with frequency from cycles
%                 If expfactor = 1, no expansion; if = 0, constant
%                 window length (as in FFT)                          {default: 0}
% Optional Coherence Type:
%       'type'  = ['coher'|'phasecoher'] Compute either linear coherence
%                 ('coher') or phase coherence ('phasecoher') also known
%                 as phase coupling factor' {default: 'phasecoher'}.
%       'subitc' = ['on'|'off'] subtract stimulus locked Inter-Trial Coherence 
%                 from X and Y. This computes the  'intrinsic' coherence
%                 X and Y not arising from common synchronization to 
%                 experimental events. See notes. {default: 'off'}
%       'shuffle' = integer indicating the number of estimates to compute
%                 bootstrap coherence based on shuffled trials. This estimates
%                 the coherence arising only from time locking of X and Y
%                 to experimental events (opposite of 'subitc')      {default: 0}
% Optional Detrend:
%       'detret' = ['on'|'off'], Linearly detrend data within epochs {def: 'off'}
%       'detrep' = ['on'|'off'], Linearly detrend data across trials {def: 'off'}
%
% Optional FFT/DFT:
%       'winsize'  = If cycles==0: data subwindow length (fastest, 2^n<frames);
%                    if cycles >0: *longest* window length to use. This
%                    determines the lowest output frequency  {default: ~frames/8}
%       'timesout' = Number of output latencies (int<frames-winsize)   {def: 200}
%       'padratio' = FFTlength/winsize (2^k)                         {default: 2}
%                    Multiplies the number of output frequencies by
%                    dividing their spacing. When cycles==0, frequency
%                    spacing is (low_frequency/padratio).
%       'maxfreq'  = Maximum frequency (Hz) to plot (& output if cycles>0) 
%                    If cycles==0, all FFT frequencies are output  {default: 50}
%       'baseline' = Coherence baseline end latency (ms). NaN -> No baseline  
%                      {default:NaN}
%       'powbase'  = Baseline spectrum to log-subtract      {default: from data}
%
% Optional Bootstrap:
%       'alpha'    = If non-0, compute two-tailed bootstrap significance prob.
%                    level. Show non-signif output values as green. {def: 0}
%       'naccu'    = Number of bootstrap replications to compute {def: 200}
%       'boottype' = ['times'|'timestrials'] Bootstrap type: Either shuffle
%                    windows ('times') or windows and trials ('timestrials')
%                    Option 'timestrials' requires more memory  {default: 'times'}
%       'memory'   = ['low'|'high'] 'low' -> decrease memory use {default: 'high'}
%       'baseboot' = Extent of bootstrap shuffling (0=to 'baseline'; 1=whole epoch) 
%                    If no baseline is given (NaN), extent of bootstrap shuffling 
%                    is the whole epoch                         {default: 0}
%       'rboot'    = Input bootstrap coherence limits (e.g., from crossf()) 
%                    The bootstrap type should be identical to that used
%                    to obtain the input limits. {default: compute from data}
% Optional Scalp Map:
%       'topovec'  = (2,nchans) matrix, plot scalp maps to plot {default: []}
%                    ELSE (c1,c2), plot two cartoons showing channel locations.
%       'elocs'    = Electrode location structure or file for scalp map  
%                    {default: none}
%       'chaninfo' = Electrode location additional information (nose position...)
%                    {default: none}
%
% Optional Plot and Compute Features:
%       'compute'   = ['matlab'|'c'] Use C subroutines to speed up the
%                     computation (currently unimplemented) {def: 'matlab'}
%       'savecoher' - [0|1] 1 --> Accumulate the individual trial coherence 
%                     vectors; output them as cohangles {default: 0 = off}
%       'plotamp'   = ['on'|'off'], Plot coherence magnitude    {def: 'on'}
%       'maxamp'    = [real] Set the maximum for the amp. scale {def: auto}
%       'plotphase' = ['on'|'off'], Plot coherence phase angle  {def: 'on'}
%       'angleunit' = Phase units: 'ms' -> msec, 'deg' -> degrees,
%                     or 'rad' -> radians                  {default: 'deg'}
%       'title'     = Optional figure title                {default:  none}
%       'vert'      = Latencies to mark with a dotted vertical line 
%                                                           {default: none}
%       'linewidth' = Line width for marktimes traces (thick=2, thin=1) 
%                                                              {default: 2}
%       'cmax'      = Maximum amplitude for color scale  {def: data limits}
%       'axesfont'  = Axes font size                          {default: 10}
%       'titlefont' = Title font size                          {default: 8}
%
% Outputs: 
%       coh         = Matrix (nfreqs,timesout) of coherence magnitudes 
%       mcoh        = Vector of mean baseline coherence at each frequency
%       timesout    = Vector of output latencies (window centers) (ms).
%       freqsout    = Vector of frequency bin centers (Hz).
%       cohboot     = Matrix (nfreqs,2) of [lower;upper] coher signif. limits
%                     if 'boottype' is 'trials',  (nfreqs,timesout, 2)
%       cohangle    = (nfreqs,timesout) matrix of coherence angles (in radians)
%       cohangles   = (nfreqs,timesout,trials) matrix of single-trial coherence 
%                      angles (in radians), saved and output only if 'savecoher',1
%
% Plot description:
%   Assuming both 'plotamp' and 'plotphase' options are 'on' (=default), the upper panel
%   presents the magnitude of either phase coherence or linear coherence, depending on 
%   the 'type' parameter (above). The lower panel presents the coherence phase difference 
%   (in degrees). Click on any plot to pop up a new window (using 'axcopy()').
%   -- The upper left marginal panel shows mean coherence during the baseline period
%      (blue), and when significance is set, the significance threshold (dotted black-green).
%   -- The horizontal panel under the coherence magnitude image indicates the maximum 
%      (green) and minimum (blue) coherence values across all frequencies. When significance 
%      is set (using option 'trials' for 'boottype'), an additional curve indicates the 
%      significance threshold (dotted black-green).
%
% Notes: 1) When cycles==0, nfreqs is total number of FFT frequencies.
%        2) As noted above: 'blue' coherence angle -> X leads Y; 'red' -> Y leads X
%        3) The 'boottype' should be ideally 'timesframes', but this creates high 
%           memory demands, so the 'times' method must be used in many cases.
%        4) If 'boottype' is 'trials', the average of the complex bootstrap
%           is subtracted from the coherence to compensate for phase differences 
%           (the average is also subtracted from the bootstrap distribution). 
%           For other bootstraps, this is not necessary since the phase is random.
%        5) If baseline is non-NaN, the baseline is subtracted from
%           the complex coherence. On the left hand side of the coherence
%           amplitude image, the baseline is displayed as a magenta line
%           (if no baseline is selected, this curve represents the average
%           coherence at every given frequency).
%        6) If a out-of-memory error occurs, set the 'memory' option to 'low'
%           (Makes computation time slower; Only the 'times' bootstrap method 
%           can be used in this mode).
%
% Authors: Arnaud Delorme, Sigurd Enghoff & Scott Makeig
%          CNL/Salk Institute 1998-2001; SCCN/INC/UCSD, La Jolla, 2002-
%
% See also: timef()

% Copyright (C) 8/1/98  Arnaud Delorme, Sigurd Enghoff & Scott Makeig, SCCN/INC/UCSD
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

% 11-20-98 defined g.linewidth constant -sm
% 04-01-99 made number of frequencies consistent -se
% 06-29-99 fixed constant-Q freq indexing -se
% 08-13-99 added cohangle plotting -sm
% 08-20-99 made bootstrap more efficient -sm
% 08-24-99 allow nan values introduced by possible eventlock() preproc. -sm
% 03-05-2007 eventlock.m deprecated to eegalign.m. -tf
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

% Note: 3 "objects" (Tf, Coher and Boot) are handled by specific functions under Matlab
%    (Tf) function Tf = tfinit(...) - create object Time Frequency (Tf) associated with some data
%    (Tf) function [Tf, itcvals] = tfitc(...) - compute itc for the selected data
%    (Tf) function [Tf, itcvals] = tfitcpost(Tf, trials) - itc normlisation 
%    (Tf) function [Tf, tmpX] = tfcomp(Tf, trials, times) - compute time freq. decomposition
%    (Coher) function Coher = coherinit(...) - initialize coherence object
%    (Coher) function Coher = cohercomp(Coher, tmpX, tmpY, trial, time) - compute coherence
%    (Coher) function Coher = cohercomppost(Coher, trials) - coherence normalization
%    (Boot) function Boot = bootinit(...) - intialize bootstrap object
%    (Boot) function Boot = bootcomp(...) - compute bootstrap
%    (Boot) function [Boot, Rbootout] = bootcomppost(...) - bootstrap normalization
% and by real objects under C++ (C++ code, incomplete)

function [R,mbase,times,freqs,Rbootout,Rangle, trialcoher, Tfx, Tfy] = crossf(X, Y, frame, tlimits, Fs, varwin, varargin)

%varwin,winsize,nwin,oversmp,maxfreq,alpha,verts,caxmax)

% ------------------------
% Commandline arg defaults:
% ------------------------
DEFAULT_ANGLEUNIT = 'deg'; % angle plotting units - 'rad', 'ms', or 'deg'
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

if (nargin < 2)
   help crossf
   return
end

if ~iscell(X)
	if (min(size(X))~=1 || length(X)<2)
		fprintf('crossf(): X must be a row or column vector.\n');
		return
	elseif (min(size(Y))~=1 || length(Y)<2)
		fprintf('crossf(): Y must be a row or column vector.\n');
		return
	elseif (length(X) ~= length(Y))
		fprintf('crossf(): X and Y must have same length.\n');
		return
	end
end

if (nargin < 3)
   frame = DEFAULT_EPOCH;
elseif (~isnumeric(frame) || length(frame)~=1 || frame~=round(frame))
   fprintf('crossf(): Value of frames must be an integer.\n');
   return
elseif (frame <= 0)
   fprintf('crossf(): Value of frames must be positive.\n');
   return
elseif ~iscell(X) && (rem(length(X),frame) ~= 0)
   fprintf('crossf(): Length of data vectors must be divisible by frames.\n');
   return
end

if (nargin < 4)
   tlimits = DEFAULT_TIMELIM;
elseif (~isnumeric(tlimits) || sum(size(tlimits))~=3)
   error('crossf(): Value of tlimits must be a vector containing two numbers.');
elseif (tlimits(1) >= tlimits(2))
   error('crossf(): tlimits interval must be [min,max].');
end

if (nargin < 5)
   Fs = DEFAULT_FS;
elseif (~isnumeric(Fs) || length(Fs)~=1)
   error('crossf(): Value of srate must be a number.');
elseif (Fs <= 0)
   error('crossf(): Value of srate must be positive.');
end

if (nargin < 6)
   varwin = DEFAULT_VARWIN;
elseif (~isnumeric(varwin) || length(varwin)>2)
   error('crossf(): Value of cycles must be a number or a (1,2) vector.');
elseif (varwin < 0)
   error('crossf(): Value of cycles must be either zero or positive.');
end

% consider structure for these arguments
% --------------------------------------
vararginori = varargin;
for index=1:length(varargin)
	if iscell(varargin{index}), varargin{index} = { varargin{index} }; end
end
if ~isempty(varargin)
   try, g = struct(varargin{:}); 
   catch, error('Argument error in the {''param'', value} sequence'); end; 
else 
	g = [];
end

try, g.shuffle;    catch, g.shuffle = 0; end
try, g.title;      catch, g.title = DEFAULT_TITLE; end
try, g.winsize;    catch, g.winsize = max(pow2(nextpow2(frame)-3),4); end
try, g.pad;        catch, g.pad = max(pow2(nextpow2(g.winsize)),4); end
try, g.timesout;   catch, g.timesout = DEFAULT_NWIN; end
try, g.padratio;   catch, g.padratio = DEFAULT_OVERSMP; end
try, g.maxfreq;    catch, g.maxfreq = DEFAULT_MAXFREQ; end
try, g.topovec;    catch, g.topovec = []; end
try, g.elocs;      catch, g.elocs = ''; end
try, g.alpha;      catch, g.alpha = DEFAULT_ALPHA; end;  
try, g.marktimes;  catch, g.marktimes = []; end; % default no vertical lines
try, g.marktimes = g.vert;       catch, g.vert = []; end; % default no vertical lines
try, g.powbase;    catch, g.powbase = nan; end
try, g.rboot;      catch, g.rboot = nan; end
try, g.plotamp;    catch, g.plotamp = 'on'; end
try, g.plotphase;  catch, g.plotphase  = 'on'; end
try, g.plotbootsub;  catch, g.plotbootsub  = 'on'; end
try, g.detrep;     catch, g.detrep = 'off'; end
try, g.detret;     catch, g.detret = 'off'; end
try, g.baseline;   catch, g.baseline = NaN; end
try, g.baseboot;   catch, g.baseboot = 0; end
try, g.linewidth;  catch, g.linewidth = 2; end
try, g.naccu;      catch, g.naccu = 200; end
try, g.angleunit;  catch, g.angleunit = DEFAULT_ANGLEUNIT; end
try, g.cmax;       catch, g.cmax = 0; end; % 0=use data limits
try, g.type;       catch, g.type = 'phasecoher'; end; 
try, g.boottype;   catch, g.boottype = 'times'; end; 
try, g.subitc;     catch, g.subitc = 'off'; end
try, g.memory;     catch, g.memory = 'high'; end
try, g.compute;    catch, g.compute = 'matlab'; end
try, g.maxamp;     catch, g.maxamp = []; end
try, g.savecoher;  catch, g.savecoher = 0; end
try, g.noinput;    catch, g.noinput = 'no'; end
try, g.chaninfo;   catch, g.chaninfo = []; end

allfields = fieldnames(g);
for index = 1:length(allfields)
	switch allfields{index}
	 case { 'shuffle' 'title' 'winsize' 'pad' 'timesout' 'padratio' 'maxfreq' 'topovec' 'elocs' 'alpha' ...
		  'marktimes' 'vert' 'powbase' 'rboot' 'plotamp' 'plotphase' 'plotbootsub' 'detrep' 'detret' ...
		  'baseline' 'baseboot' 'linewidth' 'naccu' 'angleunit' 'cmax' 'type' 'boottype' 'subitc' ...
		  'memory' 'compute' 'maxamp' 'savecoher' 'noinput' 'chaninfo' };
	  case {'plotersp' 'plotitc' }, disp(['crossf warning: timef option ''' allfields{index} ''' ignored']);
	 otherwise disp(['crossf error: unrecognized option ''' allfields{index} '''']); beep; return;
	end
end

g.tlimits = tlimits;
g.frame   = frame;
g.srate   = Fs;
g.cycles  = varwin(1);
if length(varwin)>1
	g.cyclesfact = varwin(2);
else 
	g.cyclesfact = 1;
end
g.type       = lower(g.type);
g.boottype   = lower(g.boottype);
g.detrep     = lower(g.detrep);
g.detret     = lower(g.detret);
g.plotphase  = lower(g.plotphase);
g.plotbootsub = lower(g.plotbootsub);
g.subitc     = lower(g.subitc);
g.plotamp    = lower(g.plotamp);
g.shuffle    = lower(g.shuffle);
g.compute    = lower(g.compute);
g.AXES_FONT  = 10;
g.TITLE_FONT = 14;

% testing arguments consistency
% -----------------------------
if (~ischar(g.title))
   error('Title must be a string.');
end

if (~isnumeric(g.winsize) || length(g.winsize)~=1 || g.winsize~=round(g.winsize))
   error('Value of winsize must be an integer number.');
elseif (g.winsize <= 0)
   error('Value of winsize must be positive.');
elseif (g.cycles == 0 && pow2(nextpow2(g.winsize)) ~= g.winsize)
   error('Value of winsize must be an integer power of two [1,2,4,8,16,...]');
elseif (g.winsize > g.frame)
   error('Value of winsize must be less than frame length.');
end

if (~isnumeric(g.timesout) || length(g.timesout)~=1 || g.timesout~=round(g.timesout))
   error('Value of timesout must be an integer number.');
elseif (g.timesout <= 0)
   error('Value of timesout must be positive.');
end
if (g.timesout > g.frame-g.winsize)
   g.timesout = g.frame-g.winsize;
   disp(['Value of timesout must be <= frame-winsize, timeout adjusted to ' int2str(g.timesout) ]);
end

if (~isnumeric(g.padratio) || length(g.padratio)~=1 || g.padratio~=round(g.padratio))
   error('Value of padratio must be an integer.');
elseif (g.padratio <= 0)
   error('Value of padratio must be positive.');
elseif (pow2(nextpow2(g.padratio)) ~= g.padratio)
   error('Value of padratio must be an integer power of two [1,2,4,8,16,...]');
end

if (~isnumeric(g.maxfreq) || length(g.maxfreq)~=1)
   error('Value of g.maxfreq must be a number.');
elseif (g.maxfreq <= 0)
   error('Value of g.maxfreq must be positive.');
elseif (g.maxfreq > Fs/2)
   fprintf('Warning: input value of g.maxfreq larger that Nyquist frequency %3.4 Hz\n\n',Fs/2);
end

if isempty(g.topovec)
   g.topovec = [];
elseif min(size(g.topovec))==1
   g.topovec = g.topovec(:);
   if size(g.topovec,1)~=2
      error('topovec must be a row or column vector.');
   end
end

if isempty(g.elocs)
   g.elocs = '';
elseif (~ischar(g.elocs)) && ~isstruct(g.elocs)
   error('Channel location file must be a valid text file.');
end

if (~isnumeric(g.alpha) || length(g.alpha)~=1)
   error('timef(): Value of g.alpha must be a number.\n');
elseif (round(g.naccu*g.alpha) < 2)
   fprintf('Value of g.alpha is out of the normal range [%g,0.5]\n',2/g.naccu);
   g.naccu = round(2/g.alpha);
   fprintf('  Increasing the number of bootstrap iterations to %d\n',g.naccu);
end
if g.alpha>0.5 || g.alpha<=0
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
   case { 'rad', 'ms', 'deg' },;
   otherwise error('Angleunit must be either ''rad'', ''deg'', or ''ms''');
end;    
switch g.type
   case { 'coher', 'phasecoher' 'phasecoher2' },;
   otherwise error('Type must be either ''coher'' or ''phasecoher''');
end;    
switch g.boottype
   case { 'times' 'timestrials' 'trials'},;
   otherwise error('Boot type must be either ''times'', ''trials'' or ''timestrials''');
end;    
if (~isnumeric(g.shuffle))
   error('Shuffle argument type must be numeric');
end
switch g.memory
   case { 'low', 'high' },;
   otherwise error('memory must be either ''low'' or ''high''');
end
if strcmp(g.memory, 'low') && ~strcmp(g.boottype, 'times')
   error(['Bootstrap type ''' g.boottype ''' cannot be used in low memory mode']);
end

switch g.compute
   case { 'matlab', 'c' },;
   otherwise error('compute must be either ''matlab'' or ''c''');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare 2 conditions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iscell(X)
	if length(X) ~= 2 || length(Y) ~= 2
		error('crossf: to compare conditions, X and Y input must be 2-elements cell arrays');
	end
	if ~strcmp(g.boottype, 'times')
		disp('crossf warning: The significance bootstrap type is irrelevant when comparing conditions');
	end
	for index = 1:2:length(vararginori)
		if index<=length(vararginori) % needed: if elemenets are deleted
			%if strcmp(vararginori{index}, 'alpha'), vararginori(index:index+1) = [];
			if strcmp(vararginori{index}, 'title'), vararginori(index:index+1) = []; 
			end
		end
	end
	if iscell(g.title) 
		if length(g.title) <= 2,
			g.title{3} = 'Condition 2 - condition 1';
		end
	else
		g.title = { 'Condition 1', 'Condition 2', 'Condition 2 - condition 1' };
	end
	
	fprintf('Running crossf on condition 1 *********************\n');
	fprintf('Note: If an out-of-memory error occurs, try reducing the\n');
	fprintf('      number of time points or number of frequencies\n');
	if ~strcmp(g.type, 'coher')
	   fprintf('Note: Type ''coher'' takes 3 times as much memory as other options!)\n');
        end
	figure; 
	subplot(1,3,1); title(g.title{1});
	if ~strcmp(g.type, 'coher')
		[R1,mbase,times,freqs,Rbootout1,Rangle1, savecoher1] = crossf(X{1}, Y{1}, ...
								frame, tlimits, Fs, varwin, 'savecoher', 1, 'title', ' ',vararginori{:});
	else
		[R1,mbase,times,freqs,Rbootout1,Rangle1, savecoher1, Tfx1, Tfy1] = crossf(X{1}, Y{1}, ...
								frame, tlimits, Fs, varwin, 'savecoher', 1,'title', ' ',vararginori{:});
	end
	R1 = R1.*exp(j*Rangle1); % output Rangle is in radians
	
	% Asking user for memory limitations
	% if ~strcmp(g.noinput, 'yes')
	%	  tmp = whos('Tfx1');
	%	  fprintf('This function will require an additional %d bytes, do you wish\n', ...
    %     tmp.bytes*6+size(savecoher1,1)*size(savecoher1,2)*g.naccu*8);
	%	  res = input('to continue (y/n) (use the ''noinput'' option to disable this message):', 's');
	%  	if res == 'n', return; end
	% end

	fprintf('\nRunning crossf on condition 2 *********************\n');
	subplot(1,3,2); title(g.title{2});
	if ~strcmp(g.type, 'coher')
		  [R2,mbase,times,freqs,Rbootout2,Rangle2, savecoher2] = crossf(X{2}, Y{2}, ...
								frame, tlimits, Fs, varwin,'savecoher', 1, 'title', ' ',vararginori{:});
	else
		  [R2,mbase,times,freqs,Rbootout2,Rangle2, savecoher2, Tfx2, Tfy2] = crossf(X{2}, Y{2}, ...
								frame, tlimits, Fs, varwin,'savecoher', 1, 'title', ' ',vararginori{:});
	end
	R2 = R2.*exp(j*Rangle2); % output Rangle is in radians

	subplot(1,3,3); title(g.title{3});
	if isnan(g.alpha)
		plotall(R2-R1, [], [], times, freqs, mbase,  find(freqs <= g.maxfreq), g);
	else 
		% accumulate coherence images (all arrays [nb_points * timesout * trials])
		% ---------------------------
		allsavedcoher = zeros(size(savecoher1,1), ...
                          size(savecoher1,2), ...
                          size(savecoher1,3)+size(savecoher2,3));
		allsavedcoher(:,:,1:size(savecoher1,3))     = savecoher1;
		allsavedcoher(:,:,size(savecoher1,3)+1:end) = savecoher2;
		clear savecoher1 savecoher2;
		
		if strcmp(g.type, 'coher')
			alltfx = zeros(size(Tfx1,1), size(Tfx2,2), size(Tfx1,3)+size(Tfx2,3));
			alltfx(:,:,1:size(Tfx1,3))     = Tfx1;
			alltfx(:,:,size(Tfx1,3)+1:end) = Tfx2;
			clear Tfx1 Tfx2;
			
			alltfy = zeros(size(Tfy1,1), size(Tfy2,2), size(Tfy1,3)+size(Tfy2,3));
			alltfy(:,:,1:size(Tfy1,3))   = Tfy1;
			alltfy(:,:,size(Tfy1,3)+1:end) = Tfy2;
			clear Tfy1 Tfy2;
		end
		
		coherimages = zeros(size(allsavedcoher,1), size(allsavedcoher,2), g.naccu);
		cond1trials = length(X{1})/g.frame;
		cond2trials = length(X{2})/g.frame;
		alltrials = [1:cond1trials+cond2trials];
		fprintf('Accumulating bootstrap:');
		
		% preprocess data
		% ---------------
		switch g.type
		 case 'coher', % take the square of alltfx and alltfy
		  alltfx = alltfx.^2;
		  alltfy = alltfy.^2;
		 case 'phasecoher', % normalize
		  allsavedcoher = allsavedcoher ./ abs(allsavedcoher);
		 case 'phasecoher2', % don't do anything
		end
		
		if strcmp(g.type, 'coher')
			[coherdiff coher1 coher2] = coher2conddiff( allsavedcoher, alltrials, ...
                                                        cond1trials, g.type, alltfx, alltfy);
		else
			[coherdiff coher1 coher2] = coher2conddiff( allsavedcoher, alltrials, ...
                                                        cond1trials, g.type);
		end
		%figure; g.alpha = NaN; & to check that the new images are the same as the original
		%subplot(1,3,1); plotall(coher1, [], [], times, freqs, mbase, find(freqs <= g.maxfreq), g);
		%subplot(1,3,2); plotall(coher2, [], [], times, freqs, mbase, find(freqs <= g.maxfreq), g);
		%return;

		for index=1:g.naccu
			if rem(index,10) == 0,  fprintf(' %d',index); end
			if rem(index,120) == 0, fprintf('\n'); end
			
			if strcmp(g.type, 'coher')
				coherimages(:,:,index) = coher2conddiff( allsavedcoher, shuffle(alltrials), ...
                                                        cond1trials, g.type, alltfx, alltfy);
			else
				coherimages(:,:,index) = coher2conddiff( allsavedcoher, shuffle(alltrials), ...
                                                        cond1trials, g.type);
			end
		end
		fprintf('\n');

		% create articially a Bootstrap object to compute significance
		Boot = bootinit( [], size(allsavedcoher,1), g.timesout, g.naccu, 0, g.baseboot, ...
														'noboottype', g.alpha, g.rboot);
		Boot.Coherboot.R = coherimages;
		Boot = bootcomppost(Boot, [], [], []);
		g.title = '';
		plotall(coherdiff, Boot.Coherboot.R, Boot.Rsignif, times, freqs, mbase, ...
														find(freqs <= g.maxfreq), g);
	end
	return; % ********************************** END PROCESSING TWO CONDITIONS
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shuffle trials if necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if g.shuffle ~= 0
   fprintf('x and y data trials being shuffled %d times\n',g.shuffle);
   XX = reshape(X, 1, frame, length(X)/g.frame);
   YY = Y;
   X = [];
   Y = [];
   for index = 1:g.shuffle
      XX = shuffle(XX,3);
      X = [X XX(:,:)];
      Y = [Y YY];
   end
end

% detrend over epochs (trials) if requested
% -----------------------------------------
switch g.detrep
case 'on'
   X = reshape(X, g.frame, length(X)/g.frame);
   X = X - mean(X,2)*ones(1, length(X(:))/g.frame);
   Y = reshape(Y, g.frame, length(Y)/g.frame);
   Y = Y - mean(Y,2)*ones(1, length(Y(:))/g.frame);
end;        

% time limits
wintime = 500*g.winsize/g.srate;
times = [g.tlimits(1)+wintime:(g.tlimits(2)-g.tlimits(1)-2*wintime)/(g.timesout-1):g.tlimits(2)-wintime];

%%%%%%%%%%
% baseline
%%%%%%%%%%
if ~isnan(g.baseline)
   baseln = find(times < g.baseline); % subtract means of pre-0 (centered) windows
   if isempty(baseln)
      baseln = 1:length(times); % use all times as baseline
      disp('Bootstrap baseline empty, using the whole epoch.');
   end
   baselength = length(baseln);
else
   baseln = 1:length(times); % use all times as baseline
   baselength = length(times); % used for bootstrap
end

%%%%%%%%%%%%%%%%%%%%
% Initialize objects
%%%%%%%%%%%%%%%%%%%%
tmpsaveall = (~isnan(g.alpha) & isnan(g.rboot) & strcmp(g.memory, 'high')) ...
                 | (strcmp(g.subitc, 'on') & strcmp(g.memory, 'high'));
trials = length(X)/g.frame;
if ~strcmp(g.compute, 'c')
   Tfx = tfinit(X, g.timesout, g.winsize, g.cycles, g.frame, g.padratio, g.detret, ...
									g.srate, g.maxfreq, g.subitc, g.type, g.cyclesfact, tmpsaveall);
   Tfy = tfinit(Y, g.timesout, g.winsize, g.cycles, g.frame, g.padratio, g.detret, ...
									g.srate, g.maxfreq, g.subitc, g.type, g.cyclesfact, tmpsaveall);
   Coher     = coherinit(Tfx.nb_points, trials, g.timesout, g.type);
   Coherboot = coherinit(Tfx.nb_points, trials, g.naccu   , g.type);
   Boot      = bootinit( Coherboot, Tfx.nb_points, g.timesout, g.naccu, baselength, ...
									g.baseboot, g.boottype, g.alpha, g.rboot);
   freqs = Tfx.freqs;
   dispf = find(freqs <= g.maxfreq);
   freqs = freqs(dispf);
else
   freqs = g.srate*g.cycles/g.winsize*[2:2/g.padratio:g.winsize]/2;
end
dispf     = find(Tfx.freqs <= g.maxfreq);

%-------------
% Reserve space
%-------------
% R  = zeros(tfx.nb_points,g.timesout);       % mean coherence
% RR = repmat(nan,tfx.nb_points,g.timesout); % initialize with nans
% Rboot = zeros(tfx.nb_points,g.naccu);  % summed bootstrap coher
% switch g.type
% case 'coher',
%    cumulXY = zeros(tfx.nb_points,g.timesout);
%    cumulXYboot = zeros(tfx.nb_points,g.naccu);
% end;        
% if g.bootsub > 0
%    Rboottrial = zeros(tfx.nb_points, g.timesout, g.bootsub); % summed bootstrap coher
%    cumulXYboottrial = zeros(tfx.nb_points, g.timesout, g.bootsub);
% end
% if ~isnan(g.alpha) & isnan(g.rboot)
%    tf.tmpalltimes = repmat(nan,tfx.nb_points,g.timesout);
% end
   
% --------------------
% Display text to user
% --------------------
fprintf('\nComputing Event-Related ');
switch g.type
    case 'phasecoher',  fprintf('Phase Coherence (ITC) images for %d trials.\n',length(X)/g.frame);
    case 'phasecoher2', fprintf('Phase Coherence 2 (ITC) images for %d trials.\n',length(X)/g.frame);
    case 'coher',       fprintf('Linear Coherence (ITC) images for %d trials.\n',length(X)/g.frame);
end
fprintf('The trial latency range is from %4.5g ms before to %4.5g ms after\n     the time-locking event.\n', g.tlimits(1),g.tlimits(2));
fprintf('The frequency range displayed will be %g-%g Hz.\n',min(freqs),g.maxfreq);
if ~isnan(g.baseline)
   if length(baseln) == length(times)
      fprintf('Using the full trial latency range as baseline.\n');
   else
      fprintf('Using trial latencies from %4.5g ms to %4.5g ms as baseline.\n', g.tlimits,g.baseline);
   end
else 
   fprintf('No baseline time range was specified.\n');	
end
if g.cycles==0
   fprintf('The data window size will be %d samples (%g ms).\n',g.winsize,2*wintime);
   fprintf('The FFT length will be %d samples\n',g.winsize*g.padratio);
else
   fprintf('The window size will be %2.3g cycles.\n',g.cycles);
   fprintf('The maximum window size will be %d samples (%g ms).\n',g.winsize,2*wintime);
end
fprintf('The window will be applied %d times\n',g.timesout);
fprintf('     with an average step size of %2.2g samples (%2.4g ms).\n', Tfx.stp,1000*Tfx.stp/g.srate);
fprintf('Results will be oversampled %d times.\n',g.padratio);
if ~isnan(g.alpha)
   fprintf('Bootstrap confidence limits will be computed based on alpha = %g\n', g.alpha);
else
   fprintf('Bootstrap confidence limits will NOT be computed.\n'); 
end
switch g.plotphase
case 'on', 
    if strcmp(g.angleunit,'deg')
       fprintf(['Coherence angles will be imaged in degrees.\n']);
    elseif strcmp(g.angleunit,'rad')
       fprintf(['Coherence angles will be imaged in radians.\n']);
    elseif strcmp(g.angleunit,'ms')
       fprintf(['Coherence angles will be imaged in ms.\n']);
    end
end
fprintf('\nProcessing trial (of %d): ',trials);

% firstboot = 1;
% Rn=zeros(trials,g.timesout);
% X = X(:)'; % make X and Y column vectors
% Y = Y(:)';
% tfy = tfx;

if strcmp(g.compute, 'c')
   % C PART
   filename = [ 'tmpcrossf' num2str(round(rand(1)*1000)) ];
   f = fopen([ filename '.in'], 'w');
   fwrite(f, tmpsaveall, 'int32');
   fwrite(f, g.detret, 'int32');
   fwrite(f, g.srate, 'int32');
   fwrite(f, g.maxfreq, 'int32');
   fwrite(f, g.padratio, 'int32');
   fwrite(f, g.cycles, 'int32');
   fwrite(f, g.winsize, 'int32');
   fwrite(f, g.timesout, 'int32');
   fwrite(f, g.subitc, 'int32');
   fwrite(f, g.type, 'int32');
   fwrite(f, trials, 'int32');
   fwrite(f, g.naccu, 'int32');
   fwrite(f, length(X), 'int32');
   fwrite(f, X, 'double');
   fwrite(f, Y, 'double');
   fclose(f);
   
   command = [ '!cppcrosff ' filename '.in ' filename '.out' ];
   eval(command);
   
   f = fopen([ filename '.out'], 'r');
   size1 = fread(f, 'int32', 1);
   size2 = fread(f, 'int32', 1);
   Rreal = fread(f, 'double', [size1 size2]);
   Rimg  = fread(f, 'double', [size1 size2]);
   Coher.R = Rreal + j*Rimg;
   Boot.Coherboot.R = [];
   Boot.Rsignif = [];
else
   % ------------------------
   % MATLAB PART
   % compute ITC if necessary
   % ------------------------
   if strcmp(g.subitc, 'on')
      for t=1:trials
         if rem(t,10) == 0,  fprintf(' %d',t); end
         if rem(t,120) == 0, fprintf('\n'); end
         Tfx = tfitc( Tfx, t, 1:g.timesout); 
         Tfy = tfitc( Tfy, t, 1:g.timesout); 
      end; 
	  fprintf('\n');
      Tfx = tfitcpost( Tfx, trials); 
      Tfy = tfitcpost( Tfy, trials); 
   end
   
   % ---------
   % Main loop
   % ---------
   if g.savecoher,
	   trialcoher = zeros(Tfx.nb_points, g.timesout, trials);   
   else  
	   trialcoher = [];
   end
   for t=1:trials
      if rem(t,10) == 0,  fprintf(' %d',t); end
      if rem(t,120) == 0, fprintf('\n'); end
      
      Tfx = tfcomp( Tfx, t, 1:g.timesout); 
      Tfy = tfcomp( Tfy, t, 1:g.timesout);
	  if g.savecoher
		  [Coher trialcoher(:,:,t)] = cohercomp( Coher, Tfx.tmpalltimes, ...
                                             Tfy.tmpalltimes, t, 1:g.timesout);      
      else
		  Coher = cohercomp( Coher, Tfx.tmpalltimes, Tfy.tmpalltimes, t, 1:g.timesout);      
	  end
	  
      Boot = bootcomp( Boot, Coher.Rn(t,:), Tfx.tmpalltimes, Tfy.tmpalltimes);
   end % t = trial
   [Boot Rbootout] = bootcomppost(Boot, Coher.Rn, Tfx.tmpall, Tfy.tmpall);
      % Note that the bootstrap thresholding is actually performed 
      %      in the display subfunction plotall()

   Coher  = cohercomppost(Coher, trials);
end

% ----------------------------------
% If coherence, perform the division
% ----------------------------------
% switch g.type
% case 'coher',
%   R = R ./ cumulXY;
%   if ~isnan(g.alpha) & isnan(g.rboot)
%      Rboot = Rboot ./ cumulXYboot;  
%   end
%   if g.bootsub > 0
%      Rboottrial = Rboottrial ./ cumulXYboottrial;
%   end
% case 'phasecoher',
%   Rn = sum(Rn, 1);
%   R = R ./ (ones(size(R,1),1)*Rn);               % coherence magnitude
%   if ~isnan(g.alpha) & isnan(g.rboot)
%      Rboot = Rboot / trials;  
%   end
%   if g.bootsub > 0
%      Rboottrial = Rboottrial / trials;
%   end
% end

% ----------------
% Compute baseline
% ----------------
mbase = mean(abs(Coher.R(:,baseln)'));     % mean baseline coherence magnitude

% ---------------
% Plot everything
% ---------------
plotall(Coher.R, Boot.Coherboot.R, Boot.Rsignif, times, freqs, mbase, dispf, g);

% --------------------------------------
% Convert output Rangle to degrees or ms - Disabled to keep original default: radians output
% --------------------------------------
% Rangle = angle(Coher.R); % returns radians
% if strcmp(g.angleunit,'ms')  % convert to ms
%    Rangle = (Rangle/(2*pi)).*repmat(1000./freqs(dispf)',1,length(times)); 
% elseif strcmp(g.angleunit,'deg')  % convert to deg
%    Rangle = Rangle*180/pi; % convert to degrees
% else % angleunit is 'rad'
%    % Rangle = Rangle;
% end
% Rangle(find(Rraw==0)) = 0; % mask for significance - set angle at non-signif coher points to 0

R = abs(Coher.R);
Rsignif = Boot.Rsignif;
Tfx = permute(Tfx.tmpall, [3 2 1]); % from [trials timesout nb_points] 
%                                     to   [nb_points timesout trials]
Tfy = permute(Tfy.tmpall, [3 2 1]);

return; % end crossf() *************************************************

%
% crossf() plotting functions
% ----------------------------------------------------------------------
function plotall(R, Rboot, Rsignif, times, freqs, mbase, dispf, g) 

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

%
% Compute cross-spectral angles
% -----------------------------
Rangle = angle(R); % returns radians

%
% Optionally convert Rangle to degrees or ms
% ------------------------------------------
if strcmp(g.angleunit,'ms')  % convert to ms
   Rangle = (Rangle/(2*pi)).*repmat(1000./freqs(dispf)',1,length(times)); 
   maxangle = max(max(abs(Rangle)));
elseif strcmp(g.angleunit,'deg')  % convert to degrees
   Rangle = Rangle*180/pi; % convert to degrees
   maxangle = 180; % use full-cycle plotting 
else
   maxangle = pi;  % radians
end

R = abs(R);

% if ~isnan(g.baseline)
% 	R = R - repmat(mbase',[1 g.timesout]); % remove baseline mean
% end

Rraw = R; % raw coherence (e.g., coherency) magnitude values output

if g.plot
   fprintf('\nNow plotting...\n');
   set(gcf,'DefaultAxesFontSize',g.AXES_FONT)
   colormap(jet(256));
   
   pos = get(gca,'position'); % plot relative to current axes
   q = [pos(1) pos(2) 0 0];
   s = [pos(3) pos(4) pos(3) pos(4)];
   axis('off')
end

switch lower(g.plotamp)
case 'on' 
   %
   % Image the coherence [% perturbations] 
   %
   RR = R;
   if ~isnan(g.alpha) % zero out (and 'green out') nonsignif. R values
      RR(find(RR < repmat(Rboot(:),[1 g.timesout]))) = 0;
      Rraw(find(repmat(Rsignif(:),[1,size(Rraw,2)])>=Rraw))=0;
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
   if ~isempty(g.maxamp)
	   caxis([-g.maxamp g.maxamp]);
   end
   tmpscale = caxis;
   
   hold on
   plot([0 0],[0 freqs(max(dispf))],'--m','LineWidth',g.linewidth)
   for i=1:length(g.marktimes)
      plot([g.marktimes(i) g.marktimes(i)],[0 freqs(max(dispf))],'--m','LineWidth',g.linewidth);
   end
   hold off
   set(h(6),'YTickLabel',[],'YTick',[])
   set(h(6),'XTickLabel',[],'XTick',[])
   %title('Event-Related Coherence')
   
   h(8) = axes('Position',[.95 ordinate1 .05 height].*s+q);
   cbar(h(8),151:300, [0 tmpscale(2)]); % use only positive colors (gyorv) 
   
   %
   % Plot delta-mean min and max coherence at each time point on bottom of image
   %
   h(10) = axes('Units','Normalized','Position',[.1 ordinate1-0.1 .8 .1].*s+q); 
                                                            % plot marginal means below
   Emax = max(R(dispf,:)); % mean coherence at each time point
   Emin = min(R(dispf,:)); % mean coherence at each time point
   plot(times,Emin, times, Emax, 'LineWidth',g.linewidth); hold on;
   plot([times(1) times(length(times))],[0 0],'LineWidth',0.7);
   plot([0 0],[-500 500],'--m','LineWidth',g.linewidth);
   for i=1:length(g.marktimes)
       plot([g.marktimes(i) g.marktimes(i)],[-500 500],'--m','LineWidth',g.linewidth);
   end
   if ~isnan(g.alpha) && strcmp(g.boottype, 'trials') 
       % plot bootstrap significance limits (base mean +/-)
      plot(times,mean(Rboot(dispf,:)),'g','LineWidth',g.linewidth); hold on;
      plot(times,mean(Rsignif(dispf,:)),'k:','LineWidth',g.linewidth);
      axis([min(times) max(times) 0 max([Emax(:)' Rsignif(:)'])*1.2])
   else
      axis([min(times) max(times) 0 max(Emax)*1.2])
   end
   tick = get(h(10),'YTick');
   set(h(10),'YTick',[tick(1) ; tick(length(tick))])
   set(h(10),'YAxisLocation','right')
   xlabel('Time (ms)')
   ylabel('coh.')
   
   %
   % Plot mean baseline coherence at each freq on left side of image
   %
   h(11) = axes('Units','Normalized','Position',[0 ordinate1 .1 height].*s+q); 
                                                            % plot mean spectrum
   E = abs(mbase(dispf)); % baseline mean coherence at each frequency
   plot(freqs(dispf),E,'LineWidth',g.linewidth); % plot mbase
   if ~isnan(g.alpha) % plot bootstrap significance limits (base mean +/-)
      hold on
      % plot(freqs(dispf),Rboot(:,dispf)+[E;E],'g','LineWidth',g.linewidth);
      plot(freqs(dispf),mean(Rboot  (dispf,:),2),'g','LineWidth',g.linewidth);
      plot(freqs(dispf),mean(Rsignif(dispf,:),2),'k:','LineWidth',g.linewidth);
      axis([freqs(1) freqs(max(dispf)) 0 max([E Rsignif(:)'])*1.2]);
   else             % plot marginal mean coherence only
      if ~isnan(max(E))
         axis([freqs(1) freqs(max(dispf)) 0 max(E)*1.2]);
      end
   end
   
   tick = get(h(11),'YTick');
   set(h(11),'YTick',[tick(1) ; tick(length(tick))])
   set(h(11),'View',[90 90])
   xlabel('Freq. (Hz)')
   ylabel('coh.')
end

switch lower(g.plotphase)
case 'on'
   %
   % Plot coherence phase lags in bottom panel
   %
   h(13) = axes('Units','Normalized','Position',[.1 ordinate2 .8 height].*s+q);
   Rangle(find(Rraw==0)) = 0; % when plotting, mask for significance 
                              % = set angle at non-signif coher points to 0
   
   imagesc(times,freqs(dispf),Rangle(dispf,:),[-maxangle maxangle]); % plot the 
   hold on                                             % coherence phase angles
   plot([0 0],[0 freqs(max(dispf))],'--m','LineWidth',g.linewidth); % zero-time line
   for i=1:length(g.marktimes)
      plot([g.marktimes(i) g.marktimes(i)],[0 freqs(max(dispf))],'--m','LineWidth',g.linewidth);
   end
   
   ylabel('Freq. (Hz)')
   xlabel('Time (ms)')
   
   h(14)=axes('Position',[.95 ordinate2 .05 height].*s+q);
   cbar(h(14),0,[-maxangle maxangle]); % two-sided colorbar
end

if g.plot
	try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end
    if (length(g.title) > 0) % plot title
        if h(6) ~= 0, axes(h(6)); else axes(h(13)); end
        %h = subplot('Position',[0 0  1 1].*s+q, 'Visible','Off');               
        %h(13) = text(-.05,1.01,g.title);
        h(13) = title(g.title);
        %set(h(13),'VerticalAlignment','bottom')
        %set(h(13),'HorizontalAlignment','left')
        set(h(13),'FontSize',g.TITLE_FONT);
    end
   %
   %%%%%%%%%%%%%%% plot topoplot() %%%%%%%%%%%%%%%%%%%%%%%
   %
   if (~isempty(g.topovec))
      h(15) = subplot('Position',[-.1 .43 .2 .14].*s+q);
      if size(g.topovec,2) <= 2
         topoplot(g.topovec(1),g.elocs,'electrodes','off', ...
            'style', 'blank', 'emarkersize1chan', 10, 'chaninfo', g.chaninfo);
      else
         topoplot(g.topovec(1,:),g.elocs,'electrodes','off', 'chaninfo', g.chaninfo);
      end
      axis('square')
      
      h(16) = subplot('Position',[.9 .43 .2 .14].*s+q);
      if size(g.topovec,2) <= 2
         topoplot(g.topovec(2),g.elocs,'electrodes','off', ...
            'style', 'blank', 'emarkersize1chan', 10, 'chaninfo', g.chaninfo);
      else
         topoplot(g.topovec(2,:),g.elocs,'electrodes','off', 'chaninfo', g.chaninfo);
      end
      axis('square')
   end
   
   axcopy(gcf);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  TIME FREQUENCY   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for time freq initialisation
% -------------------------------------
function Tf = tfinit(X, timesout, winsize, ...
   cycles, frame, padratio, detret, srate, maxfreq, subitc, type, cyclesfact, saveall);
Tf.X         = X(:)'; % make X column vectors
Tf.winsize   = winsize;
Tf.cycles    = cycles;
Tf.frame     = frame;
Tf.padratio  = padratio;
Tf.detret    = detret;
Tf.stp       = (frame-winsize)/(timesout-1);
Tf.subitc    = subitc; % for ITC
Tf.type      = type; % for ITC
Tf.saveall   = saveall;
if (Tf.cycles == 0) %%%%%%%%%%%%%% constant window-length FFTs %%%%%%%%%%%%%%%%
   % Tf.freqs = srate/winsize*[1:2/padratio:winsize]/2; % incorect for padratio > 2
   Tf.freqs = linspace(0, srate/2, length([1:2/padratio:winsize])+1);
   Tf.freqs = Tf.freqs(2:end);
   Tf.win   = hanning(winsize);
   Tf.nb_points = padratio*winsize/2;   
else % %%%%%%%%%%%%%%%%%% Constant-Q (wavelet) DFTs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Tf.freqs = srate*cycles/winsize*[2:2/padratio:winsize]/2;
   Tf.win = dftfilt(winsize,maxfreq/srate,cycles,padratio,cyclesfact);
   Tf.nb_points = size(Tf.win,2);
end
Tf.tmpalltimes = zeros(Tf.nb_points, timesout);
trials = length(X)/frame;
if saveall
   Tf.tmpall     = repmat(nan,[trials timesout Tf.nb_points]);
else
   Tf.tmpall = [];
end
Tf.tmpallbool = zeros(trials,timesout);
Tf.ITCdone = 0;
if Tf.subitc
   Tf.ITC  = zeros(Tf.nb_points, timesout);
   switch Tf.type,
	   case { 'coher' 'phasecoher2' }
		Tf.ITCcumul  = zeros(Tf.nb_points, timesout);
   end
end

% function for itc
% ----------------
function [Tf, itcvals] = tfitc(Tf, trials, times);
Tf = tfcomp(Tf, trials, times);
switch Tf.type
   case 'coher',
      Tf.ITC(:,times)      = Tf.ITC(:,times) + Tf.tmpalltimes; % complex coher.
      Tf.ITCcumul(:,times) = Tf.ITCcumul(:,times)+abs(Tf.tmpalltimes).^2;
   case 'phasecoher2',
      Tf.ITC(:,times)      = Tf.ITC(:,times) + Tf.tmpalltimes; % complex coher.
      Tf.ITCcumul(:,times) = Tf.ITCcumul(:,times)+abs(Tf.tmpalltimes);
   case 'phasecoher',
      Tf.ITC(:,times)      = Tf.ITC(:,times) + Tf.tmpalltimes ./ abs(Tf.tmpalltimes); 
                                                            % complex coher.
end % ~any(isnan())
return;

function [Tf, itcvals] = tfitcpost(Tf, trials);
switch Tf.type
   case 'coher',       Tf.ITC = Tf.ITC ./ sqrt(trials * Tf.ITCcumul);
   case 'phasecoher2', Tf.ITC = Tf.ITC ./ Tf.ITCcumul;
   case 'phasecoher',  Tf.ITC = Tf.ITC / trials; % complex coher.
end % ~any(isnan())

if Tf.saveall
  Tf.ITC = transpose(Tf.ITC); % do not use ' otherwise conjugate

	%imagesc(abs(Tf.ITC)); colorbar; figure;
	%squeeze(Tf.tmpall(1,1,1:Tf.nb_points))
	%squeeze(Tf.ITC   (1,1,1:Tf.nb_points))
	%Tf.ITC = shiftdim(Tf.ITC, -1);

	Tf.ITC = repmat(shiftdim(Tf.ITC, -1), [trials 1 1]);
	Tf.tmpall = (Tf.tmpall - abs(Tf.tmpall) .* Tf.ITC) ./ abs(Tf.tmpall);

  %	for index = 1:trials
  %		imagesc(squeeze(abs(Tf.tmpall(index,:,:)))); drawnow; figure;
  %		Tf.tmpall(index,:,:) = (Tf.tmpall(index,:,:) - Tf.tmpall(index,:,:) .* Tf.ITC)./Tf.tmpall(index,:,:);
  %		imagesc(squeeze(abs(Tf.tmpall(index,:,:)))); drawnow;
  %		subplot(10,10, index); imagesc(squeeze(abs(Tf.tmpall(index,:,:)))); caxis([0 1]); drawnow;
  %	end
  %	squeeze(Tf.tmpall(1,1,1:Tf.nb_points))
  %	figure; axcopy;

end
Tf.ITCdone = 1;
return;

% function for time freq decomposition
% ------------------------------------
function [Tf, tmpX] = tfcomp(Tf, trials, times);
% tf is an structure containing all the information about the decomposition
for trial = trials
   for index = times
      if ~Tf.tmpallbool(trial, index) % already computed
         tmpX = Tf.X([1:Tf.winsize]+floor((index-1)*Tf.stp)+(trial-1)*Tf.frame);
         
         if ~any(isnan(tmpX)) % perform the decomposition
            tmpX = tmpX - mean(tmpX);
            switch Tf.detret, case 'on', 
               tmpX = detrend(tmpX); 
            end
            
            if Tf.cycles == 0 % use FFTs
               tmpX = Tf.win .* tmpX(:);
               tmpX = fft(tmpX,Tf.padratio*Tf.winsize);
               tmpX = tmpX(2:Tf.padratio*Tf.winsize/2+1);
            else 
               tmpX = transpose(Tf.win) * tmpX(:);
            end
         else
            tmpX = NaN;
         end
         if Tf.ITCdone
            tmpX = (tmpX - abs(tmpX) .* Tf.ITC(:,index)) ./ abs(tmpX);
         end
         Tf.tmpalltimes(:,index) = tmpX;
         if Tf.saveall
            Tf.tmpall(trial, index,:) = tmpX;
            Tf.tmpallbool(trial, index) = 1;
         end
	  else
		  Tf.tmpalltimes(:,index) = Tf.tmpall(trial, index,:);
     end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    COHERENCE    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for coherence initialisation
% -------------------------------------
function Coher = coherinit(nb_points, trials, timesout, type);
Coher.R  = zeros(nb_points,timesout);       % mean coherence
% Coher.RR = repmat(nan,nb_points,timesout); % initialize with nans
Coher.type = type;
Coher.Rn=zeros(trials,timesout);
switch type
 case 'coher',
  Coher.cumulX = zeros(nb_points,timesout);
  Coher.cumulY = zeros(nb_points,timesout);
 case 'phasecoher2',
  Coher.cumul  = zeros(nb_points,timesout);
end

% function for coherence calculation
% -------------------------------------
% function Coher = cohercomparray(Coher, tmpX, tmpY, trial);
% switch Coher.type
%   case 'coher',
%      Coher.R = Coher.R + tmpX.*conj(tmpY); % complex coher.
%      Coher.cumulXY = Coher.cumulXY + abs(tmpX).*abs(tmpY);
%   case 'phasecoher',
%      Coher.R = Coher.R + tmpX.*conj(tmpY) ./ (abs(tmpX).*abs(tmpY)); % complex coher.
%      Coher.Rn(trial,:) = 1;
% end % ~any(isnan())

function [Coher,tmptrialcoh] = cohercomp(Coher, tmpX, tmpY, trial, time);
tmptrialcoh = tmpX.*conj(tmpY);
switch Coher.type
   case 'coher',
      Coher.R(:,time) = Coher.R(:,time) + tmptrialcoh; % complex coher.
      Coher.cumulX(:,time) = Coher.cumulX(:,time) + abs(tmpX).^2;
      Coher.cumulY(:,time) = Coher.cumulY(:,time) + abs(tmpY).^2;
 case 'phasecoher2',
      Coher.R(:,time) = Coher.R(:,time) + tmptrialcoh; % complex coher.
      Coher.cumul(:,time) = Coher.cumul(:,time) + abs(tmptrialcoh);
   case 'phasecoher',
      Coher.R(:,time) = Coher.R(:,time) + tmptrialcoh ./ abs(tmptrialcoh); % complex coher.
	  %figure; imagesc(abs(tmpX.*conj(tmpY) ./ (abs(tmpX).*abs(tmpY))));
      Coher.Rn(trial,time) = Coher.Rn(trial,time)+1;
end % ~any(isnan())

% function for post coherence calculation
% ---------------------------------------
function Coher = cohercomppost(Coher, trials);
switch Coher.type
 case 'coher',
   Coher.R = Coher.R ./ sqrt(Coher.cumulX) ./ sqrt(Coher.cumulY);
 case 'phasecoher2',
   Coher.R = Coher.R ./ Coher.cumul;
 case 'phasecoher',
   Coher.Rn = sum(Coher.Rn, 1);
   Coher.R  = Coher.R ./ (ones(size(Coher.R,1),1)*Coher.Rn); % coherence magnitude
end

% function for 2 conditions coherence calculation
% -----------------------------------------------
function [coherimage, coherimage1, coherimage2] = coher2conddiff( allsavedcoher, alltrials, cond1trials, type, tfx, tfy);
	t1s = alltrials(1:cond1trials);
	t2s = alltrials(cond1trials+1:end);
	switch type
	 case 'coher',
	  coherimage1 = sum(allsavedcoher(:,:,t1s),3) ./ sqrt(sum(tfx(:,:,t1s),3)) ./ sqrt(sum(tfy(:,:,t1s),3));
	  coherimage2 = sum(allsavedcoher(:,:,t2s),3) ./ sqrt(sum(tfx(:,:,t2s),3)) ./ sqrt(sum(tfy(:,:,t1s),3));
	 case 'phasecoher2',
	  coherimage1 = sum(allsavedcoher(:,:,t1s),3) ./ sum(abs(allsavedcoher(:,:,t1s)),3);
	  coherimage2 = sum(allsavedcoher(:,:,t2s),3) ./ sum(abs(allsavedcoher(:,:,t2s)),3);
	 case 'phasecoher',
	  coherimage1 = sum(allsavedcoher(:,:,t1s),3) / cond1trials;
	  coherimage2 = sum(allsavedcoher(:,:,t2s),3) / (size(allsavedcoher,3)-cond1trials);
	end
	coherimage = coherimage2 - coherimage1;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BOOTSTRAP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for bootstrap initialisation
% -------------------------------------
function Boot = bootinit(Coherboot, nb_points, timesout, naccu, baselength, baseboot, boottype, alpha, rboot);
Boot.Rboot       = zeros(nb_points,naccu);  % summed bootstrap coher
Boot.boottype    = boottype;
Boot.baselength  = baselength;
Boot.baseboot    = baseboot;
Boot.Coherboot   = Coherboot;
Boot.naccu       = naccu;
Boot.alpha       = alpha;
Boot.rboot       = rboot;

% function for bootstrap computation
% ----------------------------------
function Boot = bootcomp(Boot, Rn, tmpalltimesx, tmpalltimesy);
if ~isnan(Boot.alpha) && isnan(Boot.rboot)
   if strcmp(Boot.boottype, 'times') % get g.naccu bootstrap estimates for each trial
      goodbasewins = find(Rn==1);
      if Boot.baseboot % use baseline windows only
         goodbasewins = find(goodbasewins<=Boot.baselength); 
      end
      ngdbasewins = length(goodbasewins);
      j=1;
      tmpsX = zeros(size(tmpalltimesx,1), Boot.naccu);
      tmpsY = zeros(size(tmpalltimesx,1), Boot.naccu);
      if ngdbasewins > 1
         while j<=Boot.naccu
            s = ceil(rand([1 2])*ngdbasewins); % random ints [1,g.timesout]
            s = goodbasewins(s);
            if ~any(isnan(tmpalltimesx(:,s(1)))) & ~any(isnan(tmpalltimesy(:,s(2))))
               tmpsX(:,j) = tmpalltimesx(:,s(1));
               tmpsY(:,j) = tmpalltimesy(:,s(2));
               j = j+1;
            end
         end
         Boot.Coherboot = cohercomp(Boot.Coherboot, tmpsX, tmpsY, 1, 1:Boot.naccu);
      end
   end
end

% handle other trial bootstrap types
% ----------------------------------
function [Boot, Rbootout] = bootcomppost(Boot, allRn, alltmpsX, alltmpsY);
trials    = size(alltmpsX, 1);
times     = size(alltmpsX, 2);
nb_points = size(alltmpsX, 3);
if ~isnan(Boot.alpha) && isnan(Boot.rboot)
   if strcmp(Boot.boottype, 'trials') % get g.naccu bootstrap estimates for each trial
      fprintf('\nProcessing trial bootstrap (of %d):',times(end));
      tmpsX = zeros(size(alltmpsX,3), Boot.naccu);
      tmpsY = zeros(size(alltmpsY,3), Boot.naccu );
      Boot.fullcoherboot = zeros(nb_points, Boot.naccu, times);
      
      for index=1:times
         if rem(index,10) == 0,  fprintf(' %d',index); end
         if rem(index,120) == 0, fprintf('\n'); end
         for allt=1:trials
            j=1;
            while j<=Boot.naccu
               t = ceil(rand([1 2])*trials); % random ints [1,g.timesout]
               if (allRn(t(1),index) == 1) && (allRn(t(2),index) == 1)
                  tmpsX(:,j) = squeeze(alltmpsX(t(1),index,:));
                  tmpsY(:,j) = squeeze(alltmpsY(t(2),index,:));
                  j = j+1;
               end
            end
            Boot.Coherboot = cohercomp(Boot.Coherboot, tmpsX, tmpsY, 1, 1:Boot.naccu);
         end
         Boot.Coherboot = cohercomppost(Boot.Coherboot);  % CHECK IF NECSSARY FOR ALL BOOT TYPE
         Boot.fullcoherboot(:,:,index) = Boot.Coherboot.R; 
         Boot.Coherboot = coherinit(nb_points, trials, Boot.naccu, Boot.Coherboot.type);
      end
      Boot.Coherboot.R = Boot.fullcoherboot;
      Boot = rmfield(Boot, 'fullcoherboot');
   elseif strcmp(Boot.boottype, 'timestrials') % handle timestrials bootstrap
      fprintf('\nProcessing time and trial bootstrap (of %d):',trials);
      tmpsX = zeros(size(alltmpsX,3), Boot.naccu);
      tmpsY = zeros(size(alltmpsY,3), Boot.naccu );
      for allt=1:trials
         if rem(allt,10) == 0,  fprintf(' %d',allt); end
         if rem(allt,120) == 0, fprintf('\n'); end
         j=1;
         while j<=Boot.naccu
            t = ceil(rand([1 2])*trials); % random ints [1,g.timesout]
            goodbasewins = find((allRn(t(1),:) & allRn(t(2),:)) ==1);
            if Boot.baseboot % use baseline windows only
               goodbasewins = find(goodbasewins<=baselength); 
            end
            ngdbasewins = length(goodbasewins);
            
            if ngdbasewins>1
               s = ceil(rand([1 2])*ngdbasewins); % random ints [1,g.timesout]
               s=goodbasewins(s);
               
               if all(allRn(t(1),s(1)) == 1) && all(allRn(t(2),s(2)) == 1)
                  tmpsX(:,j) = squeeze(alltmpsX(t(1),s(1),:));
                  tmpsY(:,j) = squeeze(alltmpsY(t(2),s(2),:));
                  j = j+1;
               end
            end
         end
         Boot.Coherboot = cohercomp(Boot.Coherboot, tmpsX, tmpsY, 1, 1:Boot.naccu);
      end
      Boot.Coherboot = cohercomppost(Boot.Coherboot);
   elseif strcmp(Boot.boottype, 'times') % boottype is 'times'
      Boot.Coherboot = cohercomppost(Boot.Coherboot);
   end
end

% test if precomputed
if ~isnan(Boot.alpha) && isnan(Boot.rboot) % if bootstrap analysis included . . .
   % 'boottype'='times' or 'timestrials', size(R)=nb_points*naccu
   % 'boottype'='trials',                 size(R)=nb_points*naccu*times
   Boot.Coherboot.R = abs (Boot.Coherboot.R);
   Boot.Coherboot.R = sort(Boot.Coherboot.R,2);

   % compute bootstrap significance level
   i = round(Boot.naccu*Boot.alpha);
   Boot.Rsignif = mean(Boot.Coherboot.R(:,Boot.naccu-i+1:Boot.naccu),2); % significance levels for Rraw
   Boot.Coherboot.R = squeeze(mean(Boot.Coherboot.R(:,Boot.naccu-i+1:Boot.naccu),2));
   if size(Boot.Coherboot.R, 2) == 1
	   Rbootout(:,2) = Boot.Coherboot.R;
   else
	   Rbootout(:,:,2) = Boot.Coherboot.R;
   end
   % BEFORE
   %Rboot = [mean(Rboot(1:i,:)) ; mean(Rboot(g.naccu-i+1:g.naccu,:))];
elseif ~isnan(Boot.rboot)
	Boot.Coherboot.R = Boot.rboot;
	Boot.Rsignif     = Boot.rboot;
	Rbootout         = Boot.rboot;
else 
	Boot.Coherboot.R = [];
	Boot.Rsignif     = [];
	Rbootout         = [];
end % NOTE: above, mean ?????

function w = hanning(n)
if ~rem(n,2)
   w = .5*(1 - cos(2*pi*(1:n/2)'/(n+1)));
   w = [w; w(end:-1:1)];
else
   w = .5*(1 - cos(2*pi*(1:(n+1)/2)'/(n+1)));
   w = [w; w(end-1:-1:1)];
end
