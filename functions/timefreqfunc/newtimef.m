% timef() - Returns estimates and plots of event-related (log) spectral
%           perturbation (ERSP) and inter-trial coherence (ITC) changes 
%           across event-related trials (epochs) of a single input time series. 
%         * Uses either fixed-window, zero-padded FFTs (fastest), wavelet
%           0-padded DFTs (both Hanning-tapered), OR multitaper spectra ('mtaper').
%         * For the wavelet and FFT methods, output frequency spacing 
%           is the lowest frequency ('srate'/'winsize') divided by 'padratio'.
%           NaN input values (such as returned by eventlock()) are ignored.
%         * If 'alpha' is given, then bootstrap statistics are computed 
%           (from a distribution of 'naccu' surrogate data trials) and 
%           non-significant features of the output plots are zeroed out 
%           (i.e., plotted in green). 
%         * Given a 'topovec' topo vector and 'elocs' electrode location file,
%           the figure also shows a topoplot() of the specified scalp map.
%         * Note: Left-click on subplots to view and zoom in separate windows.
% Usage: 
%        >> [ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
%                timef(data,frames,tlimits,srate,cycles,...
%                        'key1',value1,'key2',value2, ... );        
% NOTE:                                        
%        >> timef details  % scrolls more detailed information about timef
%
% Required inputs:     Value                                 {default}
%       data        = Single-channel data vector (1,frames*ntrials) (required)
%                     2-D array (frames,trials) or 3-D array (1,frames,trials)
%       frames      = Frames per trial. Ignored if data is 2-D or 3-D.  {750}
%       tlimits     = [mintime maxtime] (ms) Epoch time limits {[-1000 2000]}
%       srate       = data sampling rate (Hz)                 {250}
%       cycles      = is 0 -> Use FFTs (with constant window length) {0}
%                     is >0 -> Number of cycles in each analysis wavelet 
%                     is [wavcycles fact] -> wavelet cycles increase with frequency 
%                     starting at wavcyle (0<fact<1, fact=1 no increase, fact=0
%                     same as running FFT).
%                     OR multitaper decomposition (with 'mtaper').
%
%    Optional Inter-Irial Coherence Type:
%       'type'      = ['coher'|'phasecoher'] Compute either linear coherence 
%                      ('coher') or phase coherence ('phasecoher') also known
%                      as the phase coupling factor           {'phasecoher'}.
%
%    Optional Detrending:
%       'detret'    = ['on'|'off'], Detrend data in time.       {'off'}
%       'detrep'    = ['on'|'off'], Detrend data across trialsk {'off'}
%
%    Optional FFT/DFT Parameters:
%       'winsize'   = If cycles==0: data subwindow length (fastest, 2^n<frames);
%                     If cycles >0: *longest* window length to use. This
%                      determines the lowest output frequency  {~frames/8}
%       'timesout'  = Number of output times (int<frames-winframes) {200}
%       'padratio'  = FFT-length/winframes (2^k)                    {2}
%                      Multiplies the number of output frequencies by
%                      dividing their spacing. When cycles==0, frequency
%                      spacing is (low_freq/padratio).
%       'maxfreq'   = Maximum frequency (Hz) to plot (& to output, if cycles>0) 
%                      If cycles==0, all FFT frequencies are output. {50}
%       'baseline'  = Spectral baseline end-time (in ms). NaN imply that no
%                      baseline is used                              {0}
%       'powbase'   = Baseline spectrum to log-subtract. {def|NaN->from data}
%
%    Optional Bootstrap Parameters:
%       'alpha'     = If non-0, compute two-tailed bootstrap significance prob. 
%                      level. Show non-signif. output values as green.   {0}
%       'naccu'     = Number of bootstrap replications to accumulate     {200}
%       'baseboot'  = Bootstrap baseline subtract (0 -> use 'baseline';
%                                                  1 -> use whole trial) {0}
%       'boottype'  = ['times'|'timestrials'] shuffle time only or time and
%                     trials for computing bootstrap. Both options should
%                     return identical results {'times'}.
%       'condboot'  = ['abs'|'angle'|'complex'] for comparing 2 conditions,
%                     either subtract ITC absolute vales ('abs'), angles 
%                     ('angles') or complex values ('complex').     {'abs'}
%       'pboot'     = Bootstrap power limits (e.g., from timef())   {from data}
%       'rboot'     = Bootstrap ITC limits (e.g., from timef()). Note that both
%                     pboot and rboot must be provided to avoid recomputing
%                     surogate data.                                {from data}
%
%    Optional Scalp Map:
%       'topovec'   = Scalp topography (map) to plot                     {none}
%       'elocs'     = Electrode location file for scalp map   {no default}
%                     File should be ascii in format of  >> topoplot example   
%
%    Optional Plotting Parameters:
%       'ploterps'  = ['on'|'off'] Plot power spectral perturbations    {'on'} 
%       'plotitc'   = ['on'|'off'] Plot inter trial coherence            {'on'}
%       'plotphase' = ['on'|'off'] Plot phase in the inter trial coherence {'on'}
%       'itcmax'    = [real] set the ITC maximum for the scale       { auto }
%       'title'     = Optional figure title                              {none}
%       'marktimes' = Non-0 times to mark with a dotted vertical line (ms) {none}
%       'linewidth' = Line width for 'marktimes' traces (thick=2, thin=1) {2}
%       'axesfont'  = Axes text font size                                {10}
%       'titlefont' = Title text font size                               {8}
%       'vert'      = [times_vector] -> plot vertical dashed lines at specified times
%                     in ms.
%       'outputformat' = ['old'|'new'] for compatibility with script that used the old
%                        output format, set to 'old' (mbase in absolute amplitude (not
%                        dB) and real itc instead of complex itc). Default is 'new'.
%                     
% Outputs: 
%            ersp   = Matrix (nfreqs,timesout) of log spectral diffs. from baseline (dB) 
%            itc    = Matrix of inter-trial coherencies (nfreqs,timesout) (range: [0 1])
%          powbase  = Baseline power spectrum (removed for each window to compute the ersp)
%            times  = Vector of output times (subwindow centers) (in ms).
%            freqs  = Vector of frequency bin centers (in Hz).
%         erspboot  = Matrix (2,nfreqs) of [lower;upper] ERSP significance diffs.
%          itcboot  = Matrix (2,nfreqs) of [lower;upper] ITC thresholds (not diffs).
%           tfdata  = time frequency decomposition of the data (nfreqs,timesout,trials)
%
% Author: Arnaud Delorme, Sigurd Enghoff & Scott Makeig
%          CNL / Salk Institute 1998- | SCCN/INC, UCSD 2002-
%
% See also: crossf()
 
%123456789012345678901234567890123456789012345678901234567890123456789012
%    Optional Multitaper Parameters:
%       'mtaper'    = If [N W], performs multitaper decomposition. 
%                      (N is the time resolution and W the frequency resolution; 
%                      maximum taper number is 2NW-1). Overwrites 'winsize' and 'padratio'. 
%                     If [N W K], forces the use of K Slepian tapers (if possible).
%                      Phase is calculated using standard methods.
%                      The use of mutitaper with wavelets (cycles>0) is not 
%                      recommended (as multiwavelets are not implemented). 
%                      Uses Matlab functions DPSS, PMTM.      {no multitaper}
%

% Copyright (C) 1998 Arnaud Delorme, Sigurd Enghoff, Scott Makeig, Arnaud Delorme 
% CNL / Salk Institute 8/1/98-8/28/01
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
% Revision 1.5  2002/10/15 00:11:39  arno
% title diff bug
%
% Revision 1.4  2002/10/09 18:31:50  arno
% axcopy problem -> only try to execute it
%
% Revision 1.3  2002/10/07 22:44:46  arno
% creating outputs
%
% Revision 1.2  2002/10/02 00:35:53  arno
% update condstat, debug
%
% Revision 1.1  2002/10/01 16:06:55  arno
% Initial revision
%
% Revision 1.49  2002/08/14 21:07:05  arno
% hanning debug
%
% Revision 1.48  2002/08/14 21:02:19  arno
% implementing hanning funciton
%
% Revision 1.47  2002/08/12 01:47:29  arno
% color
%
% Revision 1.46  2002/08/11 22:30:20  arno
% color
%
% Revision 1.45  2002/08/09 22:29:44  arno
% implementing wavelet factor
%
% Revision 1.44  2002/07/22 14:28:56  arno
% debugging input baseline spectrum
%
% Revision 1.43  2002/07/11 18:19:30  arno
% header typo
%
% Revision 1.42  2002/07/11 15:27:28  arno
% debuging linear coherence
%
% Revision 1.41  2002/07/11 15:18:40  arno
% programing phase coherence 2
%
% Revision 1.40  2002/07/11 00:04:57  arno
% same
%
% Revision 1.39  2002/07/11 00:03:09  arno
% debugging itcmax
%
% Revision 1.38  2002/07/10 23:59:57  arno
% implement itcmax
%
% Revision 1.37  2002/07/10 21:46:03  arno
% adding phase argument
%
% Revision 1.36  2002/05/03 03:17:45  arno
% diffdlt 0.5 -> 1
%
% Revision 1.35  2002/05/01 22:25:21  arno
% no modif
%
% Revision 1.34  2002/05/01 21:23:24  arno
% no changes
%
% Revision 1.33  2002/04/30 04:27:54  arno
% trying to debug phsamp, still unsucessfull
%
% Revision 1.32  2002/04/29 15:13:50  scott
% debugging PA -sm
%
% Revision 1.31  2002/04/29 15:06:51  scott
% same -sm
%
% Revision 1.30  2002/04/29 15:03:31  scott
% same -sm
%
% Revision 1.29  2002/04/29 15:01:58  scott
% debugging -sm
%
% Revision 1.28  2002/04/29 14:58:56  scott
% debugging -sm
%
% Revision 1.27  2002/04/29 14:32:56  scott
% removing debugging statements -sm
%
% Revision 1.26  2002/04/29 14:32:31  scott
% adding ; -sm
%
% Revision 1.25  2002/04/29 14:29:53  scott
% debugging cumulX/PA -sm
%
% Revision 1.24  2002/04/29 14:24:13  scott
% debugging -sm
%
% Revision 1.23  2002/04/29 14:21:23  scott
% same -sm
%
% Revision 1.22  2002/04/29 14:15:58  scott
% insured cumulX sized -sm
%
% Revision 1.21  2002/04/29 14:12:33  scott
% used switches to test g.phsamp -sm
%
% Revision 1.20  2002/04/29 14:07:08  scott
% fixed typo -sm
%
% Revision 1.19  2002/04/29 14:02:34  scott
% made sure cumulX is computed -sm
%
% Revision 1.18  2002/04/29 13:57:49  scott
% modified PC->PA, 'phasecouple'->'phsamp', made output format (phs,amp,time) -sm
%
% Revision 1.17  2002/04/27 21:26:24  scott
% debugging PC -sm
%
% Revision 1.16  2002/04/27 21:19:02  scott
% debugging PC -sm
%
% Revision 1.15  2002/04/27 21:17:27  scott
% debugging PC -sm
%
% Revision 1.14  2002/04/27 21:13:57  scott
% added undocumented arg 'phasecouple',{'on'|'off'} -sm
%
% Revision 1.13  2002/04/25 02:56:03  arno
% redebugging topovec
%
% Revision 1.12  2002/04/25 02:54:33  arno
% improved topovec check
%
% Revision 1.11  2002/04/23 18:34:29  arno
% modified baseline way of computation
%
% Revision 1.10  2002/04/23 18:28:02  arno
% correcting coher computation
%
% Revision 1.9  2002/04/11 19:56:12  arno
% debuging baseboot -ad & lf
%
% Revision 1.8  2002/04/11 02:10:27  arno
% correcting typo
%
% Revision 1.7  2002/04/09 00:24:55  arno
% editing latest modifications
%
% Revision 1.6  2002/04/09 00:21:17  arno
% adding vertical bars
%
% Revision 1.5  2002/04/08 19:57:19  arno
% Editing, maxfreq-> Niquist
%
% Revision 1.4  2002/04/07 02:49:35  scott
% clarified hlpe message, changed default srate and winsize -sm
%
% Revision 1.3  2002/04/06 03:48:07  arno
% changing input for 1 channel for topoplot
%
% Revision 1.2  2002/04/06 03:40:58  arno
% modifying location file check
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 10-19-98 avoided division by zero (using MIN_ABS) -sm
% 10-19-98 improved usage message and commandline info printing -sm
% 10-19-98 made valid [] values for tvec and g.elocs -sm
% 04-01-99 added missing freq in freqs and plots, fixed log scaling bug -se & -tpj
% 06-29-99 fixed frequency indexing for constant-Q -se
% 08-24-99 reworked to handle NaN input values -sm
% 12-07-99 adjusted ERPtimes to plot ERP under ITC -sm
% 12-22-99 debugged ERPtimes, added BASE_BOOT -sm 
% 01-10-00 debugged BASE_BOOT=0 -sm
% 02-28-00 added NOTE on formula derivation below -sm
% 03-16-00 added axcopy() feature -sm & tpj
% 04-16-00 added multiple marktimes loop -sm
% 04-20-00 fixed ITC cbar limits when spcified in input -sm
% 07-29-00 changed frequencies displayed msg -sm
% 10-12-00 fixed bug in freqs when cycles>0 -sm
% 02-07-01 fixed inconsistency in BASE_BOOT use -sm
% 08-28-01 matlab 'key' value arguments -ad
% 08-28-01 multitaper decomposition -ad
% 01-25-02 reformated help & license -ad 
% 03-08-02 debug & compare to old timef function -ad 
% 03-16-02 timeout automatically adjusted if too high -ad 
% 04-02-02 added 'coher' option -ad 

function [P,R,mbase,times,freqs,Pboot,Rboot,alltfX,PA] = timef( X, frame, tlimits, Fs, varwin, varargin);

% Note: PA is output of 'phsamp','on' 

%varwin,winsize,g.timesout,g.padratio,g.maxfreq,g.topovec,g.elocs,g.alpha,g.marktimes,g.powbase,g.pboot,g.rboot)

% ITC:   Normally, R = |Sum(Pxy)| / (Sum(|Pxx|)*Sum(|Pyy|)) is coherence.
%        But here, we consider    Phase(PPy) = 0 and |Pyy| = 1 -> Pxy = Pxx
%        Giving, R = |Sum(Pxx)|/Sum(|Pxx|), the inter-trial coherence (ITC)
%        Also called 'phase-locking factor' by Tallon-Baudry et al. (1996)

% Constants set here:
ERSP_CAXIS_LIMIT = 0;           % 0 -> use data limits; else positive value
                                % giving symmetric +/- caxis limits.
ITC_CAXIS_LIMIT  = 0;           % 0 -> use data limits; else positive value
                                % giving symmetric +/- caxis limits.
MIN_ABS          = 1e-8;        % avoid division by ~zero

% Commandline arg defaults:
DEFAULT_EPOCH	= 750;			% Frames per trial
DEFAULT_TIMLIM = [-1000 2000];	% Time range of g.frames (ms)
DEFAULT_FS	= 250;			% Sampling frequency (Hz)
DEFAULT_NWIN	= 200;			% Number of windows = horizontal resolution
DEFAULT_VARWIN	= 0;			% Fixed window length or fixed number of cycles.
								% =0: fix window length to that determined by nwin
								% >0: set window length equal to varwin cycles
								%     Bounded above by winsize, which determines
								%     the min. freq. to be computed.
DEFAULT_OVERSMP	= 2;			% Number of times to oversample frequencies 
DEFAULT_MAXFREQ = 50;			% Maximum frequency to display (Hz)
DEFAULT_TITLE	= '';			% Figure title
DEFAULT_ELOC    = 'chan.locs';	% Channel location file
DEFAULT_ALPHA   = NaN;			% Percentile of bins to keep
DEFAULT_MARKTIME= NaN;

% Font sizes:
AXES_FONT       = 10;           % axes text FontSize
TITLE_FONT      = 8;

if (nargin < 1)
	help timef
	return
end

if isstr(X) & strcmp(X,'details')
   more on
   help timefdetails
   more off
   return
end
if ~iscell(X)
    [X, frame] = reshapeX(X, frame);
    trials = size(X,2);
else 
    [X{1}, frame] = reshapeX(X{1}, frame);
    [X{2}, frame] = reshapeX(X{2}, frame);
    trials = size(X{1},2);
end;

if (nargin < 2)
	frame = DEFAULT_EPOCH;
elseif (~isnumeric(frame) | length(frame)~=1 | frame~=round(frame))
	error('Value of frames must be an integer.');
elseif (frame <= 0)
	error('Value of frames must be positive.');
end;

if (nargin < 3)
	tlimits = DEFAULT_TIMLIM;
elseif (~isnumeric(tlimits) | sum(size(tlimits))~=3)
	error('Value of tlimits must be a vector containing two numbers.');
elseif (tlimits(1) >= tlimits(2))
	error('tlimits interval must be ascending.');
end

if (nargin < 4)
	Fs = DEFAULT_FS;
elseif (~isnumeric(Fs) | length(Fs)~=1)
	error('Value of srate must be a number.');
elseif (Fs <= 0)
	error('Value of srate must be positive.');
end

if (nargin < 5)
	varwin = DEFAULT_VARWIN;
elseif (~isnumeric(varwin) | length(varwin)>2)
	error('Value of cycles must be a number.');
elseif (varwin < 0)
	error('Value of cycles must be zero or positive.');
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
g.cycles  = varwin(1);
if length(varwin)>1
	g.cyclesfact = varwin(2);
else 
	g.cyclesfact = 1;
end;
try, g.boottype;   catch, g.boottype = 'times'; end;
try, g.condboot;   catch, g.condboot = 'abs'; end;
try, g.title;      catch, g.title = DEFAULT_TITLE; end;
try, g.winsize;    catch, g.winsize = max(pow2(nextpow2(g.frame)-3),4); end;
try, g.pad;        catch, g.pad = max(pow2(nextpow2(g.winsize)),4); end;
try, g.timesout;   catch, g.timesout = DEFAULT_NWIN; end;
try, g.padratio;   catch, g.padratio = DEFAULT_OVERSMP; end;
try, g.maxfreq;    catch, g.maxfreq = DEFAULT_MAXFREQ; end;
try, g.topovec;    catch, g.topovec = []; end;
try, g.elocs;      catch, g.elocs = DEFAULT_ELOC; end;
try, g.alpha;      catch, g.alpha = DEFAULT_ALPHA; end;  
try, g.marktimes;  catch, g.marktimes = DEFAULT_MARKTIME; end;
try, g.powbase;    catch, g.powbase = NaN; end;
try, g.pboot;      catch, g.pboot = NaN; end;
try, g.rboot;      catch, g.rboot = NaN; end;
try, g.plotersp;   catch, g.plotersp = 'on'; end;
try, g.plotitc;    catch, g.plotitc  = 'on'; end;
try, g.detrep;     catch, g.detrep = 'off'; end;
try, g.detret;     catch, g.detret = 'off'; end;
try, g.baseline;   catch, g.baseline = 0; end;
try, g.baseboot;   catch, g.baseboot = 0; end;
try, g.linewidth;  catch, g.linewidth = 2; end;
try, g.naccu;      catch, g.naccu = 200; end;
try, g.mtaper;     catch, g.mtaper = []; end;
try, g.vert;       catch, g.vert = []; end;
try, g.type;       catch, g.type = 'phasecoher'; end;
try, g.phsamp;     catch, g.phsamp = 'off'; end;
try, g.plotphase;  catch, g.plotphase = 'on'; end;
try, g.outputformat;  catch, g.outputformat = 'new'; end;
try, g.itcmax;     catch, g.itcmax = []; end;
g.AXES_FONT       = AXES_FONT;           % axes text FontSize
g.TITLE_FONT      = TITLE_FONT;
g.ERSP_CAXIS_LIMIT = ERSP_CAXIS_LIMIT;         
g.ITC_CAXIS_LIMIT  = ITC_CAXIS_LIMIT;        

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
	error('Value of maxfreq must be a number.');
elseif (g.maxfreq <= 0)
	error('Value of maxfreq must be positive.');
elseif (g.maxfreq > Fs/2)
	fprintf(['Warning: value of maxfreq reduced to Nyquist rate' ...
		 ' (%3.2f)\n\n'], Fs/2);
	g.maxfreq = Fs/2;
end

if isempty(g.topovec)
	g.topovec = [];
	if isempty(g.elocs)
		error('Channel location file must be specified.');
	end;
end
if isempty(g.elocs)
	g.elocs = DEFAULT_ELOC;
elseif (~ischar(g.elocs)) & ~isstruct(g.elocs)
	error('Channel location file must be a valid text file.');
end
if ~strcmpi(g.boottype, 'times') & ~strcmpi(g.boottype, 'timestrials')
	error('Boottype must be either ''times'' or ''timestrials''.');
end;	
if ~strcmpi(g.condboot, 'abs') & ~strcmpi(g.condboot, 'angle') ...
		& ~strcmpi(g.condboot, 'complex')
	error('Condboot must be either ''abs'', ''angle'' or ''complex''.');
end;
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
if ~isnumeric(g.vert)
    error('vertical line(s) option must be a vector');
else
	if min(g.vert) < g.tlimits(1) | max(g.vert) > g.tlimits(2)
        error('vertical line(s) time out-of-bound');
	end;
end;

% Multitape not used any more
% ---------------------------
if ~isempty(g.mtaper) % mutitaper, inspired from Bijan Pesaran matlab function
  if length(g.mtaper) < 3
        %error('mtaper arguement must be [N W] or [N W K]');
    
    if g.mtaper(1) * g.mtaper(2) < 1
        error('mtaper 2 first arguments'' product must be higher than 1');
    end;
    if length(g.mtaper) == 2
        g.mtaper(3) = floor( 2*g.mtaper(2)*g.mtaper(1) - 1);
    end
    if length(g.mtaper) == 3
        if g.mtaper(3) > 2 * g.mtaper(1) * g.mtaper(2) -1
            error('mtaper number too high (maximum (2*N*W-1))');
        end;
    end
    disp(['Using ' num2str(g.mtaper(3)) ' tapers.']);
    NW = g.mtaper(1)*g.mtaper(2);   % product NW
    N  = g.mtaper(1)*g.srate;     
    [e,v] = dpss(N, NW, 'calc');
    e=e(:,1:g.mtaper(3));
    g.alltapers = e;
  else    
    g.alltapers = g.mtaper;
    disp('mtaper argument not [N W] or [N W K]; considering raw taper matrix');
  end;
  g.winsize = size(g.alltapers, 1);
  g.pad = max(pow2(nextpow2(g.winsize)),256); % pad*nextpow
  nfk = floor([0 g.maxfreq]./g.srate.*g.pad);
  g.padratio = 2*nfk(2)/g.winsize;
 
  %compute number of frequencies
  %nf = max(256, g.pad*2^nextpow2(g.winsize+1)); 
  %nfk = floor([0 g.maxfreq]./g.srate.*nf);
  
  %freqs = linspace( 0, g.maxfreq, diff(nfk)); % this also work in the case of a FFT
  
end;           

switch lower(g.plotphase)
    case { 'on', 'off' }, ;
    otherwise error('plotphase must be either on or off');
end;
switch lower(g.plotersp)
    case { 'on', 'off' }, ;
    otherwise error('plotersp must be either on or off');
end;
switch lower(g.plotitc)
    case { 'on', 'off' }, ;
    otherwise error('plotitc must be either on or off');
end;
switch lower(g.detrep)
    case { 'on', 'off' }, ;
    otherwise error('detrep must be either on or off');
end;
switch lower(g.detret)
    case { 'on', 'off' }, ;
    otherwise error('detret must be either on or off');
end;
switch lower(g.phsamp)
    case { 'on', 'off' }, ;
    otherwise error('phsamp must be either on or off');
end;
if ~isnumeric(g.linewidth)
    error('linewidth must be numeric');
end;
if ~isnumeric(g.naccu)
    error('naccu must be numeric');
end;
if ~isnumeric(g.baseline)
    error('baseline must be numeric');
end;
switch g.baseboot
    case {0,1}, ;
    otherwise, error('baseboot must be 0 or 1');
end;
switch g.type
    case { 'coher', 'phasecoher', 'phasecoher2' },;
    otherwise error('Type must be either ''coher'' or ''phasecoher''');
end;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare 2 conditions part
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iscell(X)
    vararginori = varargin;
	if length(X) ~= 2
		error('timef: to compare conditions, data must be 2-elements cell arrays');
	end;
	if ~strcmp(g.boottype, 'times')
		disp('timef warning: significance bootstrap type is only use for subcondition plots');
	end;
	
    % deal with titles
    % ----------------
    for index = 1:2:length(vararginori)
		if index<=length(vararginori) % needed: if elemenets are deleted
			if strcmp(vararginori{index}, 'title'), vararginori(index:index+1) = []; 
			end;
		end;
	end;
	if iscell(g.title) 
		if length(g.title) <= 2,
			g.title{3} = 'Condition 1 - condition 2';
		end;
	else
		g.title = { 'Condition 1', 'Condition 2', 'Condition 2 - condition 1' };
	end;
	
	fprintf('Running timef on condition 1 *********************\n');
	fprintf('Note: if an out-of-memory error occurs, try reducing the\n');
	fprintf('      number of time points or number of frequencies\n');
	fprintf('      (the ''coher'' options takes 3 times more memory than other options)\n');
    [P1,R1,mbase1,times,freqs,Pboot1,Rboot1,alltfX1] = newtimef( X{1}, frame, tlimits, Fs, varwin, ...
                                                      'plotitc', 'off', 'plotersp', 'off', vararginori{:});
    
	fprintf('\nRunning crossf on condition 2 *********************\n');
    [P2,R2,mbase2,times,freqs,Pboot2,Rboot2,alltfX2] = newtimef( X{2}, frame, tlimits, Fs, varwin,  ...
                                                      'plotitc', 'off', 'plotersp', 'off', vararginori{:});
    
    % recompute baselines for power
    % -----------------------------
    if ~isnan( g.baseline ) & ~isnan( mbase1 )
        mbase = (mbase1 + mbase2)/2;
        P1 = P1 + repmat(mbase1(1:size(P1,1))',[1 size(P1,2)]); 
        P2 = P2 + repmat(mbase2(1:size(P1,1))',[1 size(P1,2)]); 
        P1 = P1 - repmat(mbase (1:size(P1,1))',[1 size(P1,2)]); 
        P2 = P2 - repmat(mbase (1:size(P1,1))',[1 size(P1,2)]);        
        if ~isnan(g.alpha)
			Pboot1 = Pboot1 + repmat(mbase1(1:size(Pboot1,1))',[1 size(Pboot1,2) size(Pboot1,3)]); 
			Pboot2 = Pboot2 + repmat(mbase2(1:size(Pboot1,1))',[1 size(Pboot1,2) size(Pboot1,3)]); 
			Pboot1 = Pboot1 - repmat(mbase (1:size(Pboot1,1))',[1 size(Pboot1,2) size(Pboot1,3)]); 
			Pboot2 = Pboot2 - repmat(mbase (1:size(Pboot1,1))',[1 size(Pboot1,2) size(Pboot1,3)]);        
        end;
        fprintf('\nSubtracting common baseline\n');
    end;
	
    % plotting
    % --------
    g.titleall = g.title;
    figure; subplot(1,3,1); g.title = g.titleall{1}; 
    plottimef(P1, R1, Pboot1, Rboot1, mean(X{1},2), freqs, times, mbase, g);
    subplot(1,3,2); g.title = g.titleall{2}; 
	plottimef(P2, R2, Pboot2, Rboot2, mean(X{2},2), freqs, times, mbase, g);
    
	subplot(1,3,3); g.title = g.titleall{3};
    if isnan(g.alpha)
        Rdiff = R2-R1;
		plottimef(P1-P2, Rdiff, [], [], mean(X{1},2)-mean(X{2},2), freqs, times, mbase, g);
	else 		
		% preprocess data and run compstat
		% --------------------------------
		alltfX1power = alltfX1.*conj(alltfX1);
		alltfX2power = alltfX2.*conj(alltfX2);
		formula = {'log10(mean(arg1(:,:,X),3))'};
		switch g.type
		 case 'coher', % take the square of alltfx and alltfy first to speed up
		  formula = { formula{1} ['mean(arg2(:,:,X),3)./sqrt(sum(arg1(:,:,X),3))'] };
		  [resdiff resimages res1 res2] = condstat(formula, g.naccu, g.alpha, {'both' 'upper'}, { '' g.condboot}, ...
                                    { alltfX1power alltfX2power }, {alltfX1 alltfX2});
		 case 'phasecoher2', % normalize first to speed up
		  formula = { formula{1} ['sum(arg2(:,:,X),3)./sum(arg3(:,:,X),3)'] }; 
		  alltfX1abs = sqrt(alltfX1power); % these 2 lines can be suppressed 
		  alltfX2abs = sqrt(alltfX2power); % by inserting sqrt(arg1(:,:,X)) instead of arg3(:,:,X))
		  [resdiff resimages res1 res2] = condstat(formula, g.naccu, g.alpha, {'both' 'upper'}, { '' g.condboot}, ...
                                    { alltfX1power alltfX2power }, {alltfX1 alltfX2}, { alltfX1abs alltfX2abs });
		 case 'phasecoher',
		  alltfX1norm = alltfX1./sqrt(alltfX1.*conj(alltfX1));
		  alltfX2norm = alltfX2./sqrt(alltfX2.*conj(alltfX2)); % maybe have to suppress preprocessing -> lot of memory
		  formula = { formula{1} ['mean(arg2(:,:,X),3)'] }; 
		  [resdiff resimages res1 res2] = condstat(formula, g.naccu, g.alpha, {'both' 'both'}, { '' g.condboot}, ...
                                                   { alltfX1power alltfX2power }, { alltfX1norm alltfX2norm });
		end;
        
		% same as below: plottimef(P1-P2, R2-R1, 10*resimages{1}, resimages{2}, mean(X{1},2)-mean(X{2},2), freqs, times, mbase, g);
        plottimef(10*resdiff{1}, resdiff{2}, 10*resimages{1}, resimages{2}, ...
                  mean(X{1},2)-mean(X{2},2), freqs, times, mbase, g);
		R1 = res1{2};
		R2 = res2{2};
        Rdiff = resdiff{2};
        Pboot = { Pboot1 Pboot2 10*resimages{1} };
        Rboot = { Rboot1 Rboot2 resimages{2} };
	end;
    P = { P1 P2 P1-P2 };
    R = { R1 R2 Rdiff };
    if nargout >= 8, alltfX = { alltfX1 alltfX2 }; end;
    
	return; % ********************************** END FOR SEVERAL CONDITIONS
end;
 
%%%%%%%%%%%%%%%%%%%%%%
% display text to user (computation perfomed only for display)
%%%%%%%%%%%%%%%%%%%%%%
wintime = 1000/g.srate*(g.winsize/2); % (1000/g.srate)*(g.winsize/2);
times = [g.tlimits(1)+wintime:(g.tlimits(2)-g.tlimits(1)-2*wintime)/(g.timesout-1):g.tlimits(2)-wintime];
if (g.cycles == 0) % FFT
    freqs = g.srate/g.winsize*[1:2/g.padratio:g.winsize]/2;
else % wavelet
    freqs = g.srate*g.cycles/g.winsize*[2:2/g.padratio:g.winsize]/2;
end;
dispf = find(freqs <= g.maxfreq);
freqs = freqs(dispf);
stp = (g.frame-g.winsize)/(g.timesout-1);
fprintf('Computing Event-Related Spectral Perturbation (ERSP) and\n');
switch g.type
    case 'phasecoher',  fprintf('  Inter-Trial Phase Coherence (ITC) images based on %d trials\n',trials);
    case 'phasecoher2', fprintf('  Inter-Trial Phase Coherence 2 (ITC) images based on %d trials\n',trials);
    case 'coher',       fprintf('  Linear Inter-Trial Coherence (ITC) images based on %d trials\n',trials);
end;
fprintf('  of %d frames sampled at %g Hz.\n',g.frame,g.srate);
fprintf('Each trial contains samples from %d ms before to\n',g.tlimits(1));
fprintf('  %d ms after the timelocking event.\n',g.tlimits(2));
fprintf('The window size used is %d samples (%g ms) wide.\n',g.winsize,2*wintime);
fprintf('The window is applied %d times at an average step\n',g.timesout);
fprintf('  size of %g samples (%gms).\n',stp,1000*stp/g.srate);
fprintf('Results are oversampled %d times; the %d frequencies\n',g.padratio,length(dispf));
fprintf('  displayed are from %2.1f Hz to %3.1f Hz.\n',freqs(dispf(1)),freqs(dispf(end)));
if ~isnan(g.alpha)
  fprintf('Only significant values (bootstrap p<%g) will be colored;\n',g.alpha) 
  fprintf('  non-significant values will be plotted in green\n');
end
fprintf('\nOf %d trials total, processing trial:',trials);

% -----------------------------------------
% detrend over epochs (trials) if requested
% -----------------------------------------
if strcmpi(g.detrep, 'on')
    X = X - mean(X,2)*ones(1, length(X(:))/g.frame);
end;        

% ----------------------------------------------------
% compute time frequency decompositions, power and ITC
% ----------------------------------------------------
g.subitc = 'off';
fprintf('\nProcessing trial (of %d):',trials);
[alltfX freqs times R] = timefreq(X, g.srate, 'timesout', g.timesout, 'winsize', g.winsize, ...
                                'tlimits', g.tlimits, 'detrend', g.detret, 'itctype', ...
                                g.type, 'subitc', g.subitc, 'wavelet', g.cycles, 'padratio', g.padratio); 
nb_points = size(alltfX,1);
dispf     = find(freqs <= g.maxfreq);
freqs = freqs(dispf);
P  = mean(alltfX.*conj(alltfX), 3); % power

% ----------------
% phase amp option
% ----------------
if strcmpi(g.phsamp, 'on')
%  switch g.phsamp
%  case 'on'
% $$$     PA = zeros(size(P,1),size(P,1),g.timesout); % NB: (freqs,freqs,times)
% $$$ end                                             %       phs   amp
    %PA (freq x freq x time)
    PA(:,:,j) = PA(:,:,j)  + (tmpX ./ abs(tmpX)) * ((PP(:,j)))';
    % x-product: unit phase column
    % times amplitude row

    tmpcx(1,:,:) = cumulX; % allow ./ below
    for j=1:g.timesout
        PA(:,:,j) = PA(:,:,j) ./ repmat(PP(:,j)', [size(PP,1) 1]);
    end
end

% --------
% baseline
% --------
if ~isempty(find(times < g.baseline))
   baseln = find(times < g.baseline); % subtract means of pre-0 (centered) windows
else
   baseln = 1:length(times); % use all times as baseline
end
if ~isnan(g.alpha) & length(baseln)==0
  fprintf('timef(): no window centers in baseline (times<%g) - shorten (max) window length.\n', g.baseline)
  return
elseif ~isnan(g.alpha) & g.baseboot
  fprintf('   %d bootstrap windows in baseline (times<%g).\n',...
          g.baseline,length(baseln))
end
if isnan(g.powbase)
  fprintf('\nComputing the mean baseline spectrum\n');
  mbase = mean(P(:,baseln),2)';
else
  fprintf('\nUsing the input baseline spectrum\n');
  mbase = 10.^(g.powbase/10);
end
baselength = length(baseln);
if ~isnan( g.baseline ) & ~isnan( mbase )
    P = 10 * (log10(P) - repmat(log10(mbase(1:size(P,1)))',[1 g.timesout])); % convert to (10log10) dB
else
    P = 10 * log10(P);
end;

% ---------
% bootstrap
% ---------
if ~isnan(g.alpha) % if bootstrap analysis included . . .
    if ~isnan(g.pboot) & ~isnan(g.rboot)
        Rboot = g.rboot;
        Pboot = g.pboot;
	else	
		formulaout = { 'power' 'itc' };
		switch g.type
		 case 'coher',
		  formulainit = [ 'power   = zeros(' int2str(nb_points) ',' int2str(g.naccu) ');' ...
						  'itc     = zeros(' int2str(nb_points) ',' int2str(g.naccu) ');' ...
						  'cumul   = zeros(' int2str(nb_points) ',' int2str(g.naccu) ');' ];
		  formula     =   'power   = power + arg1.*conj(arg1); itc   = itc + arg2; cumul = cumul + arg2.*conj(arg2);';
		  formulapost = [ 'power   = power /' int2str(trials) ';' ...
						  'itc     = itc ./ sqrt(cumul) /' int2str(trials) ];
		 case 'phasecoher',
		  formulainit = [ 'power   = zeros(' int2str(nb_points) ',' int2str(g.naccu) ');' ...
						  'itc     = zeros(' int2str(nb_points) ',' int2str(g.naccu) ');' ];
		  formula     =   'power   = power + arg1.*conj(arg1); itc   = itc + arg2;';
		  formulapost = [ 'power   = power /' int2str(trials) ';' ...
						  'itc     = itc /' int2str(trials) ];
		 case 'phasecoher2',
		  formulainit = [ 'power   = zeros(' int2str(nb_points) ',' int2str(g.naccu) ');' ...
						  'itc     = zeros(' int2str(nb_points) ',' int2str(g.naccu) ');' ...
						  'cumul   = zeros(' int2str(nb_points) ',' int2str(g.naccu) ');' ];
		  formula     =   'power   = power + arg1.*conj(arg1); itc   = itc + arg2; cumul = cumul + sqrt(arg2.*conj(arg2));';
		  formulapost = [ 'power   = power /' int2str(trials) ';' ...
						  'itc     = itc ./ cumul;' ];
		end;
		if g.baseboot == 0, baselntmp = [];
		else                baselntmp = baseln; 
		end;
		resboot = bootstat(alltfX, alltfX./sqrt(alltfX.*conj(alltfX)), formula, 'boottype', g.boottype, ...
						   'formulapost', formulapost, 'formulainit', formulainit, ...
						   'formulaout', formulaout, 'bootside', {'both' 'upper'}, 'naccu', g.naccu, 'basevect', baselntmp );
		Pboot = resboot{1};
		Rboot = resboot{2};
		
        if ~isnan( g.baseline ) 
			Pboot = 10*(log10(Pboot) - repmat(log10(mbase)', [1 size(Pboot,2) size(Pboot,3)])); 
        else
            Pboot = 10 * log10(Pboot);
        end;  
	end;
else 
    Pboot = []; Rboot = [];
end

% --------
% plotting
% --------
ERP = mean(X,2);
plottimef(P, R, Pboot, Rboot, ERP, freqs, times, mbase, g);
if strcmpi(g.outputformat, 'old')
    R = abs(R); % convert coherence vector to magnitude
else 
    mbase = log10(mbase)*10;
end;
return;

% -----------------
% plotting function
% -----------------
function plottimef(P, R, Pboot, Rboot, ERP, freqs, times, mbase, g);
    %
    % compute ERP
    %
    ERPtimes = [g.tlimits(1):(g.tlimits(2)-g.tlimits(1))/(g.frame-1):g.tlimits(2)+0.000001];
    ERPindices = zeros(1, length(times));
    for ti=1:length(times)
        [tmp ERPindices(ti)] = min(abs(ERPtimes-times(ti)));
    end
    ERPtimes = ERPtimes(ERPindices); % subset of ERP frames on t/f window centers
    ERP = ERP(ERPindices);
    
    dispf = find(freqs <= g.maxfreq);
	if ~isreal(R)
		Rsign = sign(imag(R));
		R = abs(R); % convert coherence vector to magnitude
		setylim = 1;
    else 
		Rsign = ones(size(R));
		setylim = 0;
	end;
	switch lower(g.plotitc)
     case 'on',  
      switch lower(g.plotersp), 
       case 'on', ordinate1 = 0.67; ordinate2 = 0.1; height = 0.33; g.plot = 1;
       case 'off', ordinate2 = 0.1; height = 0.9; g.plot = 1;
      end;     
     case 'off', ordinate1 = 0.1; height = 0.9; 
      switch lower(g.plotersp), 
       case 'on', ordinate1 = 0.1; height = 0.9;  g.plot = 1;
       case 'off', g.plot = 0;
      end;     
    end;    

    if g.plot
        fprintf('\nNow plotting...\n');
        set(gcf,'DefaultAxesFontSize',g.AXES_FONT)
        colormap(jet(256));
        pos = get(gca,'position');
        q = [pos(1) pos(2) 0 0];
        s = [pos(3) pos(4) pos(3) pos(4)];
    end;

    switch lower(g.plotersp)
     case 'on' 
      %
      %%%%%%% image the ERSP %%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      
      h(1) = subplot('Position',[.1 ordinate1 .9 height].*s+q);
      
      PP = P;
      if ~isnan(g.alpha) % zero out nonsignif. power differences
          if size(PP,1) == size(Pboot,1) & size(PP,2) == size(Pboot,2)
              PP(find(PP > Pboot(:,:,1) & (PP < Pboot(:,:,2)))) = 0;
              Pboot = squeeze(mean(Pboot,2));
          else
              PP(find((PP > repmat(Pboot(:,1),[1 g.timesout])) ...
                  & (PP < repmat(Pboot(:,2),[1 g.timesout])))) = 0;
          end;
      end

      if g.ERSP_CAXIS_LIMIT == 0
          ersp_caxis = [-1 1]*1.1*max(max(abs(P(dispf,:))));
      else
          ersp_caxis = g.ERSP_CAXIS_LIMIT*[-1 1];
      end

      if ~isnan( g.baseline ) 
          imagesc(times,freqs(dispf),PP(dispf,:),ersp_caxis); 
      else
          imagesc(times,freqs(dispf),PP(dispf,:));
      end;
      
      hold on
      plot([0 0],[0 freqs(max(dispf))],'--m','LineWidth',g.linewidth); % plot time 0
      if ~isnan(g.marktimes) % plot marked time
          for mt = g.marktimes(:)'
              plot([mt mt],[0 freqs(max(dispf))],'--k','LineWidth',g.linewidth);
          end
      end
      hold off
      set(h(1),'YTickLabel',[],'YTick',[])
      set(h(1),'XTickLabel',[],'XTick',[])
      if ~isempty(g.vert)
          for index = 1:length(g.vert)
              line([g.vert(index), g.vert(index)], [min(freqs(dispf)) max(freqs(dispf))], 'linewidth', 1, 'color', 'm');
          end;
      end;

      h(2) = gca;
      h(3) = cbar('vert'); % ERSP colorbar axes
      set(h(2),'Position',[.1 ordinate1 .8 height].*s+q)
      set(h(3),'Position',[.95 ordinate1 .05 height].*s+q)
      title('ERSP (dB)')

      E = [min(P(dispf,:));max(P(dispf,:))];
      h(4) = subplot('Position',[.1 ordinate1-0.1 .8 .1].*s+q); % plot marginal ERSP means
                                                                % below the ERSP image
      plot(times,E,[0 0],...
           [min(E(1,:))-max(max(abs(E)))/3 max(E(2,:))+max(max(abs(E)))/3], ...
           '--m','LineWidth',g.linewidth)
      axis([min(times) max(times) ...
            min(E(1,:))-max(max(abs(E)))/3 max(E(2,:))+max(max(abs(E)))/3])
      
      tick = get(h(4),'YTick');
      set(h(4),'YTick',[tick(1) ; tick(end)])
      set(h(4),'YAxisLocation','right')
      set(h(4),'TickLength',[0.020 0.025]);
      xlabel('Time (ms)')
      ylabel('dB')

      E = 10 * log10(mbase(dispf));
      h(5) = subplot('Position',[0 ordinate1 .1 height].*s+q); % plot mean spectrum
                                                               % to left of ERSP image
      if ~isnan(g.alpha)
          plot(freqs(dispf),Pboot(dispf,:)'+[E;E],'LineWidth',g.linewidth)
      else
          plot(freqs(dispf),E,'LineWidth',g.linewidth)
      end

      axis([freqs(1) freqs(max(dispf)) min(E)-max(abs(E))/3 max(E)+max(abs(E))/3])
      tick = get(h(5),'YTick');
      if (length(tick)>1)
          set(h(5),'YTick',[tick(1) ; tick(end-1)])
      end
      set(h(5),'TickLength',[0.020 0.025]);
      set(h(5),'View',[90 90])
      xlabel('Frequency (Hz)')
      ylabel('dB')
    end;

    switch lower(g.plotitc)
     case 'on'
      %
      %%%%%%%%%%%% Image the ITC %%%%%%%%%%%%%%%%%%
      %
      h(6) = subplot('Position',[.1 ordinate2 .9 height].*s+q); % ITC image

      RR = R;
      if ~isnan(g.alpha)
          if size(RR,1) == size(Rboot,1) & size(RR,2) == size(Rboot,2)
              tmp = gcf;
			  if size(Rboot,3) == 2	 RR(find(RR > Rboot(:,:,1) & RR < Rboot(:,:,2))) = 0;			  
			  else                   RR(find(RR < Rboot)) = 0;				  
			  end;
			  Rboot = mean(Rboot(:,:,end),2);
		  else
              RR(find(RR < repmat(Rboot(:),[1 g.timesout]))) = 0;
          end;
      end

      if g.ITC_CAXIS_LIMIT == 0
          coh_caxis = min(max(max(R(dispf,:))),1)*[-1 1]; % 1 WAS 0.4 !
      else
          coh_caxis = g.ITC_CAXIS_LIMIT*[-1 1];
      end

      if exist('Rsign') & strcmp(g.plotphase, 'on')
          imagesc(times,freqs(dispf),Rsign(dispf,:).*RR(dispf,:),coh_caxis); % <---
      else
          imagesc(times,freqs(dispf),RR(dispf,:),coh_caxis); % <---
      end
      if ~isempty(g.itcmax)
          caxis([-g.itcmax g.itcmax]);
      end;
      tmpcaxis = caxis;

      hold on
      plot([0 0],[0 freqs(max(dispf))],'--m','LineWidth',g.linewidth);
      if ~isnan(g.marktimes)
          for mt = g.marktimes(:)'
              plot([mt mt],[0 freqs(max(dispf))],'--k','LineWidth',g.linewidth);
          end
      end
      hold off
      set(h(6),'YTickLabel',[],'YTick',[])
      set(h(6),'XTickLabel',[],'XTick',[])
      if ~isempty(g.vert)
          for index = 1:length(g.vert)
              line([g.vert(index), g.vert(index)], [min(freqs(dispf)) max(freqs(dispf))], 'linewidth', 1, 'color', 'm');
          end;
      end;

      h(7) = gca;
      h(8) = cbar('vert');
      h(9) = get(h(8),'Children');
      set(h(7),'Position',[.1 ordinate2 .8 height].*s+q)
      set(h(8),'Position',[.95 ordinate2 .05 height].*s+q)
	  if setylim
		  set(h(8),'YLim',[0 tmpcaxis(2)]); 
      end;
	  title('ITC')

      %
      %%%%% plot the ERP below the ITC image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      % E = mean(R(dispf,:));

      ERPmax = max(ERP);
      ERPmin = min(ERP);
      ERPmax = ERPmax + 0.1*(ERPmax-ERPmin);
      ERPmin = ERPmin - 0.1*(ERPmax-ERPmin);
      h(10) = subplot('Position',[.1 ordinate2-0.1 .8 .1].*s+q); % ERP

      if ~isnan(g.alpha)

          % plot(times,E,[times(1) times(length(times))],...
          %       mean(Rboot(dispf))*[1 1],[0 0],...
          %    [min(E)-max(E)/3 max(E)+max(E)/3],'--m','LineWidth',g.linewidth)
          % axis([min(times) max(times) min(E)-max(E)/3 ...
          %       max([E mean(Rboot(dispf))])+max(E)/3]);

          plot(ERPtimes,ERP,...
               [times(1) times(length(times))],[0 0],...
               [0 0],[ERPmin ERPmin],'--m','LineWidth',g.linewidth);
          axis([min(ERPtimes) max(ERPtimes) ERPmin ERPmax]);
      else

          % plot(times,E,[0 0],[min(E)-max(E)/3 max(E)+max(E)/3],'--m',...
          %      'LineWidth',g.linewidth)
          % axis([min(times) max(times) min(E)-max(E)/3 max(E)+max(E)/3]);

          plot(ERPtimes,ERP,...
               [times(1) times(length(times))],[0 0],...
               [0 0],[ERPmin ERPmax],'--m','LineWidth',g.linewidth);
          axis([min(ERPtimes) max(ERPtimes) ERPmin ERPmax]);
      end

      tick = get(h(10),'YTick');
      set(h(10),'YTick',[tick(1) ; tick(end)])
      set(h(10),'TickLength',[0.02 0.025]);
      set(h(10),'YAxisLocation','right')
      xlabel('Time (ms)')
      ylabel('uV')

      E = mean(R(dispf,:)');
      
      h(11) = subplot('Position',[0 ordinate2 .1 height].*s+q); % plot the marginal mean
                                                                % ITC left of the ITC image
      if ~isnan(g.alpha)
          plot(freqs(dispf),E,freqs(dispf),Rboot(dispf),'LineWidth',g.linewidth)
          axis([freqs(1) freqs(max(dispf)) 0 max([E Rboot(dispf)'])+max(E)/3])
      else
          plot(freqs(dispf),E,'LineWidth',g.linewidth)
          axis([freqs(1) freqs(max(dispf)) min(E)-max(E)/3 max(E)+max(E)/3])
      end

      tick = get(h(11),'YTick');
      set(h(11),'YTick',[tick(1) ; tick(length(tick))])
      set(h(11),'View',[90 90])
      set(h(11),'TickLength',[0.020 0.025]);
      xlabel('Frequency (Hz)')
      ylabel('ERP')
      %
      %%%%%%%%%%%%%%% plot a topoplot() %%%%%%%%%%%%%%%%%%%%%%%
      %
      if (~isempty(g.topovec))
          h(12) = subplot('Position',[-.1 .43 .2 .14].*s+q);
          if length(g.topovec) == 1
              topoplot(g.topovec,g.elocs,'electrodes','off', ...
                       'style', 'blank', 'emarkersize1chan', 10);
          else
              topoplot(g.topovec,g.elocs,'electrodes','off');
          end;
          axis('square')
      end
    end; %switch

    if g.plot
        try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;
        if (length(g.title) > 0) & ~iscell(g.title)
            axes('Position',pos,'Visible','Off');               
            h(13) = text(-.05,1.01,g.title);
            set(h(13),'VerticalAlignment','bottom')     
            set(h(13),'HorizontalAlignment','left') 
            set(h(13),'FontSize',g.TITLE_FONT);
        end

        try, axcopy(gcf); catch, end;
    end;

% reshaping X
% -----------
function [X, frame] = reshapeX(X, frame)
    X = squeeze(X);
    if min(size(X)) == 1
        if (rem(length(X),frame) ~= 0)
            error('Length of data vector must be divisible by frames.');
        end
        X = reshape(X, frame, length(X)/frame);
    else
        frame = size(X,1);
    end 
