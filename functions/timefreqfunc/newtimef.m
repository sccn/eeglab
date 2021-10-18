% newtimef() - Return estimates and plots of mean event-related (log) spectral
%           perturbation (ERSP) and inter-trial coherence (ITC) events across
%           event-related trials (epochs) of a single input channel time series.
%
%         * Also can compute and statistically compare transforms for two time
%           series. Use this to compare ERSP and ITC means in two conditions.
%
%         * Uses either fixed-window, zero-padded FFTs (fastest), or wavelet
%           0-padded DFTs. FFT uses Hanning tapers; wavelets use (similar) Morlet
%           tapers.
%
%         * For the wavelet and FFT methods, output frequency spacing
%           is the lowest frequency ('srate'/'winsize') divided by 'padratio'.
%           NaN input values (such as returned by eventlock()) are ignored.
%
%         * If 'alpha' is given (see below), permutation statistics are computed
%           (from a distribution of 'naccu' surrogate data trials) and
%           non-significant features of the output plots are zeroed out
%           and plotted in green.
%
%         * Given a 'topovec' topo vector and 'elocs' electrode location file,
%           the figure also shows a topoplot() view of the specified scalp map.
%
%         * Note: Left-click on subplots to view and zoom in separate windows.
%
% Usage with single dataset:
%        >> [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = ...
%                  newtimef(data, frames, epochlim, srate, cycles,...
%                       'key1',value1, 'key2',value2, ... );
%
% Example to compare two condition (channel 1 EEG versus ALLEEG(2)):
%        >> [ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
%                  newtimef({EEG.data(1,:,:) ALLEEG(2).data(1,:,:)},
%                       EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, cycles);
% NOTE:
%        >> timef details  % presents more detailed argument information
%                          % Note: version timef() also computes multitaper transforms
%
% Required inputs:    Value                                 {default}
%       data        = Single-channel data vector (1,frames*ntrials), else 
%                     2-D array (frames,trials) or 3-D array (1,frames,trials).
%                     To compare two conditions (data1 versus data2), in place of 
%                     a single data matrix enter a cell array {data1 data2}
%       frames      = Frames per trial. Ignored if data are 2-D or 3-D.  {750}
%       tlimits     = [mintime maxtime] (ms).  Note that these are the time limits 
%                     of the data epochs themselves, NOT A SUB-WINDOW TO EXTRACT 
%                     FROM THE EPOCHS as is the case for pop_newtimef(). {[-1000 2000]}
%       Fs          = data sampling rate (Hz)  {default: read from icadefs.m or 250}
%       varwin      = [real] indicates the number of cycles for the time-frequency 
%                        decomposition {default: 0}
%                     If 0, use FFTs and Hanning window tapering.  
%                     If [real positive scalar], the number of cycles in each Morlet 
%                        wavelet, held constant across frequencies.
%                     If [cycles cycles(2)] wavelet cycles increase with 
%                        frequency beginning at cycles(1) and, if cycles(2) > 1, 
%                        increasing to cycles(2) at the upper frequency,
%                      If cycles(2) = 0, use same window size for all frequencies 
%                        (similar to FFT when cycles(1) = 1)
%                      If cycles(2) = 1, cycles do not increase (same as giving
%                         only one value for 'cycles'). This corresponds to a pure
%                         wavelet decomposition, same number of cycles at each frequency.
%                      If 0 < cycles(2) < 1, cycles increase linearly with frequency:
%                         from 0 --> FFT (same window width at all frequencies) 
%                         to 1 --> wavelet (same number of cycles at all frequencies).
%                     The exact number of cycles in the highest frequency window is 
%                     indicated in the command line output. Typical value: 'cycles', [3 0.5]
%
%    Optional inter-trial coherence (ITC) Type:
%       'itctype'   = ['coher'|'phasecoher'|'phasecoher2'] Compute either linear
%                     coherence ('coher') or phase coherence ('phasecoher').
%                     Originally called 'phase-locking factor' {default: 'phasecoher'}
%
%    Optional detrending:
%       'detrend'   = ['on'|'off'], Linearly detrend each data epoch   {'off'}
%       'rmerp'     = ['on'|'off'], Remove epoch mean from data epochs {'off'}
%
%    Optional FFT/DFT parameters:
%       'winsize'   = If cycles==0: data subwindow length (fastest, 2^n<frames);
%                     If cycles >0: The *longest* window length to use. This
%                     determines the lowest output frequency. Note: this parameter 
%                     is overwritten when the minimum frequency requires
%                     a longer time window {default: ~frames/8}
%       'timesout'  = Number of output times (int<frames-winframes). Enter a
%                     negative value [-S] to subsample original times by S.
%                     Enter an array to obtain spectral decomposition at
%                     specific times (Note: The algorithm finds the closest time
%                     point in data; this could give a slightly unevenly spaced
%                     time array                                    {default: 200}
%       'padratio'  = FFT-length/winframes (2^k)                    {default: 2}
%                     Multiplies the number of output frequencies by dividing
%                     their spacing (standard FFT padding). When cycles~=0,
%                     frequency spacing is divided by padratio.
%       'maxfreq'   = Maximum frequency (Hz) to plot (& to output, if cycles>0)
%                     If cycles==0, all FFT frequencies are output. {default: 50}
%                     DEPRECATED, use 'freqs' instead,and never both.
%       'freqs'     = [min max] frequency limits. {default [minfreq 50],
%                     minfreq being determined by the number of data points,
%                     cycles and sampling frequency.
%       'nfreqs'    = number of output frequencies. For FFT, closest computed
%                     frequency will be returned. Overwrite 'padratio' effects
%                     for wavelets. {default: use 'padratio'}
%       'freqscale' = ['log'|'linear'] frequency scale. {default: 'linear'}
%                     Note that for obtaining 'log' spaced freqs using FFT,
%                     closest correspondent frequencies in the 'linear' space
%                     are returned.
%       'verbose'   = ['on'|'off'] print text {'on'}
%       'subitc'    = ['on'|'off'] subtract stimulus locked Inter-Trial Coherence
%                     (ITC) from x and y. This computes an 'intrinsic' coherence
%                     of x and y not arising directly from common phase locking 
%                     to experimental events. See notes.    {default: 'off'}
%       'wletmethod' = ['dftfilt'|'dftfilt2'|'dftfilt3'] Wavelet type to use.
%                     'dftfilt2' -> Morlet-variant wavelets, or Hanning DFT.
%                     'dftfilt3' -> Morlet wavelets.  See the timefreq() function 
%                     for more detials {default: 'dftfilt3'}
%       'cycleinc'    ['linear'|'log'] mode of cycles increase when [min max] cycles 
%                     are provided in 'cycle' parameter. Applies only to 
%                     'wletmethod','dftfilt'  {default: 'linear'}
%       
%   Optional baseline parameters:
%       'baseline'  = Spectral baseline end-time (in ms). NaN --> no baseline is used. 
%                     A [min max] range may also be entered
%                     You may also enter one row per region for baseline
%                     e.g. [0 100; 300 400] considers the window 0 to 100 ms and
%                     300 to 400 ms This parameter validly defines all baseline types 
%                     below. Again, [NaN] Prevent baseline subtraction.
%                     {default: 0 -> all negative time values}. 
%       'powbase'   = Baseline spectrum to log-subtract {default|NaN -> from data}
%       'commonbase' = ['on'|'off'] use common baseline when comparing two 
%                     conditions {default: 'on'}.
%       'basenorm'  = ['on'|'off'] 'on' normalize baseline in the power spectral
%                     average; else 'off', divide by the average power across 
%                     trials at each frequency (gain model). {default: 'off'}
%       'trialbase' = ['on'|'off'|'full'] perform baseline (normalization or division 
%                     above in single trial instead of the trial average. Default
%                     if 'off'. 'full' is an option that perform single
%                     trial normalization (or simple division based on the
%                     'basenorm' input over the full trial length) before
%                     performing standard baseline removal. It has been
%                     shown to be less sensitive to noisy trials in Grandchamp R, 
%                     Delorme A. (2011) Single-trial normalization for event-related 
%                     spectral decomposition reduces sensitivity to noisy trials. 
%                     Front Psychol. 2:236.
%
%    Optional time warping parameter: 
%       'timewarp'  = [eventms matrix] Time-warp amplitude and phase time-
%                     courses(following time/freq transform but before 
%                     smoothing across trials). 'eventms' is a matrix 
%                     of size (all_trials,epoch_events) whose columns
%                     specify the epoch times (latencies) (in ms) at which 
%                     the same series of successive events occur in each 
%                     trial. If two data conditions, eventms should be 
%                     [eventms1;eventms2] --> all trials stacked vertically.
%      'timewarpms' = [warpms] optional vector of event times (latencies) (in ms) 
%                     to which the series of events should be warped.
%                     (Note: Epoch start and end should not be declared
%                     as eventms or warpms}. If 'warpms' is absent or [], 
%                     the median of each 'eventms' column will be used;
%                     If two datasets, the grand medians of the two are used.
%     'timewarpidx' = [plotidx] is an vector of indices telling which of 
%                     the time-warped 'eventms' columns (above) to show with 
%                     vertical lines. If undefined, all columns are plotted. 
%                     Overwrites the 'vert' argument (below) if any.
%
%    Optional permutation parameters:
%       'alpha'     = If non-0, compute two-tailed permutation significance 
%                      probability level. Show non-signif. output values 
%                      as green.                              {default: 0}
%       'mcorrect'  = ['none'|'fdr'] correction for multiple comparison
%                     'fdr' uses false detection rate (see function fdr()).
%                     Not available for condition comparisons. {default:'none'} 
%       'pcontour'  = ['on'|'off'] draw contour around significant regions
%                     instead of masking them. Not available for condition 
%                     comparisons. {default:'off'} 
%       'naccu'     = Number of permutation replications to accumulate {200}
%       'baseboot'  = permutation baseline subtract (1 -> use 'baseline';
%                                                    0 -> use whole trial
%                                            [min max] -> use time range) 
%                     You may also enter one row per region for baseline,
%                     e.g. [0 100; 300 400] considers the window 0 to 100 ms 
%                     and 300 to 400 ms. {default: 1}
%       'boottype'  = ['shuffle'|'rand'|'randall'] 'shuffle' -> shuffle times 
%                     and trials; 'rand' -> invert polarity of spectral data 
%                     (for ERSP) or randomize phase (for ITC); 'randall' -> 
%                     compute significances by accumulating random-polarity 
%                     inversions for each time/frequency point (slow!). Note
%                     that in the previous revision of this function, this
%                     method was called 'bootstrap' though it is actually 
%                     permutation {default: 'shuffle'}
%       'condboot'  = ['abs'|'angle'|'complex'] to compare two conditions,
%                     either subtract ITC absolute values ('abs'), angles
%                     ('angles'), or complex values ('complex'). {default: 'abs'}
%       'pboot'     = permutation power limits (e.g., from newtimef()) {def: from data}
%       'rboot'     = permutation ITC limits (e.g., from newtimef()). 
%                     Note: Both 'pboot' and 'rboot' must be provided to avoid 
%                     recomputing the surrogate data! {default: from data}
%
%    Optional Scalp Map:
%       'topovec'   = Scalp topography (map) to plot              {none}
%       'elocs'     = Electrode location file for scalp map       {none}
%                     Value should be a string array containing the path
%                     and name of the file.  For file format, see
%                         >> topoplot example
%       'chaninfo'    Passed to topoplot, if called.
%                     [struct] optional structure containing fields 
%                     'nosedir', 'plotrad', and/or 'chantype'. See these 
%                     field definitions above, below.
%                     {default: nosedir +X, plotrad 0.5, all channels}
%
%     Optional Plotting Parameters:
%       'scale'     = ['log'|'abs'] visualize power in log scale (dB) or absolute
%                     scale. {default: 'log'}
%       'plottype'  = ['image'|'curve'] plot time/frequency images or traces
%                     (curves, one curve per frequency). {default: 'image'}
%       'plotmean'  = ['on'|'off'] For 'curve' plots only. Average all
%                     frequencies given as input. {default: 'on'}
%       'highlightmode'  = ['background'|'bottom'] For 'curve' plots only,
%                     display significant time regions either in the plot background
%                     or under the curve.
%       'plotersp'  = ['on'|'off'] Plot power spectral perturbations    {'on'}
%       'plotitc'   = ['on'|'off'] Plot inter-trial coherence           {'on'}
%       'plotphasesign' = ['on'|'off'] Plot phase sign in the inter trial coherence {'on'}
%       'plotphaseonly' = ['on'|'off'] Plot ITC phase instead of ITC amplitude {'off'}
%       'erspmax'   = [real] set the ERSP max. For the color scale (min= -max) {auto}
%       'itcmax'    = [real] set the ITC image maximum for the color scale {auto}
%       'hzdir'     = ['up' or 'normal'|'down' or 'reverse'] Direction of
%                     the frequency axes {default: as in icadefs.m, or 'up'}
%       'ydir'      = ['up' or 'normal'|'down' or 'reverse'] Direction of
%                     the ERP axis plotted below the ITC {as in icadefs.m, or 'up'}
%       'erplim'    = [min max] ERP limits for ITC (below ITC image)       {auto}
%       'itcavglim' = [min max] average ITC limits for all freq. (left of ITC) {auto}
%       'speclim'   = [min max] average spectrum limits (left of ERSP image)   {auto}
%       'erspmarglim' = [min max] average marginal ERSP limits (below ERSP image) {auto}
%       'title'     = Optional figure or (brief) title {none}. For multiple conditions
%                     this must contain a cell array of 2 or 3 title strings.
%       'marktimes' = Non-0 times to mark with a dotted vertical line (ms)     {none}
%       'linewidth' = Line width for 'marktimes' traces (thick=2, thin=1)      {2}
%       'axesfont'  = Axes text font size                                      {10}
%       'titlefont' = Title text font size                                     {8}
%       'vert'      = [times_vector] -> plot vertical dashed lines at specified times
%                     in ms. {default: none}
%       'newfig'    = ['on'|'off'] Create new figure for difference plots {'on'}
%       'caption'   = Caption of the figure {none}
%       'outputformat' = ['old'|'plot'] for compatibility with script that used the 
%                        old output format, set to 'old' (mbase in absolute amplitude (not
%                        dB) and real itc instead of complex itc). 'plot' returns
%                        the plotted result {default: 'plot'}
% Outputs:
%            ersp   = (nfreqs,timesout) matrix of log spectral diffs from baseline
%                     (in dB log scale or absolute scale). Use the 'plot' output format
%                     above to output the ERSP as shown on the plot.
%            itc    = (nfreqs,timesout) matrix of complex inter-trial coherencies.
%                     itc is complex -- ITC magnitude is abs(itc); ITC phase in radians
%                     is angle(itc), or in deg phase(itc)*180/pi.
%          powbase  = baseline power spectrum. Note that even, when selecting the 
%                     the 'trialbase' option, the average power spectrum is
%                     returned (not trial based). To obtain the baseline of
%                     each trial, recompute it manually using the tfdata
%                     output described below.
%            times  = vector of output times (spectral time window centers) (in ms).
%            freqs  = vector of frequency bin centers (in Hz).
%         erspboot  = (nfreqs,2) matrix of [lower upper] ERSP significance.
%          itcboot  = (nfreqs) matrix of [upper] abs(itc) threshold.
%           tfdata  = optional (nfreqs,timesout,trials) time/frequency decomposition 
%                      of the single data trials. Values are complex.
%
% Plot description:
%   Assuming both 'plotersp' and 'plotitc' options are 'on' (= default). 
%   The upper panel presents the data ERSP (Event-Related Spectral Perturbation) 
%   in dB, with mean baseline spectral activity (in dB) subtracted. Use 
%   "'baseline', NaN" to prevent timef() from removing the baseline. 
%   The lower panel presents the data ITC (Inter-Trial Coherence). 
%   Click on any plot axes to pop up a new window (using 'axcopy()')
%   -- Upper left marginal panel presents the mean spectrum during the baseline 
%      period (blue), and when significance is set, the significance threshold 
%      at each frequency (dotted green-black trace).
%   -- The marginal panel under the ERSP image shows the maximum (green) and 
%      minimum (blue) ERSP values relative to baseline power at each frequency.
%   -- The lower left marginal panel shows mean ITC across the imaged time range 
%      (blue), and when significance is set, the significance threshold (dotted 
%      green-black).  
%   -- The marginal panel under the ITC image shows the ERP (which is produced by 
%      ITC across the data spectral pass band).
%
% Authors: Arnaud Delorme, Jean Hausser from timef() by Sigurd Enghoff, Scott Makeig
%          CNL / Salk Institute 1998- | SCCN/INC, UCSD 2002-
%
% See also: timefreq(), condstat(), newcrossf(), tftopo()

%    Deprecated Multitaper Parameters: [not included here]
%       'mtaper'    = If [N W], performs multitaper decomposition.
%                      (N is the time resolution and W the frequency resolution;
%                      maximum taper number is 2NW-1). Overwrites 'winsize' and 'padratio'.
%                     If [N W K], forces the use of K Slepian tapers (if possible).
%                      Phase is calculated using standard methods.
%                      The use of mutitaper with wavelets (cycles>0) is not
%                      recommended (as multiwavelets are not implemented).
%                      Uses Matlab functions DPSS, PMTM.      {no multitaper}

%    Deprecated time warp keywords (working?)
%      'timewarpfr' = {{[events], [warpfr], [plotidx]}} Time warp amplitude and phase
%                     time-courses (after time/freq transform but before smoothingtimefreqfunc
%                     across trials). 'events' is a matrix whose columns specify the
%                     epoch frames [1 ... end] at which a series of successive events
%                     occur in each trial. 'warpfr' is an optional vector of event
%                     frames to which the series of events should be time locked.
%                     (Note: Epoch start and end should not be declared as events or
%                     warpfr}. If 'warpfr' is absent or [], the median of each 'events'
%                     column will be used. [plotidx] is an optional vector of indices
%                     telling which of the warpfr to plot with vertical lines. If
%                     undefined, all marks are plotted. Overwrites 'vert' argument,
%                     if any. [Note: In future releases, 'timewarpfr' will be deprecated
%                     in favor of 'timewarp' using latencies in ms instead of frames].

%    Deprecated original time warp keywords (working?)
%       'timeStretchMarks' = [(marks,trials) matrix] Each trial data will be
%                     linearly warped (after time/freq. transform) so that the
%                     event marks are time locked to the reference frames
%                     (see timeStretchRefs). Marks must be specified in frames
%       'timeStretchRefs' = [1 x marks] Common reference frames to all trials.
%                     If empty or undefined, median latency for each mark will be used.boottype
%       'timeStretchPlot' = [vector] Indicates the indices of the reference frames
%                     (in StretchRefs) should be overplotted on the ERSP and ITC.
%
%
% Copyright (C) University of California San Diego, La Jolla, CA
%
% First built as timef.m at CNL / Salk Institute 8/1/98-8/28/01 by
% Sigurd Enghoff and Scott Makeig, edited by Arnaud Delorme
% SCCN/INC/UCSD/ reprogrammed as newtimef -Arnaud Delorme 2002-
% SCCN/INC/UCSD/ added time warping capabilities -Jean Hausser 2005
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

% 10-19-98 avoided division by zero (using MIN_ABS) -sm
% 10-19-98 improved usage message and commandline info printing -sm
% 10-19-98 made valid [] values for tvec and g.elocs -sm
% 04-01-99 added missing freq in freqs and plots, fixed log scaling bug -se && -tpj
% 06-29-99 fixed frequency indexing for constant-Q -se
% 08-24-99 reworked to handle NaN input values -sm
% 12-07-99 adjusted ERPtimes to plot ERP under ITC -sm
% 12-22-99 debugged ERPtimes, added BASE_BOOT -sm
% 01-10-00 debugged BASE_BOOT=0 -sm
% 02-28-00 added NOTE on formula derivation below -sm
% 03-16-00 added axcopy() feature -sm && tpj
% 04-16-00 added multiple marktimes loop -sm
% 04-20-00 fixed ITC cbar limits when spcified in input -sm
% 07-29-00 changed frequencies displayed msg -sm
% 10-12-00 fixed bug in freqs when cycles>0 -sm
% 02-07-01 fixed inconsistency in BASE_BOOT use -sm
% 08-28-01 matlab 'key' value arguments -ad
% 08-28-01 multitaper decomposition -ad
% 01-25-02 reformated help && license -ad
% 03-08-02 debug && compare to old timef function -ad
% 03-16-02 timeout automatically adjusted if too high -ad
% 04-02-02 added 'coher' option -ad

function [P,R,mbase,timesout,freqs,Pboot,Rboot,alltfX,PA] = newtimef( data, frames, tlimits, Fs, varwin, varargin);

% Note: Above, PA is output of 'phsamp','on'

% For future 'timewarp' keyword help: 'timewarp' 3rd element {colors} contains a
%               list of Matlab linestyles to use for vertical lines marking the occurrence
%               of the time warped events. If '', no line will be drawn for this event
%               column. If fewer colors than event columns, cycles through the given color
%               labels.  Note: Not compatible with 'vert' (below).

%varwin,winsize,g.timesout,g.padratio,g.maxfreq,g.topovec,g.elocs,g.alpha,g.marktimes,g.powbase,g.pboot,g.rboot)

% ITC:   Normally, R = |Sum(Pxy)| / (Sum(|Pxx|)*Sum(|Pyy|)) is coherence.
%        But here, we consider    Phase(Pyy) = 0 and |Pyy| = 1 -> Pxy = Pxx
%        Giving, R = |Sum(Pxx)|/Sum(|Pxx|), the inter-trial coherence (ITC)
%        Also called 'phase-locking factor' by Tallon-Baudry et al. (1996)

if nargin < 1
    help newtimef;
    return;
end

% Read system (or directory) constants and preferences:
% ------------------------------------------------------
icadefs % read local EEGLAB constants: HZDIR, YDIR, DEFAULT_SRATE, DEFAULT_TIMLIM

if ~exist('HZDIR'), HZDIR = 'up'; end; % ascending freqs
if ~exist('YDIR'), YDIR = 'up'; end;   % positive up

if YDIR == 1, YDIR = 'up'; end;        % convert from [-1|1] as set in icadefs.m  
if YDIR == -1, YDIR = 'down'; end;     % and read by other plotting functions

if ~exist('DEFAULT_SRATE'), DEFAULT_SRATE = 250; end;            % 250 Hz
if ~exist('DEFAULT_TIMLIM'), DEFAULT_TIMLIM = [-1000 2000]; end; % [-1 2] s epochs

if ~exist('DEFAULT_COLORMAP'), DEFAULT_COLORMAP = 'jet(256)'; end; % Default colormap

% Constants set here:
% ------------------
ERSP_CAXIS_LIMIT = 0;           % 0 -> use data limits; else positive value
% giving symmetric +/- caxis limits.
ITC_CAXIS_LIMIT  = 0;           % 0 -> use data limits; else positive value
% giving symmetric +/- caxis limits.
MIN_ABS          = 1e-8;        % avoid division by ~zero

% Command line argument defaults:
% ------------------------------
DEFAULT_NWIN	= 200;		% Number of windows = horizontal resolution
DEFAULT_VARWIN	= 0;		% Fixed window length or fixed number of cycles.
% =0: fix window length to that determined by nwin
% >0: set window length equal to varwin cycles
%     Bounded above by winsize, which determines
%     the min. freq. to be computed.

DEFAULT_OVERSMP	= 2;		% Number of times to oversample frequencies
DEFAULT_MAXFREQ = 50;		% Maximum frequency to display (Hz)
DEFAULT_TITLE	= '';		% Figure title (no default)
DEFAULT_ELOC    = 'chan.locs';	% Channel location file
DEFAULT_ALPHA   = NaN;		% Percentile of bins to keep
DEFAULT_MARKTIME= NaN;

% Font sizes:
AXES_FONT       = 10;           % axes text FontSize
TITLE_FONT      =  8;

if (nargin < 2)
    frames = floor((DEFAULT_TIMLIN(2)-DEFAULT_TIMLIM(1))/DEFAULT_SRATE);
elseif (~isnumeric(frames) || length(frames)~=1 || frames~=round(frames))
    error('Value of frames must be an integer.');
elseif (frames <= 0)
    error('Value of frames must be positive.');
end

DEFAULT_WINSIZE = max(pow2(nextpow2(frames)-3),4);
DEFAULT_PAD = max(pow2(nextpow2(DEFAULT_WINSIZE)),4);

if (nargin < 1)
    help newtimef
    return
end

if ischar(data) && strcmp(data,'details')
    more on
    help timefdetails
    more off
    return
end
if ~iscell(data)
    data = reshape_data(data, frames);
    trials = size(data,ndims(data));
else
    if ndims(data) == 3 && size(data,1) == 1
        error('Cannot process multiple channel component in compare mode');
    end
    [data{1}, frames] = reshape_data(data{1}, frames);
    [data{2}, frames] = reshape_data(data{2}, frames);
    trials = size(data{1},2);
end

if (nargin < 3)
    tlimits = DEFAULT_TIMLIM;
elseif (~isnumeric(tlimits) || sum(size(tlimits))~=3)
    error('Value of tlimits must be a vector containing two numbers.');
elseif (tlimits(1) >= tlimits(2))
    error('tlimits interval must be ascending.');
end

if (nargin < 4)
    Fs = DEFAULT_SRATE;
elseif (~isnumeric(Fs) || length(Fs)~=1)
    error('Value of srate must be a number.');
elseif (Fs <= 0)
    error('Value of srate must be positive.');
end

if (nargin < 5)
    varwin = DEFAULT_VARWIN;
elseif ~isnumeric(varwin) && strcmpi(varwin, 'cycles')
    varwin = varargin{1};
    varargin(1) = [];
elseif (varwin < 0)
    error('Value of cycles must be zero or positive.');
end

% build a structure for keyword arguments
% --------------------------------------
if ~isempty(varargin)
    [tmp indices] = unique_bc(varargin(1:2:end));
    varargin = varargin(sort(union(indices*2-1, indices*2))); % these 2 lines remove duplicate arguments
    try, g = struct(varargin{:});
    catch, error('Argument error in the {''param'', value} sequence'); end
end

[ g timefreqopts ] = finputcheck(varargin, ...
    {'boottype'      'string'    {'shuffle','rand','randall'}    'shuffle'; ...
    'condboot'      'string'    {'abs','angle','complex'}       'abs'; ...
    'title'         { 'string','cell' }   { [] [] }         DEFAULT_TITLE; ...
    'title2'        'string'    []          DEFAULT_TITLE; ...
    'winsize'       'integer'      [0 Inf]  DEFAULT_WINSIZE; ...
    'pad'           'real'      []          DEFAULT_PAD; ...
    'timesout'      'integer'   []          DEFAULT_NWIN; ...
    'padratio'      'integer'   [0 Inf]     DEFAULT_OVERSMP; ...
    'topovec'       'real'      []          []; ...
    'elocs'         {'string','struct'} []  DEFAULT_ELOC; ...
    'alpha'         'real'      [0 0.5]     DEFAULT_ALPHA; ...
    'marktimes'     'real'      []          DEFAULT_MARKTIME; ...
    'powbase'       'real'      []          NaN; ...
    'pboot'         'real'      []          NaN; ...
    'rboot'         'real'      []          NaN; ...
    'plotersp'      'string'    {'on','off'} 'on'; ...
    'plotamp'       'string'    {'on','off'} 'on'; ...
    'plotitc'       'string'    {'on','off'} 'on'; ...
    'detrend'       'string'    {'on','off'} 'off'; ...
    'rmerp'         'string'    {'on','off'} 'off'; ...
    'basenorm'      'string'    {'on','off'} 'off'; ...
    'commonbase'    'string'    {'on','off'} 'on'; ...
    'baseline'      'real'      []           0; ...
    'baseboot'      'real'      []           1; ...
    'linewidth'     'integer'   [1 2]        2; ...
    'naccu'         'integer'   [1 Inf]      200; ...
    'mtaper'        'real'      []           []; ...
    'maxfreq'       'real'      [0 Inf]      DEFAULT_MAXFREQ; ...
    'freqs'         'real'      [0 Inf]      [0 DEFAULT_MAXFREQ]; ...
    'cycles'        'integer'   []           []; ...
    'nfreqs'        'integer'   []           []; ...
    'freqscale'     'string'    []           'linear'; ...
    'vert'          'real'      []           [];  ...
    'newfig'        'string'    {'on','off'} 'on'; ...
    'type'          'string'    {'coher','phasecoher','phasecoher2'}  'phasecoher'; ...
    'itctype'       'string'    {'coher','phasecoher','phasecoher2'}  'phasecoher'; ...
    'outputformat'  'string'    {'old','new','plot' } 'plot'; ...
    'phsamp'        'string'    {'on','off'} 'off'; ...  % phsamp not completed - Toby 9.28.2006
    'plotphaseonly' 'string'    {'on','off'} 'off'; ...
    'plotphasesign' 'string'    {'on','off'} 'on'; ...
    'plotphase'     'string'    {'on','off'} 'on'; ... % same as above for backward compatibility
    'pcontour'      'string'    {'on','off'} 'off'; ... 
    'precomputed'   'struct'    []           struct([]); ...
    'itcmax'        'real'      []           []; ...
    'erspmax'       'real'      []           []; ...
    'lowmem'        'string'    {'on','off'} 'off'; ...
    'verbose'       'string'    {'on','off'} 'on'; ...
    'plottype'      'string'    {'image','curve'}   'image'; ...
    'mcorrect'      'string'    {'fdr','none'}      'none'; ...
    'plotmean'      'string'    {'on','off'} 'on'; ...
    'plotmode'      'string'    {}           ''; ... % for metaplottopo
    'highlightmode' 'string'    {'background','bottom'}     'background'; ...
    'chaninfo'      'struct'    []           struct([]); ...
    'erspmarglim'   'real'      []           []; ...
    'itcavglim'     'real'      []           []; ...
    'erplim'        'real'      []           []; ...
    'speclim'       'real'      []           []; ...
    'ntimesout'     'real'      []           []; ...
    'scale'         'string'    { 'log','abs'} 'log'; ...
    'timewarp'      'real'      []           []; ...
    'timewarpms'    'real'      []           []; ...
    'timewarpfr'    'real'      []           []; ...
    'timewarpidx'   'real'      []           []; ...
    'timewarpidx'   'real'      []           []; ...
    'timeStretchMarks'  'real'  []           []; ...
    'timeStretchRefs'   'real'  []           []; ...
    'timeStretchPlot'   'real'  []           []; ...
    'trialbase'     'string'    {'on','off','full'} 'off'; 
    'caption'       'string'    []           ''; ...
    'hzdir'         'string'    {'up','down','normal','reverse'}   HZDIR; ...
    'ydir'          'string'    {'up','down','normal','reverse'}   YDIR; ...
    'cycleinc'      'string'    {'linear','log'}        'linear'
    'colormap'      {'string' 'float' }    []            DEFAULT_COLORMAP;...
    }, 'newtimef', 'ignore');
if ischar(g), error(g); end
if strcmpi(g.plotamp, 'off'), g.plotersp = 'off'; end;    
if strcmpi(g.basenorm, 'on'), g.scale = 'abs'; end
if ~strcmpi(g.itctype , 'phasecoher'), g.type = g.itctype; end

g.tlimits = tlimits;
g.frames  = frames;
g.srate   = Fs;
if isempty(g.cycles)
    g.cycles  = varwin;
end
g.AXES_FONT        = AXES_FONT;      % axes text FontSize
g.TITLE_FONT       = TITLE_FONT;
g.ERSP_CAXIS_LIMIT = ERSP_CAXIS_LIMIT;
g.ITC_CAXIS_LIMIT  = ITC_CAXIS_LIMIT;
if ~strcmpi(g.plotphase, 'on'), g.plotphasesign = g.plotphase; end

% unpack 'timewarp' (and undocumented 'timewarpfr') arguments
%------------------------------------------------------------
if isfield(g,'timewarpfr')
    if iscell(g.timewarpfr) && length(g.timewarpfr) > 3
        error('undocumented ''timewarpfr'' cell array may have at most 3 elements');
    end
end

if ~isempty(g.nfreqs)
    verboseprintf(g.verbose, 'Warning: ''nfreqs'' input overwrite ''padratio''\n');
end
if strcmpi(g.basenorm, 'on')
    verboseprintf(g.verbose, 'Baseline normalization is on (results will be shown as z-scores)\n');
end

if isfield(g,'timewarp') && ~isempty(g.timewarp)
    if ndims(data) == 3
        error('Cannot perform time warping on 3-D data input');
    end
    if ~isempty(g.timewarp) % convert timewarp ms to timewarpfr frames -sm
        fprintf('\n')
        if iscell(g.timewarp)
           error('timewarp argument must be a (total_trials,epoch_events) matrix');
        end
        evntms = g.timewarp;
        warpfr = round((evntms - g.tlimits(1))/1000*g.srate)+1;
        g.timewarpfr{1} = warpfr';

        if isfield(g,'timewarpms')
           refms = g.timewarpms;
           reffr = round((refms - g.tlimits(1))/1000*g.srate)+1;
           g.timewarpfr{2} = reffr';
        end
        if isfield(g,'timewarpidx')
           g.timewarpfr{3} = g.timewarpidx;
        end
    end

    % convert again to timeStretch parameters
    % ---------------------------------------
    if ~isempty(g.timewarpfr)
        g.timeStretchMarks = g.timewarpfr{1};
        if length(g.timewarpfr) > 1
            g.timeStretchRefs = g.timewarpfr{2};
        end

        if length(g.timewarpfr) > 2
          if isempty(g.timewarpfr{3})
            stretchevents = size(g.timeStretchMarks,1);
            g.timeStretchPlot = [1:stretchevents]; % default to plotting all lines
          else
            g.timeStretchPlot = g.timewarpfr{3};
          end
        end

        if max(max(g.timeStretchMarks)) > frames-2 || min(min(g.timeStretchMarks)) < 3
            error('Time warping events must be inside the epochs.');
        end
        if ~isempty(g.timeStretchRefs)
            if max(g.timeStretchRefs) > frames-2 || min(g.timeStretchRefs) < 3
                error('Time warping reference latencies must be within the epochs.');
            end
        end
    end
end

% Determining source of the call 
% --------------------------------------% 'guicall'= 1 if newtimef is called 
callerstr = dbstack(1);                 % from EEGLAB GUI, otherwise 'guicall'= 0
if isempty(callerstr)                   % 7/3/2014, Ramon
    guicall = 0;
elseif strcmp(callerstr(end).name,'pop_newtimef')     
    guicall = 1;
else
    guicall = 0;
end

% test argument consistency
% --------------------------
if g.tlimits(2)-g.tlimits(1) < 30
    verboseprintf(g.verbose, 'newtimef(): WARNING: Specified time range is very small (< 30 ms)???\n');
    verboseprintf(g.verbose, '                     Epoch time limits should be in msec, not seconds!\n');
end

if (g.winsize > g.frames)
    error('Value of winsize must be smaller than epoch frames.');
end

if length(g.timesout) == 1 && g.timesout > 0
    if g.timesout > g.frames-g.winsize
        g.timesout = g.frames-g.winsize;
        disp(['Value of timesout must be <= frames-winsize, timeout adjusted to ' int2str(g.timesout) ]);
    end
end

if (pow2(nextpow2(g.padratio)) ~= g.padratio)
    error('Value of padratio must be an integer power of two [1,2,4,8,16,...]');
end

if (g.maxfreq > Fs/2)
    verboseprintf(g.verbose, ['Warning: value of maxfreq reduced to Nyquist rate' ...
        ' (%3.2f)\n\n'], Fs/2);
    g.maxfreq = Fs/2;
end
if g.maxfreq ~= DEFAULT_MAXFREQ, g.freqs(2) = g.maxfreq; end

if isempty(g.topovec)
    g.topovec = [];
    if isempty(g.elocs)
        error('Channel location file must be specified.');
    end
end

% naccu adjustment for FDR
% ------------------------
if (round(g.naccu*g.alpha) < 10)
    verboseprintf(g.verbose, 'Value of alpha is outside its normal range [%g,0.5]\n',10/g.naccu);
    g.naccu = round(10/g.alpha);
    verboseprintf(g.verbose, '  Increasing the number of iterations to %d\n',g.naccu);
end

if ~isnan(g.alpha)
    if length(g.baseboot) == 2
        verboseprintf(g.verbose, 'Permutation analysis will use data from %3.2g to %3.2g ms.\n', ...
            g.baseboot(1),  g.baseboot(2))
    elseif g.baseboot > 0
        verboseprintf(g.verbose, 'Permutation analysis will use data in (pre-0) baseline subwindows only.\n')
    else
        verboseprintf(g.verbose, 'Permutation analysis will use data in all subwindows.\n')
    end
end

if ~isempty(g.timeStretchMarks) % timeStretch code by Jean Hauser
    if isempty(g.timeStretchRefs)
        verboseprintf(g.verbose, ['Using median event latencies as reference event times for time warping.\n']);
        g.timeStretchRefs = median(g.timeStretchMarks,2); 
                                          % Note: Uses (grand) median latencies for two conditions
    else
        verboseprintf(g.verbose, ['Using supplied latencies as reference event times for time warping.\n']);
    end
    if isempty(g.timeStretchPlot)
        verboseprintf(g.verbose, 'Will not overplot the reference event times on the ERSP.\n');
    elseif length(g.timeStretchPlot) > 0
        g.vert = ((g.timeStretchRefs(g.timeStretchPlot)-1) ...
            /g.srate+g.tlimits(1)/1000)*1000;
        fprintf('Plotting timewarp markers at ')
           for li = 1:length(g.vert), fprintf('%d ',g.vert(li)); end
        fprintf(' ms.\n')
    end
end 

if ~isempty(g.vert)
    if min(g.vert(:)) < g.tlimits(1) || max(g.vert(:)) > g.tlimits(2)
        error('vertical line (''vert'') latency outside of epoch boundaries');
    end
end

if strcmp(g.hzdir,'up') || strcmp(g.hzdir,'normal')
    g.hzdir = 'normal'; % convert to Matlab graphics constants
elseif strcmp(g.hzdir,'down') || strcmp(g.hzdir,'reverse') || g.hzdir==-1
    g.hzdir = 'reverse';
else
    error('unknown ''hzdir'' argument'); 
end

if strcmp(g.ydir,'up') || strcmp(g.ydir,'normal')
    g.ydir = 'normal'; % convert to Matlab graphics constants
elseif strcmp(g.ydir,'down') || strcmp(g.ydir,'reverse')
    g.ydir = 'reverse';
else
    error('unknown ''ydir'' argument'); 
end

% -----------------
% ERSP scaling unit
% -----------------
if strcmpi(g.scale, 'log')
    if strcmpi(g.basenorm, 'on')
        g.unitpower = '10*log(std.)'; % impossible
    elseif isnan(g.baseline)
        g.unitpower = '10*log10(\muV^{2}/Hz)';
    else
        g.unitpower = 'dB';
    end
else
    if strcmpi(g.basenorm, 'on')
        g.unitpower = 'std.';
    elseif isnan(g.baseline)
        g.unitpower = '\muV^{2}/Hz';
    else
        g.unitpower = '% of baseline';
    end
end

% Multitaper - used in timef
% --------------------------
if ~isempty(g.mtaper) % multitaper, inspired from a Bijan Pesaran matlab function
    if length(g.mtaper) < 3
        %error('mtaper argument must be [N W] or [N W K]');

        if g.mtaper(1) * g.mtaper(2) < 1
            error('mtaper 2 first arguments'' product must be larger than 1');
        end
        if length(g.mtaper) == 2
            g.mtaper(3) = floor( 2*g.mtaper(2)*g.mtaper(1) - 1);
        end
        if length(g.mtaper) == 3
            if g.mtaper(3) > 2 * g.mtaper(1) * g.mtaper(2) -1
                error('mtaper number too high (maximum (2*N*W-1))');
            end
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
    end
    g.winsize = size(g.alltapers, 1);
    g.pad = max(pow2(nextpow2(g.winsize)),256); % pad*nextpow
    nfk = floor([0 g.maxfreq]./g.srate.*g.pad);
    g.padratio = 2*nfk(2)/g.winsize;

    %compute number of frequencies
    %nf = max(256, g.pad*2^nextpow2(g.winsize+1));
    %nfk = floor([0 g.maxfreq]./g.srate.*nf);

    %freqs = linspace( 0, g.maxfreq, diff(nfk)); % this also works in the case of a FFT
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute frequency by frequency if low memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(g.lowmem, 'on') && numel(data) ~= g.frames && isempty(g.nfreqs) && ~iscell(data)
    disp('Lowmem is a deprecated option that is not functional any more');
    return;
    
    % NOTE: the code below is functional but the graphical output is
    % different when the 'lowmem' option is used compared to when it is not
    % used - AD, 29 April 2011
    
    % compute for first 2 trials to get freqsout
    XX = reshape(data, 1, frames, prod(size(data))/g.frames);
    [P,R,mbase,timesout,freqsout] = newtimef(XX(1,:,1), frames, tlimits, Fs, g.cycles, 'plotitc', 'off', 'plotamp', 'off',varargin{:}, 'lowmem', 'off');

    % scan all frequencies
    for index = 1:length(freqsout)
        if nargout < 8
            [P(index,:),R(index,:),mbase(index),timesout,tmpfreqs(index),Pboottmp,Rboottmp] = ...
                newtimef(data, frames, tlimits, Fs, g.cycles, ...
                          'freqs', [freqsout(index) freqsout(index)], 'nfreqs', 1, ...
                             'plotamp', 'off', 'plotitc', 'off', 'plotphasesign', 'off',varargin{:}, ...
                                  'lowmem', 'off', 'timesout', timesout);
            if ~isempty(Pboottmp)
                Pboot(index,:) = Pboottmp;
                Rboot(index,:) = Rboottmp;
            else
                Pboot = [];
                Rboot = [];
            end
        else
            [P(index,:),R(index,:),mbase(index),timesout,tmpfreqs(index),Pboot(index,:),Rboot(index,:), ...
                alltfX(index,:,:)] = ...
                newtimef(data, frames, tlimits, Fs, g.cycles, ...
                            'freqs', [freqsout(index) freqsout(index)], 'nfreqs', 1, ...
                                  'plotamp', 'off', 'plotphasesign', 'off',varargin{:}, ...
                                          'lowmem', 'off', 'timesout', timesout);
        end
    end

    % compute trial-average ERP 
    % -------------------------
    ERP = mean(data,2);

    % plot results 
    %-------------
    plottimef(P, R, Pboot, Rboot, ERP, freqsout, timesout, mbase, [], [], g);

    return; % finished
end


%%%%%%%%%%%%%%%%%%%%%%%
% compare 2 conditions 
%%%%%%%%%%%%%%%%%%%%%%%
if iscell(data)
    if ~guicall && (strcmp(g.basenorm, 'on') || strcmp(g.trialbase, 'on'))  % ------------------------------------- Temporary fix for error when using
        error('EEGLAB error: basenorm and/or trialbase options cannot be used when processing 2 conditions');     % basenorm or trialbase with two conditions
    end
    Pboot = [];
    Rboot = [];
    if ~strcmpi(g.mcorrect, 'none')
        error('Correction for multiple comparison not implemented for comparing conditions');
    end
    
    vararginori = varargin;
    if length(data) ~= 2
        error('newtimef: to compare two conditions, data must be a length-2 cell array');
    end
    
    % deal with titles
    % ----------------
    for index = 1:2:length(vararginori)
        if index<=length(vararginori) % needed if elements are deleted
            
            %  if      strcmp(vararginori{index}, 'title') | ... % Added by Jean Hauser
            %          strcmp(vararginori{index}, 'title2') | ...
            if strcmp(vararginori{index}, 'timeStretchMarks') || ...
                    strcmp(vararginori{index}, 'timeStretchRefs') || ...
                    strcmp(vararginori{index}, 'timeStretchPlots')
                vararginori(index:index+1) = [];
            end
        end
    end
    if iscell(g.title) && length(g.title) >= 2 % Changed that part because providing titles
        % as cells caused the function to crash (why?)
        % at line 704 (g.tlimits = tlimits) -Jean
        if length(g.title) == 2,
            g.title{3} = [ g.title{1} ' - '  g.title{2} ];
        end
    else
        disp('Warning: title must be a cell array');
        g.title = { 'Condition 1' 'Condition 2' 'Condition 1 minus Condition 2' };
    end
    
    verboseprintf(g.verbose, '\nRunning newtimef() on Condition 1 **********************\n\n');
    
    verboseprintf(g.verbose, 'Note: If an out-of-memory error occurs, try reducing the\n');
    verboseprintf(g.verbose, '      the number of time points or number of frequencies\n');
    verboseprintf(g.verbose, '(''coher'' options take 3 times the memory of other options)\n\n');
    
    cond_1_epochs = size(data{1},2);
    
    if ~isempty(g.timeStretchMarks)
        [P1,R1,mbase1,timesout,freqs,Pboot1,Rboot1,alltfX1] = ...
            newtimef( data{1}, frames, tlimits, Fs, g.cycles, 'plotitc', 'off', ...
            'plotersp', 'off', vararginori{:}, 'lowmem', 'off', ...
            'timeStretchMarks', g.timeStretchMarks(:,1:cond_1_epochs), ...
            'timeStretchRefs', g.timeStretchRefs);
    else
        [P1,R1,mbase1,timesout,freqs,Pboot1,Rboot1,alltfX1] = ...
            newtimef( data{1}, frames, tlimits, Fs, g.cycles, 'plotitc', 'off', ...
            'plotersp', 'off', vararginori{:}, 'lowmem', 'off');
    end
    
    verboseprintf(g.verbose,'\nRunning newtimef() on Condition 2 **********************\n\n');
    
    [P2,R2,mbase2,timesout,freqs,Pboot2,Rboot2,alltfX2] = ...
        newtimef( data{2}, frames, tlimits, Fs, g.cycles, 'plotitc', 'off', ...
        'plotersp', 'off', vararginori{:}, 'lowmem', 'off', ...
        'timeStretchMarks', g.timeStretchMarks(:,cond_1_epochs+1:end), ...
        'timeStretchRefs', g.timeStretchRefs);
    
    verboseprintf(g.verbose,'\nComputing difference **********************\n\n');
    
    % recompute power baselines
    % -------------------------
    if ~isnan( g.baseline(1) ) && ~isnan( mbase1(1) ) && isnan(g.powbase(1)) && strcmpi(g.commonbase, 'on')
        disp('Recomputing baseline power: using the grand mean of both conditions ...');
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
        end
        verboseprintf(g.verbose, '\nSubtracting the common power baseline ...\n');
        meanmbase = mbase;
        mbase = { mbase mbase };
    elseif strcmpi(g.commonbase, 'on')
        mbase = { NaN NaN };
        meanmbase = mbase{1}; %Ramon :for bug 1657 
    else
        meanmbase = (mbase1 + mbase2)/2;
        mbase = { mbase1 mbase2 };
    end
    
    % plotting
    % --------
    if strcmpi(g.plotersp, 'on') || strcmpi(g.plotitc, 'on')
        g.titleall = g.title;
        if strcmpi(g.newfig, 'on'), figure; end; % declare a new figure
        
        % using same color scale
        % ----------------------
        if ~isfield(g, 'erspmax')
            g.erspmax = max( max(max(abs(Pboot1))), max(max(abs(Pboot2))) );
        end
        if ~isfield(g, 'itcmax')
            g.itcmax  = max( max(max(abs(Rboot1))), max(max(abs(Rboot2))) );
        end
        
        subplot(1,3,1); % plot Condition 1
        g.title = g.titleall{1};
        g = plottimef(P1, R1, Pboot1, Rboot1, mean(data{1},2), freqs, timesout, mbase{1}, [], [], g);
        g.itcavglim = [];
        
        subplot(1,3,2); % plot Condition 2
        g.title = g.titleall{2};
        plottimef(P2, R2, Pboot2, Rboot2, mean(data{2},2), freqs, timesout, mbase{2}, [], [], g);
        
        subplot(1,3,3); % plot Condition 1 - Condition 2
        g.title =  g.titleall{3};
    end
    
    if isnan(g.alpha)
        switch(g.condboot)
            case 'abs',  Rdiff = abs(R1)-abs(R2);
            case 'angle',  Rdiff = angle(R1)-angle(R2);
            case 'complex',  Rdiff = R1-R2;
        end
        if strcmpi(g.plotersp, 'on') || strcmpi(g.plotitc, 'on')
            g.erspmax = []; g.itcmax  = []; % auto scale inserted for diff
            plottimef(P1-P2, Rdiff, [], [], mean(data{1},2)-mean(data{2},2), freqs, timesout, meanmbase, [], [], g);
        end
    else
        % preprocess data and run compstat() function
        % -------------------------------------------
        alltfX1power = alltfX1.*conj(alltfX1);
        alltfX2power = alltfX2.*conj(alltfX2);
        
        if ~isnan(mbase{1}(1))
            mbase1 = 10.^(mbase{1}(1:size(alltfX1,1))'/20);
            mbase2 = 10.^(mbase{2}(1:size(alltfX1,1))'/20);
            alltfX1 = alltfX1./repmat(mbase1/2,[1 size(alltfX1,2) size(alltfX1,3)]);
            alltfX2 = alltfX2./repmat(mbase2/2,[1 size(alltfX2,2) size(alltfX2,3)]);
            alltfX1power = alltfX1power./repmat(mbase1,[1 size(alltfX1power,2) size(alltfX1power,3)]);
            alltfX2power = alltfX2power./repmat(mbase2,[1 size(alltfX2power,2) size(alltfX2power,3)]);
        end
        
        %formula = {'log10(mean(arg1,3))'};              % toby 10.02.2006
        %formula = {'log10(mean(arg1(:,:,data),3))'};
        
        formula = {'log10(mean(arg1(:,:,X),3))'};
        switch g.type
            case 'coher', % take the square of alltfx and alltfy first to speed up
                formula = { formula{1} ['sum(arg2(:,:,data),3)./sqrt(sum(arg1(:,:,data),3)*length(data) )'] };
                if strcmpi(g.lowmem, 'on')
                    for ind = 1:2:size(alltfX1power,1)
                        if ind == size(alltfX1,1), indarr = ind; else indarr = [ind:ind+1]; end
                        [resdifftmp resimagestmp res1tmp res2tmp] = ...
                            condstat(formula, g.naccu, g.alpha, {'both' 'upper'}, { '' g.condboot}, ...
                            { alltfX1power(indarr,:,:) alltfX2power(indarr,:,:) }, {alltfX1(indarr,:,:) alltfX2(indarr,:,:)});
                        resdiff{1}(indarr,:)     = resdifftmp{1};   resdiff{2}(indarr,:)     = resdifftmp{2};
                        resimages{1}(indarr,:,:) = resimagestmp{1}; resimages{2}(indarr,:,:) = resimagestmp{2};
                        res1{1}(indarr,:)        = res1tmp{1};      res1{2}(indarr,:)        = res1tmp{2};
                        res2{1}(indarr,:)        = res2tmp{1};      res2{2}(indarr,:)        = res2tmp{2};
                    end
                else
                    alltfXpower = { alltfX1power alltfX2power };
                    alltfX      = { alltfX1 alltfX2 };
                    alltfXabs   = { alltfX1abs alltfX2abs };
                    [resdiff resimages res1 res2] = condstat(formula, g.naccu, g.alpha, {'both' 'upper'}, { '' g.condboot}, alltfXpower, alltfX, alltfXabs);
                end
            case 'phasecoher2', % normalize first to speed up
                
                %formula = { formula{1} ['sum(arg2(:,:,data),3)./sum(arg3(:,:,data),3)'] };
                % toby 10/3/2006
                
                formula = { formula{1} ['sum(arg2(:,:,X),3)./sum(arg3(:,:,X),3)'] };
                alltfX1abs = sqrt(alltfX1power); % these 2 lines can be suppressed
                alltfX2abs = sqrt(alltfX2power); % by inserting sqrt(arg1(:,:,data)) instead of arg3(:,:,data))
                if strcmpi(g.lowmem, 'on')
                    for ind = 1:2:size(alltfX1abs,1)
                        if ind == size(alltfX1,1), indarr = ind; else indarr = [ind:ind+1]; end
                        [resdifftmp resimagestmp res1tmp res2tmp] = ...
                            condstat(formula, g.naccu, g.alpha, {'both' 'upper'}, { '' g.condboot}, ...
                            { alltfX1power(indarr,:,:) alltfX2power(indarr,:,:) }, {alltfX1(indarr,:,:) ...
                            alltfX2(indarr,:,:)}, { alltfX1abs(indarr,:,:) alltfX2abs(indarr,:,:) });
                        resdiff{1}(indarr,:)     = resdifftmp{1};   resdiff{2}(indarr,:)     = resdifftmp{2};
                        resimages{1}(indarr,:,:) = resimagestmp{1}; resimages{2}(indarr,:,:) = resimagestmp{2};
                        res1{1}(indarr,:)        = res1tmp{1};      res1{2}(indarr,:)        = res1tmp{2};
                        res2{1}(indarr,:)        = res2tmp{1};      res2{2}(indarr,:)        = res2tmp{2};
                    end
                else
                    alltfXpower = { alltfX1power alltfX2power };
                    alltfX      = { alltfX1 alltfX2 };
                    alltfXabs   = { alltfX1abs alltfX2abs };
                    [resdiff resimages res1 res2] = condstat(formula, g.naccu, g.alpha, {'both' 'upper'}, { '' g.condboot}, alltfXpower, alltfX, alltfXabs);
                end
            case 'phasecoher',
                
                %formula = { formula{1} ['mean(arg2,3)'] };              % toby 10.02.2006
                %formula = { formula{1} ['mean(arg2(:,:,data),3)'] };
                
                formula = { formula{1} ['mean(arg2(:,:,X),3)'] };
                if strcmpi(g.lowmem, 'on')
                    for ind = 1:2:size(alltfX1,1)
                        if ind == size(alltfX1,1), indarr = ind; else indarr = [ind:ind+1]; end
                        alltfX1norm = alltfX1(indarr,:,:)./sqrt(alltfX1(indarr,:,:).*conj(alltfX1(indarr,:,:)));
                        alltfX2norm = alltfX2(indarr,:,:)./sqrt(alltfX2(indarr,:,:).*conj(alltfX2(indarr,:,:)));
                        alltfXpower = { alltfX1power(indarr,:,:) alltfX2power(indarr,:,:) };
                        alltfXnorm  = { alltfX1norm alltfX2norm };
                        [resdifftmp resimagestmp res1tmp res2tmp] = ...
                            condstat(formula, g.naccu, g.alpha, {'both' 'both'}, { '' g.condboot}, ...
                            alltfXpower, alltfXnorm);
                        resdiff{1}(indarr,:)     = resdifftmp{1};   resdiff{2}(indarr,:)     = resdifftmp{2};
                        resimages{1}(indarr,:,:) = resimagestmp{1}; resimages{2}(indarr,:,:) = resimagestmp{2};
                        res1{1}(indarr,:)        = res1tmp{1};      res1{2}(indarr,:)        = res1tmp{2};
                        res2{1}(indarr,:)        = res2tmp{1};      res2{2}(indarr,:)        = res2tmp{2};
                    end
                else
                    alltfX1norm = alltfX1./sqrt(alltfX1.*conj(alltfX1));
                    alltfX2norm = alltfX2./sqrt(alltfX2.*conj(alltfX2)); % maybe have to suppress preprocessing -> lot of memory
                    alltfXpower = { alltfX1power alltfX2power };
                    alltfXnorm  = { alltfX1norm alltfX2norm };
                    [resdiff resimages res1 res2] = condstat(formula, g.naccu, g.alpha, {'both' 'both'}, { '' g.condboot}, ...
                        alltfXpower, alltfXnorm);
                end
        end
        
        % same as below: plottimef(P1-P2, R2-R1, 10*resimages{1}, resimages{2}, mean(data{1},2)-mean(data{2},2), freqs, times, mbase, g);
        if strcmpi(g.plotersp, 'on') || strcmpi(g.plotitc, 'on')
            g.erspmax = []; % auto scale
            g.itcmax  = []; % auto scale
            plottimef(10*resdiff{1}, resdiff{2}, 10*resimages{1}, resimages{2}, ...
                mean(data{1},2)-mean(data{2},2), freqs, timesout, meanmbase, [], [], g);
        end
        R1 = res1{2};
        R2 = res2{2};
        Rdiff = resdiff{2};
        Pboot = { Pboot1 Pboot2 10*resimages{1} };
        Rboot = { Rboot1 Rboot2 resimages{2} };
    end
    P = { P1 P2 P1-P2 };
    R = { R1 R2 Rdiff };
    
    if nargout >= 8, alltfX = { alltfX1 alltfX2 }; end
    
    return; % ********************************** END FOR MULTIPLE CONDITIONS
end

%%%%%%%%%%%%%%%%%%%%%%
% display text to user (computation performed only for display)
%%%%%%%%%%%%%%%%%%%%%%
verboseprintf(g.verbose, 'Computing Event-Related Spectral Perturbation (ERSP) and\n');
switch g.type
    case 'phasecoher',  verboseprintf(g.verbose, '  Inter-Trial Phase Coherence (ITC) images based on %d trials\n',trials);
    case 'phasecoher2', verboseprintf(g.verbose, '  Inter-Trial Phase Coherence 2 (ITC) images based on %d trials\n',trials);
    case 'coher',       verboseprintf(g.verbose, '  Linear Inter-Trial Coherence (ITC) images based on %d trials\n',trials);
end
verboseprintf(g.verbose, '  of %d frames sampled at %g Hz.\n',g.frames,g.srate);
verboseprintf(g.verbose, 'Each trial contains samples from %1.0f ms before to\n',g.tlimits(1));
verboseprintf(g.verbose, '  %1.0f ms after the timelocking event.\n',g.tlimits(2));
if ~isnan(g.alpha)
    verboseprintf(g.verbose, 'Only significant values (permutation statistics p<%g) will be colored;\n',g.alpha)
    verboseprintf(g.verbose, '  non-significant values will be plotted in green\n');
end
verboseprintf(g.verbose,'  Image frequency direction: %s\n',g.hzdir);

if isempty(g.precomputed)
    % -----------------------------------------
    % detrend over epochs (trials) if requested
    % -----------------------------------------
    if strcmpi(g.rmerp, 'on')
        if ndims(data) == 2
             data = data - mean(data,2)*ones(1, length(data(:))/g.frames);
        else data = data - repmat(mean(data,3), [1 1 trials]);
        end
    end

    % ----------------------------------------------------
    % compute time frequency decompositions, power and ITC
    % ----------------------------------------------------
    if length(g.timesout) > 1,   tmioutopt = { 'timesout' , g.timesout };
    elseif ~isempty(g.ntimesout) tmioutopt = { 'ntimesout', g.ntimesout };
    else                         tmioutopt = { 'ntimesout', g.timesout };
    end

    [alltfX freqs timesout R] = timefreq(data, g.srate, tmioutopt{:}, ...
        'winsize', g.winsize, 'tlimits', g.tlimits, 'detrend', g.detrend, ...
        'itctype', g.type, 'wavelet', g.cycles, 'verbose', g.verbose, ...
        'padratio', g.padratio, 'freqs', g.freqs, 'freqscale', g.freqscale, ...
        'nfreqs', g.nfreqs, 'timestretch', {g.timeStretchMarks', g.timeStretchRefs}, timefreqopts{:});
else
    alltfX   = g.precomputed.tfdata;
    timesout = g.precomputed.times;
    freqs    = g.precomputed.freqs;
    R = [];
    if ~isfield(g.precomputed, 'recompute') || strcmpi(g.precomputed.recompute, 'itc')
        switch g.itctype
            case 'coher',       R = alltfX ./ repmat(sqrt(sum(alltfX .* conj(alltfX),3) * size(alltfX,3)), [1 1 size(alltfX,3)]);
            case 'phasecoher2', R = alltfX ./ repmat(sum(sqrt(alltfX .* conj(alltfX)),3), [1 1 size(alltfX,3)]);
            case 'phasecoher',  R = alltfX ./ sqrt(alltfX .* conj(alltfX));
        end
        R = mean(R,3);
        if isfield(g.precomputed, 'recompute') && strcmpi(g.precomputed.recompute, 'itc')
            P = []; mbase = []; return;
        end
    end
end

if g.cycles(1) == 0
    alltfX = 2/0.375*alltfX/g.winsize; % TF and MC (12/11/2006): normalization, divide by g.winsize
    P  = alltfX.*conj(alltfX); % power    
    % TF and MC (12/14/2006): multiply by 2 account for negative frequencies,
    % and ounteract the reduction by a factor 0.375 that occurs as a result of 
    % cosine (Hann) tapering. Refer to Bug 446
    % Modified again 04/29/2011 due to comment in bug 1032
else 
    P  = alltfX.*conj(alltfX); % power for wavelets
end

% ----------------
% remove baseline
% ----------------
if strcmpi(g.scale, 'log') && ~any(isnan(g.powbase)), g.powbase = 10.^(g.powbase/10); end; 
P = newtimeftrialbaseln(P, timesout, 'baseline', g.baseline, 'basenorm', g.basenorm, 'trialbase', g.trialbase);
[P, baseln, mbase] = newtimefbaseln(P, timesout, 'baseline', g.baseline, 'basenorm', g.basenorm, ...
                                   'verbose', g.verbose, 'powbase', g.powbase, 'trialbase', g.trialbase, 'singletrials','on');
% ----------------
% phase amp option
% ----------------
if strcmpi(g.phsamp, 'on')
    disp( 'phsamp option is deprecated');
    %  switch g.phsamp
    %  case 'on'
    %PA = zeros(size(P,1),size(P,1),g.timesout); % NB: (freqs,freqs,times)
    % $$$ end                                             %       phs   amp
    %PA (freq x freq x time)
    %PA(:,:,j) = PA(:,:,j)  + (tmpX ./ abs(tmpX)) * ((P(:,j)))';
    % x-product: unit phase column
    % times amplitude row

    %tmpcx(1,:,:) = cumulX; % allow ./ below
    %for jj=1:g.timesout
    %    PA(:,:,jj) = PA(:,:,jj) ./ repmat(P(:,jj)', [size(P,1) 1]);
    %end
end

% ---------
% bootstrap
% --------- % this ensures that if bootstrap limits provided that no
% 'alpha' won't prevent application of the provided limits
if ~isnan(g.alpha) || ~isempty(find(~isnan(g.pboot))) || ~isempty(find(~isnan(g.rboot)))% if bootstrap analysis requested . . .
    
    % ERSP bootstrap
    % --------------
    if ~isempty(find(~isnan(g.pboot))) % if ERSP bootstrap limits provided already
        Pboot = g.pboot(:);
    else
        if size(g.baseboot,2) == 1
            if g.baseboot == 0, baselntmp = [];
            elseif ~isnan(g.baseline(1))
                baselntmp = baseln;
            else baselntmp = find(timesout <= 0); % if it is empty use whole epoch
            end
        else
            baselntmp = [];
            for index = 1:size(g.baseboot,1)
                tmptime   = find(timesout >= g.baseboot(index,1) & timesout <= g.baseboot(index,2));
                if isempty(tmptime),
                    fprintf('Warning: empty baseline interval [%3.2f %3.2f]\n', g.baseboot(index,1), g.baseboot(index,2));
                end
                baselntmp = union_bc(baselntmp, tmptime);
            end
        end
        if prod(size(g.baseboot)) > 2
            fprintf('Permutation statistics will use data in multiple selected windows.\n');
        elseif size(g.baseboot,2) == 2
            fprintf('Permutation statistics will use data in range %3.2g-%3.2g ms.\n', g.baseboot(1),  g.baseboot(2));
        elseif g.baseboot
            fprintf('   %d permutation statistics windows in baseline (times<%g).\n', length(baselntmp), g.baseboot)
        end
        
        % power significance
        % ------------------
        if strcmpi(g.boottype, 'shuffle')
            formula = 'mean(arg1,3);';
            [ Pboot Pboottrialstmp Pboottrials] = bootstat(P, formula, 'boottype', 'shuffle', ...
                'label', 'ERSP', 'bootside', 'both', 'naccu', g.naccu, ...
                'basevect', baselntmp, 'alpha', g.alpha, 'dimaccu', 2 );
            clear Pboottrialstmp;
        else
            center = 0;
            if strcmpi(g.basenorm, 'off'), center = 1; end
            
            % bootstrap signs
            Pboottmp    = P;
            Pboottrials = zeros([ size(P,1) size(P,2) g.naccu ]);
            for index = 1:g.naccu
                Pboottmp = (Pboottmp-center).*(ceil(rand(size(Pboottmp))*2-1)*2-1)+center;
                Pboottrials(:,:,index) = mean(Pboottmp,3);
            end
            Pboot = [];
        end
        if size(Pboot,2) == 1, Pboot = Pboot'; end
    end
    
    % ITC bootstrap
    % -------------
    if ~isempty(find(~isnan(g.rboot))) % if itc bootstrap provided
        Rboot = g.rboot;
    else
        if ~isempty(find(~isnan(g.pboot))) % if ERSP limits were provided (but ITC not)
            if size(g.baseboot,2) == 1
                if g.baseboot == 0, baselntmp = [];
                elseif ~isnan(g.baseline(1))
                    baselntmp = baseln;
                else baselntmp = find(timesout <= 0); % if it is empty use whole epoch
                end
            else
                baselntmp = [];
                for index = 1:size(g.baseboot,1)
                    tmptime   = find(timesout >= g.baseboot(index,1) && timesout <= g.baseboot(index,2));
                    if isempty(tmptime),
                        fprintf('Warning: empty baseline interval [%3.2f %3.2f]\n', g.baseboot(index,1), g.baseboot(index,2));
                    end
                    baselntmp = union_bc(baselntmp, tmptime);
                end
            end
            if prod(size(g.baseboot)) > 2
                fprintf('Permutation statistics will use data in multiple selected windows.\n');
            elseif size(g.baseboot,2) == 2
                fprintf('Permutation statistics will use data in range %3.2g-%3.2g ms.\n', g.baseboot(1),  g.baseboot(2));
            elseif g.baseboot
                fprintf('   %d permutation statistics windows in baseline (times<%g).\n', length(baselntmp), g.baseboot)
            end
        end;        
        % ITC significance
        % ----------------
        inputdata = alltfX;
        switch g.type
            case 'coher',       formula = [ 'sum(arg1,3)./sqrt(sum(arg1.*conj(arg1),3))/ sqrt(' int2str(size(alltfX,3)) ');' ];
            case 'phasecoher',  formula = [ 'mean(arg1,3);' ]; inputdata = alltfX./sqrt(alltfX.*conj(alltfX));
            case 'phasecoher2', formula = [ 'sum(arg1,3)./sum(sqrt(arg1.*conj(arg1)),3);' ];
        end
        if strcmpi(g.boottype, 'randall'), dimaccu = []; g.boottype = 'rand';
        else										 dimaccu = 2;
        end
        [Rboot Rboottmp Rboottrials] = bootstat(inputdata, formula, 'boottype', g.boottype, ...
            'label', 'ITC', 'bootside', 'upper', 'naccu', g.naccu, ...
            'basevect', baselntmp, 'alpha', g.alpha, 'dimaccu', 2 );
        fprintf('\n');
        clear Rboottmp;        
    end
else
    Pboot = []; Rboot = [];
end

% average the power
% -----------------
PA = P;
if ndims(P) == 4,     P = mean(P, 4);
elseif ndims(P) == 3, P = mean(P, 3);
end

% correction for multiple comparisons
% -----------------------------------
maskersp = [];
maskitc  = []; 
if ~isnan(g.alpha)
    if isempty(find(~isnan(g.pboot))) % if ERSP lims not provided
        if ndims(Pboottrials) < 3, Pboottrials = Pboottrials'; end
        exactp_ersp = compute_pvals(P, Pboottrials);
        if strcmpi(g.mcorrect, 'fdr')
            alphafdr = fdr(exactp_ersp, g.alpha);
            if alphafdr ~= 0
                fprintf('ERSP correction for multiple comparisons using FDR, alpha_fdr = %3.6f\n', alphafdr);
            else fprintf('ERSP correction for multiple comparisons using FDR, nothing significant\n', alphafdr);
            end
            maskersp = exactp_ersp <= alphafdr;
        else
            maskersp = exactp_ersp <= g.alpha;
        end
    end;    
    if isempty(find(~isnan(g.rboot))) % if ITC lims not provided
        exactp_itc  = compute_pvals(abs(R), abs(Rboottrials'));        
        if strcmpi(g.mcorrect, 'fdr')
            alphafdr = fdr(exactp_itc, g.alpha);
            if alphafdr ~= 0
                fprintf('ITC  correction for multiple comparisons using FDR, alpha_fdr = %3.6f\n', alphafdr);
            else fprintf('ITC  correction for multiple comparisons using FDR, nothing significant\n', alphafdr);
            end
            maskitc = exactp_itc <= alphafdr;
        else
            maskitc = exactp_itc  <= g.alpha;
        end
    end
end

% convert to log if necessary
% ---------------------------
if strcmpi(g.scale, 'log')
    if ~isnan( g.baseline(1) ) && ~isnan( mbase(1) ) && strcmpi(g.trialbase, 'off'), mbase = log10(mbase)*10; end
    P = 10 * log10(P);
    if ~isempty(Pboot)
        Pboot = 10 * log10(Pboot);
    end
end
if isempty(Pboot) && exist('maskersp')
    Pboot = maskersp;
end

% auto scalling
% -------------
if isempty(g.erspmax)
    g.erspmax = [max(max(abs(P)))]/2;
    if strcmpi(g.scale, 'abs') && strcmpi(g.basenorm, 'off') % % of baseline
        g.erspmax = [max(max(abs(P)))];
        if g.erspmax > 1
             g.erspmax = [1-(g.erspmax-1) g.erspmax];
        else g.erspmax = [g.erspmax 1+(1-g.erspmax)];
     	end
    end
    %g.erspmax = [-g.erspmax g.erspmax]+1;
end

% --------
% plotting
% --------
if strcmpi(g.plotersp, 'on') || strcmpi(g.plotitc, 'on')
    if ndims(P) == 3
        P = squeeze(P(2,:,:,:));
        R = squeeze(R(2,:,:,:));
        mbase = squeeze(mbase(2,:));
        ERP = mean(squeeze(data(1,:,:)),2);
    else      
        ERP = mean(data,2);
    end
    if strcmpi(g.plottype, 'image')
        plottimef(P, R, Pboot, Rboot, ERP, freqs, timesout, mbase, maskersp, maskitc, g);
    else
        plotallcurves(P, R, Pboot, Rboot, ERP, freqs, timesout, mbase, g);
    end
end

% --------------
% format outputs
% --------------
if strcmpi(g.outputformat, 'old')
    R = abs(R); % convert coherence vector to magnitude
    if strcmpi(g.scale, 'log'), mbase = 10.^(mbase/10); end
end
if strcmpi(g.verbose, 'on')
    disp('Note: Add output variables to command line call in history to');
    disp('      retrieve results and use the tftopo function to replot them');
end
mbase = mbase';

if ~isempty(g.caption)
    h = textsc(g.caption, 'title');
    set(h, 'FontWeight', 'bold');
end

return;

% -----------------
% plotting function
% -----------------
function g = plottimef(P, R, Pboot, Rboot, ERP, freqs, times, mbase, maskersp, maskitc, g);

persistent showwarning;

if isempty(showwarning)
    warning( [ 'Some versions of Matlab crash on this function. If this is' 10 ...
               'the case, simply comment the code line 1655-1673 in newtimef.m' 10 ...
               'which aims at "plotting marginal ERSP mean below ERSP image"' ]);
    showwarning = 1;
end;    

%
% compute ERP
%
ERPtimes = [g.tlimits(1):(g.tlimits(2)-g.tlimits(1))/(g.frames-1):g.tlimits(2)+0.000001];
ERPindices = zeros(1, length(times));
for ti=1:length(times)
    [tmp ERPindices(ti)] = min(abs(ERPtimes-times(ti)));
end
ERPtimes = ERPtimes(ERPindices); % subset of ERP frames on t/f window centers
ERP = ERP(ERPindices);

if ~isreal(R)
    Rangle = angle(R);
    Rsign = sign(imag(R));
    R = abs(R); % convert coherence vector to magnitude
    setylim = 1;
else
    Rangle = zeros(size(R)); % Ramon: if isreal(R) then we get an error because Rangle does not exist
    Rsign = ones(size(R));
    setylim = 0;
end
switch lower(g.plotitc)
    case 'on',
        switch lower(g.plotersp),
            case 'on', ordinate1 = 0.67; ordinate2 = 0.1; height = 0.33; g.plot = 1;
            case 'off', ordinate2 = 0.1; height = 0.9; g.plot = 1;
        end
    case 'off', ordinate1 = 0.1; height = 0.9;
        switch lower(g.plotersp),
            case 'on', ordinate1 = 0.1; height = 0.9;  g.plot = 1;
            case 'off', g.plot = 0;
        end
end

if g.plot
    % verboseprintf(g.verbose, '\nNow plotting...\n');
    set(gcf,'DefaultAxesFontSize',g.AXES_FONT)
    pos = get(gca,'position');
    q = [pos(1) pos(2) 0 0];
    s = [pos(3) pos(4) pos(3) pos(4)];
    axis off;
end

switch lower(g.plotersp)
    case 'on'
        %
        %%%%%%% image the ERSP %%%%%%%%%%%%%%%%%%%%%%%%%%
        %

        h(1) = axes('Position',[.1 ordinate1 .9 height].*s+q);
        set(h(1), 'tag', 'ersp');

        PP = P;
        if strcmpi(g.scale, 'abs') && strcmpi(g.basenorm, 'off')
             baseval = 1;
        else baseval = 0;
        end
        if ~isnan(g.alpha)
            if strcmpi(g.pcontour, 'off') && ~isempty(maskersp) % zero out nonsignif. power differences
                PP(~maskersp) = baseval;
                %PP = PP .* maskersp;
            elseif isempty(maskersp)
                if size(PP,1) == size(Pboot,1) && size(PP,2) == size(Pboot,2)
                    PP(find(PP > Pboot(:,:,1) & (PP < Pboot(:,:,2)))) = baseval;
                    Pboot = squeeze(mean(Pboot,2));
                    if size(Pboot,2) == 1, Pboot = Pboot'; end
                else
                    PP(find((PP > repmat(Pboot(:,1),[1 length(times)])) ...
                        & (PP < repmat(Pboot(:,2),[1 length(times)])))) = baseval;
                end
            end
        end
 
        % find color limits
        % -----------------
        if isempty(g.erspmax)
            if g.ERSP_CAXIS_LIMIT == 0
                g.erspmax = [-1 1]*1.1*max(max(abs(P(:,:))));
            else
                g.erspmax = g.ERSP_CAXIS_LIMIT*[-1 1];
            end
        elseif length(g.erspmax) == 1
            g.erspmax = [ -g.erspmax g.erspmax];
        end
        if isnan( g.baseline(1) ) && g.erspmax(1) < 0
            g.erspmax = [ min(min(P(:,:))) max(max(P(:,:)))];
        end

        % plot image
        % ----------
        if ~strcmpi(g.freqscale, 'log')
            imagesc(times,freqs,PP(:,:), g.erspmax);
        else
            imagesclogy(times,freqs,PP(:,:),g.erspmax);
        end
        set(gca,'ydir',g.hzdir);  % make frequency ascend or descend

        % put contour for multiple comparison masking
        if ~isempty(maskersp) && strcmpi(g.pcontour, 'on')
            hold on; [tmpc tmph] = contour(times, freqs, maskersp);
            set(tmph, 'linecolor', 'k', 'linewidth', 0.25)
        end
        
        hold on
        plot([0 0],[0 freqs(end)],'--m','LineWidth',g.linewidth); % plot time 0
        if ~isnan(g.marktimes) % plot marked time
            for mt = g.marktimes(:)'
                plot([mt mt],[0 freqs(end)],'--k','LineWidth',g.linewidth);
            end
        end
        hold off
        set(h(1),'YTickLabel',[],'YTick',[])
        set(h(1),'XTickLabel',[],'XTick',[])
        if ~isempty(g.vert)
            for index = 1:length(g.vert)
                line([g.vert(index), g.vert(index)], [min(freqs) max(freqs)], 'linewidth', 1, 'color', 'm');
            end
        end

        h(2) = gca;
        h(3) = cbar('vert'); % ERSP colorbar axes
        set(h(2),'Position',[.1 ordinate1 .8 height].*s+q)
        set(h(3),'Position',[.95 ordinate1 .05 height].*s+q)
        title([ 'ERSP(' g.unitpower ')' ])

        %
        %%%%% plot marginal ERSP mean below ERSP image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %

        h(4) = axes('Position',[.1 ordinate1-0.1 .8 .1].*s+q);

        E = [min(P(:,:),[],1);max(P(:,:),[],1)];

        % plotting limits
        if isempty(g.erspmarglim)
            g.erspmarglim = [min(E(1,:))-max(max(abs(E)))/3 max(E(2,:))+max(max(abs(E)))/3];
        end

        plot(times,E,[0 0],g.erspmarglim, '--m','LineWidth',g.linewidth)
        xlim([min(times) max(times)])
        ylim(g.erspmarglim)

        tick = get(h(4),'YTick');
        set(h(4),'YTick',[tick(1) ; tick(end)])
        set(h(4),'YAxisLocation','right')
        set(h(4),'TickLength',[0.020 0.025]);
        xlabel('Time (ms)')
        ylabel(g.unitpower)

        %
        %%%%% plot mean spectrum to left of ERSP image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %

        h(5) = axes('Position',[0 ordinate1 .1 height].*s+q);

        if isnan(g.baseline)              % Ramon :for bug 1657
            E = zeros(size(freqs));
        else
            E = mbase;
        end

        if ~isnan(E(1))

            % plotting limits
            if isempty(g.speclim)
               % g.speclim = [min(E)-max(abs(E))/3 max(E)+max(abs(E))/3];
               if all(~isnan(mbase))
                   g.speclim = [min(mbase)-max(abs(mbase))/3 max(mbase)+max(abs(mbase))/3]; % RMC: Just for plotting
               else
                   g.speclim = [min(E)-max(abs(E))/3 max(E)+max(abs(E))/3];
               end
            end

            % plot curves
            if ~strcmpi(g.freqscale, 'log')
                plot(freqs,E,'LineWidth',g.linewidth); hold on;
                if ~isnan(g.alpha) && size(Pboot,2) == 2
                    try
                        plot(freqs,Pboot(:,:)'+[E;E], 'g', 'LineWidth',g.linewidth)
                        plot(freqs,Pboot(:,:)'+[E;E], 'k:','LineWidth',g.linewidth)
                    catch
                        plot(freqs,Pboot(:,:)+[E E], 'g', 'LineWidth',g.linewidth)
                        plot(freqs,Pboot(:,:)+[E E], 'k:','LineWidth',g.linewidth)
                    end
                end
                if freqs(1) ~= freqs(end), xlim([freqs(1) freqs(end)]); end
                if g.speclim(1) ~= g.speclim(2), ylim(g.speclim); end; % Ramon :for bug 1657 

            else % 'log'
                semilogx(freqs,E,'LineWidth',g.linewidth); hold on;
                if ~isnan(g.alpha)
                    try
                        semilogx(freqs,Pboot(:,:)'+[E;E],'g', 'LineWidth',g.linewidth)
                        semilogx(freqs,Pboot(:,:)'+[E;E],'k:','LineWidth',g.linewidth)
                    catch
                        semilogx(freqs,Pboot(:,:)+[E E],'g', 'LineWidth',g.linewidth)
                        semilogx(freqs,Pboot(:,:)+[E E],'k:','LineWidth',g.linewidth)
                    end
                end
                if freqs(1) ~= freqs(end), xlim([freqs(1) freqs(end)]); end
                if g.speclim(1) ~= g.speclim(2), ylim(g.speclim); end; %RMC
                set(h(5),'View',[90 90])
                divs = linspace(log(freqs(1)), log(freqs(end)), 10);
                set(gca, 'xtickmode', 'manual');
                divs = ceil(exp(divs)); divs = unique_bc(divs); % ceil is critical here, round might misalign
                set(gca, 'xtick', divs);
            end
            set(h(5),'TickLength',[0.020 0.025]);
            set(h(5),'View',[90 90])
            xlabel('Frequency (Hz)')
            if strcmp(g.hzdir,'normal')
                set(gca,'xdir','reverse');
            else
                set(gca,'xdir','normal');
            end
            ylabel(g.unitpower)
            tick = get(h(5),'YTick');
            if (length(tick)>2)
                set(h(5),'YTick',[tick(1) ; tick(end-1)])
            end
        end
end

switch lower(g.plotitc)
    case 'on'
        %
        %%%%%%%%%%%% Image the ITC %%%%%%%%%%%%%%%%%%
        %
        h(6) = axes('Position',[.1 ordinate2 .9 height].*s+q); % ITC image
        if ishandle(h(1));set(h(1), 'tag', 'itc');end

        if abs(R(1,1)-1) < 0.0001, g.plotphaseonly = 'on'; end
        if strcmpi(g.plotphaseonly, 'on')
            RR = Rangle/pi*180;
        else
            RR = R;
        end
        if ~isnan(g.alpha)
            if ~isempty(maskitc) && strcmpi(g.pcontour, 'off')
                RR = RR .* maskitc;
            elseif isempty(maskitc)
                if size(RR,1) == size(Rboot,1) && size(RR,2) == size(Rboot,2)
                    tmp = gcf;
                    if size(Rboot,3) == 2	 RR(find(RR > Rboot(:,:,1) & RR < Rboot(:,:,2))) = 0;
                    else                   RR(find(RR < Rboot)) = 0;
                    end
                    Rboot = mean(Rboot(:,:,end),2);
                else
                    RR(find(RR < repmat(Rboot(:),[1 length(times)]))) = 0;
                end
            end
        end

        if g.ITC_CAXIS_LIMIT == 0
            coh_caxis = min(max(max(R(:,:))),1)*[-1 1]; % 1 WAS 0.4 !
        else
            coh_caxis = g.ITC_CAXIS_LIMIT*[-1 1];
        end

        if strcmpi(g.plotphaseonly, 'on')
            if ~strcmpi(g.freqscale, 'log')
                imagesc(times,freqs,RR(:,:)); % <---
            else
                imagesclogy(times,freqs,RR(:,:)); % <---
            end
            g.itcmax = [-180 180];
            setylim = 0;
        else
            if max(coh_caxis) == 0,              % toby 10.02.2006
                coh_caxis = [-1 1];
            end
            if ~strcmpi(g.freqscale, 'log')
                if exist('Rsign') && strcmp(g.plotphasesign, 'on')
                    imagesc(times,freqs,Rsign(:,:).*RR(:,:),coh_caxis); % <---
                else
                    imagesc(times,freqs,RR(:,:),coh_caxis); % <---
                end
            else
                if exist('Rsign') && strcmp(g.plotphasesign, 'on')
                    imagesclogy(times,freqs,Rsign(:,:).*RR(:,:),coh_caxis); % <---
                else
                    imagesclogy(times,freqs,RR(:,:),coh_caxis); % <---
                end
            end
        end
        set(gca,'ydir',g.hzdir);  % make frequency ascend or descend

        % plot contour if necessary
        if ~isempty(maskitc) && strcmpi(g.pcontour, 'on')
            hold on; [tmpc tmph] = contour(times, freqs, maskitc);
            set(tmph, 'linecolor', 'k', 'linewidth', 0.25)
        end

        if isempty(g.itcmax)
            g.itcmax = caxis;
        elseif length(g.itcmax) == 1
            g.itcmax = [ -g.itcmax g.itcmax ];
        end
        caxis(g.itcmax);

        hold on
        plot([0 0],[0 freqs(end)],'--m','LineWidth',g.linewidth);
        if ~isnan(g.marktimes)
            for mt = g.marktimes(:)'
                plot([mt mt],[0 freqs(end)],'--k','LineWidth',g.linewidth);
            end
        end
        hold off
        set(h(6),'YTickLabel',[],'YTick',[])
        set(h(6),'XTickLabel',[],'XTick',[])
        if ~isempty(g.vert)
            for index = 1:length(g.vert)
                line([g.vert(index), g.vert(index)], [min(freqs) max(freqs)], 'linewidth', 1, 'color', 'm');
            end
        end

        h(7) = gca;
        h(8) = cbar('vert');
        %h(9) = get(h(8),'Children'); % make the function crash
        set(h(7),'Position',[.1 ordinate2 .8 height].*s+q)
        set(h(8),'Position',[.95 ordinate2 .05 height].*s+q)
        if setylim
            set(h(8),'YLim',[0 g.itcmax(2)]);
        end
        if strcmpi(g.plotphaseonly, 'on')
            title('ITC phase')
        else
            title('ITC')
        end

        %
        %%%%% plot the ERP below the ITC image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %

        h(10) = axes('Position',[.1 ordinate2-0.1 .8 .1].*s+q); % ERP

        if isempty(g.erplim)
            ERPmax = max(ERP);
            ERPmin = min(ERP);
            g.erplim = [ ERPmin - 0.1*(ERPmax-ERPmin) ERPmax + 0.1*(ERPmax-ERPmin) ];
        end

        plot(ERPtimes,ERP, [0 0],g.erplim,'--m','LineWidth',g.linewidth);
        hold on;
        plot([times(1) times(length(times))],[0 0], 'k');
        xlim([min(ERPtimes) max(ERPtimes)]);
        ylim(g.erplim)
        set(gca,'ydir',g.ydir);

        tick = get(h(10),'YTick');
        set(h(10),'YTick',[tick(1) ; tick(end)])
        set(h(10),'TickLength',[0.02 0.025]);
        set(h(10),'YAxisLocation','right')
        xlabel('Time (ms)')
        ylabel('\muV')
        if (~isempty(g.topovec))
            if length(g.topovec) ~= 1, ylabel(''); end; % ICA component
        end
        E = nan_mean(R(:,:)'); % don't let a few NaN's crash this

        %
        %%%%% plot the marginal mean left of the ITC image %%%%%%%%%%%%%%%%%%%%%
        %

        h(11) = axes('Position',[0 ordinate2 .1 height].*s+q); % plot the marginal mean
        % ITC left of the ITC image
        % set plotting limits
        if isempty(g.itcavglim)
            if ~isnan(g.alpha)
                g.itcavglim = [ min(E)-max(E)/3 max(Rboot)+max(Rboot)/3];
            else
                g.itcavglim = [ min(E)-max(E)/3 max(E)+max(E)/3];
            end
        end
        if max(g.itcavglim) == 0 || any(isnan(g.itcavglim))
            g.itcavglim = [-1 1];
        end
        
        % plot marginal ITC
        if ~strcmpi(g.freqscale, 'log')
            plot(freqs,E,'LineWidth',g.linewidth); hold on;
            if ~isnan(g.alpha)
                plot(freqs,Rboot,'g', 'LineWidth',g.linewidth)
                plot(freqs,Rboot,'k:','LineWidth',g.linewidth)
            end
            if freqs(1) ~= freqs(end), xlim([freqs(1) freqs(end)]); end
            ylim(g.itcavglim)
        else
            semilogx(freqs,E,'LineWidth',g.linewidth); hold on;
            if ~isnan(g.alpha)
                semilogx(freqs,Rboot(:),'g', 'LineWidth',g.linewidth)
                semilogx(freqs,Rboot(:),'k:','LineWidth',g.linewidth)
            end
            if freqs(1) ~= freqs(end), xlim([freqs(1) freqs(end)]); end
            ylim(g.itcavglim)
            divs = linspace(log(freqs(1)), log(freqs(end)), 10);
            set(gca, 'xtickmode', 'manual');
            divs = ceil(exp(divs)); divs = unique_bc(divs); % ceil is critical here, round might misalign
            set(gca, 'xtick', divs);
         end

        % ITC plot details
        tick = get(h(11),'YTick');
        if length(tick) > 1
            set(h(11),'YTick',[tick(1) ; tick(length(tick))])
        end
        set(h(11),'View',[90 90])
        %set(h(11),'TickLength',[0.020 0.025]);
        xlabel('Frequency (Hz)')
        if strcmp(g.hzdir,'normal')
            set(gca,'xdir','reverse');
        else
            set(gca,'xdir','normal');
        end
        ylabel('ERP')

end %switch

%
%%%%%%%%%%%%%%% plot a topoplot() %%%%%%%%%%%%%%%%%%%%%%%
%
if (~isempty(g.topovec)) && strcmpi(g.plotitc, 'on') && strcmpi(g.plotersp, 'on')
    
    if strcmp(g.plotersp,'off')
        h(12) = axes('Position',[-.207 .95 .2 .14].*s+q); % place the scalp map at top-left
    else
        h(12) = axes('Position',[-.1 .43 .2 .14].*s+q);   % place the scalp map at middle-left
    end
    if length(g.topovec) == 1
        topoplot(g.topovec,g.elocs,'electrodes','off', ...
                 'style', 'blank', 'emarkersize1chan', 10, 'chaninfo', g.chaninfo);
    else
        topoplot(g.topovec,g.elocs,'electrodes','off', 'chaninfo', g.chaninfo);
    end
    axis('square')
end

if g.plot
    try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end
    if (length(g.title) > 0) && ~iscell(g.title)
        axes('Position',pos,'Visible','Off');
        h(13) = text(-.05,1.01,g.title);
        set(h(13),'VerticalAlignment','bottom')
        set(h(13),'HorizontalAlignment','left')
        set(h(13),'FontSize',g.TITLE_FONT);
    end

    try, axcopy(gcf); catch, end
end
colormap(g.colormap);

% ---------------
% Plotting curves
% ---------------
function plotallcurves(P, R, Pboot, Rboot, ERP, freqs, times, mbase, g);

if ~isreal(R)
    Rangle = angle(R);
    R = abs(R); % convert coherence vector to magnitude
    setylim = 1;
else
    Rangle = zeros(size(R)); % Ramon: if isreal(R) then we get an error because Rangle does not exist
    Rsign = ones(size(R));
    setylim = 0;
end

if strcmpi(g.plotitc, 'on') || strcmpi(g.plotersp, 'on')
    verboseprintf(g.verbose, '\nNow plotting...\n');
    pos = get(gca,'position');
    q = [pos(1) pos(2) 0 0];
    s = [pos(3) pos(4) pos(3) pos(4)];
end

% time unit
% ---------
if times(end) > 10000
    times = times/1000;
    timeunit = 's';
else
    timeunit = 'ms';
end

if strcmpi(g.plotersp, 'on')
    %
    %%%%%%% image the ERSP %%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if strcmpi(g.plotitc, 'on'), subplot(2,1,1); end
    set(gca, 'tag', 'ersp');
    alllegend = {};

    for index = 1:length(freqs)
        alllegend{index} = [ num2str(freqs(index)) 'Hz baseline ' num2str(mbase(index)) ' dB' ];
    end
    if strcmpi(g.plotmean, 'on') && freqs(1) ~= freqs(end)
        alllegend = { alllegend{:} [ num2str(freqs(1)) '-' num2str(freqs(end)) ...
            'Hz mean baseline ' num2str(mean(mbase)) ' dB' ] };
    end
    plotcurve(times, P, 'maskarray', Pboot, 'title', 'ERSP', ...
        'xlabel', [ 'Time (' timeunit ')' ], 'ylabel', 'dB', 'ylim', [-g.erspmax g.erspmax], ...
        'vert', g.vert, 'marktimes', g.marktimes, 'legend', alllegend, ...
        'linewidth', g.linewidth, 'highlightmode', g.highlightmode, 'plotmean', g.plotmean);
end

if strcmpi(g.plotitc, 'on')
    %
    %%%%%%%%%%%% Image the ITC %%%%%%%%%%%%%%%%%%
    %
    if strcmpi(g.plotersp, 'on'), subplot(2,1,2); end
    set(gca, 'tag', 'itc');
    if abs(R(1,1)-1) < 0.0001, g.plotphaseonly = 'on'; end
    if strcmpi(g.plotphaseonly, 'on') % plot ITC phase instead of amplitude (e.g. for continuous data)
        RR = Rangle/pi*180;
    else RR = R;
    end

    % find regions of significance
    % ----------------------------
    alllegend = {};
    for index = 1:length(freqs)
        alllegend{index} = [ num2str(freqs(index)) 'Hz baseline ' num2str(mbase(index)) ' dB' ];
    end
    if strcmpi(g.plotmean, 'on') && freqs(1) ~= freqs(end)
        alllegend = { alllegend{:} [ num2str(freqs(1)) '-' num2str(freqs(end)) ...
            'Hz mean baseline ' num2str(mean(mbase)) ' dB' ] };
    end
    plotcurve(times, RR, 'maskarray', Rboot, 'val2mask', R, 'title', 'ITC', ...
        'xlabel', [ 'Time (' timeunit ')' ], 'ylabel', 'dB', 'ylim', g.itcmax, ...
        'vert', g.vert, 'marktimes', g.marktimes, 'legend', alllegend, ...
        'linewidth', g.linewidth, 'highlightmode', g.highlightmode, 'plotmean', g.plotmean);
end

if strcmpi(g.plotitc, 'on') || strcmpi(g.plotersp, 'on')
    %
    %%%%%%%%%%%%%%% plot a topoplot() %%%%%%%%%%%%%%%%%%%%%%%
    %
    if (~isempty(g.topovec))
        h(12) = axes('Position',[-.1 .43 .2 .14].*s+q);
        if length(g.topovec) == 1
            topoplot(g.topovec,g.elocs,'electrodes','off', ...
                'style', 'blank', 'emarkersize1chan', 10);
        else
            topoplot(g.topovec,g.elocs,'electrodes','off');
        end
        axis('square')
    end

    try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end
    if (length(g.title) > 0) && ~iscell(g.title)
        axes('Position',pos,'Visible','Off');
        h(13) = text(-.05,1.01,g.title);
        set(h(13),'VerticalAlignment','bottom')
        set(h(13),'HorizontalAlignment','left')
        set(h(13),'FontSize',g.TITLE_FONT);
    end

    try, axcopy(gcf); catch, end
end

%
%%%%%%%%%%%%%%%%%%%%%%% Highlight regions %%%%%%%%%%%%%%%%%%%%%%%%%%
%
function highlight(ax, times, regions, highlightmode);
color1 = [0.75 0.75 0.75];
color2 = [0 0 0];
yl  = ylim;

if ~strcmpi(highlightmode, 'background')
    yl2 = [ yl(1)-(yl(2)-yl(1))*0.15   yl(1)-(yl(2)-yl(1))*0.1 ];
    tmph = patch([times(1) times(end) times(end) times(1)], ...
        [yl2(1) yl2(1) yl2(2) yl2(2)], [1 1 1]); hold on;
    ylim([ yl2(1) yl(2)]);
    set(tmph, 'edgecolor', [1 1 1]);
end

if ~isempty(regions)
    axes(ax);
    in_a_region = 0;
    for index=1:length(regions)
        if regions(index) && ~in_a_region
            tmpreg(1) = times(index);
            in_a_region = 1;
        end
        if ~regions(index) && in_a_region
            tmpreg(2) = times(index);
            in_a_region = 0;
            if strcmpi(highlightmode, 'background')
                tmph = patch([tmpreg(1) tmpreg(2) tmpreg(2) tmpreg(1)], ...
                    [yl(1) yl(1) yl(2) yl(2)], color1); hold on;
                set(tmph, 'edgecolor', color1);
            else
                tmph = patch([tmpreg(1) tmpreg(2) tmpreg(2) tmpreg(1)], ...
                    [yl2(1) yl2(1) yl2(2) yl2(2)], color2); hold on;
                set(tmph, 'edgecolor', color2);
            end
        end
    end
end

% reshaping data
% -----------
function [data, frames] = reshape_data(data, frames)
data = squeeze(data);
if min(size(data)) == 1
    if (rem(length(data),frames) ~= 0)
        error('Length of data vector must be divisible by frames.');
    end
    data = reshape(data, frames, length(data)/frames);
else
    frames = size(data,1);
end

function verboseprintf(verbose, varargin)
if strcmpi(verbose, 'on')
    fprintf(varargin{:});
end

% reshaping data
% -----------
function pvals = compute_pvals(oridat, surrog, tail)
    
    if nargin < 3
        tail = 'both';
    end
    
    if myndims(oridat) > 1        
        if size(oridat,2) ~= size(surrog, 2) || myndims(surrog) == 2
            if size(oridat,1) == size(surrog, 1)
                surrog = repmat( reshape(surrog, [size(surrog,1) 1 size(surrog,2)]), [1 size(oridat,2) 1]);
            elseif size(oridat,2) == size(surrog, 1)
                surrog = repmat( reshape(surrog, [1 size(surrog,1) size(surrog,2)]), [size(oridat,1) 1 1]);
            else
                error('Permutation statistics array size error');
            end
        end
    end

    surrog = sort(surrog, myndims(surrog)); % sort last dimension
    
    if myndims(surrog) == 1    
        surrog(end+1) = oridat;        
    elseif myndims(surrog) == 2
        surrog(:,end+1) = oridat;        
    elseif myndims(surrog) == 3
        surrog(:,:,end+1) = oridat;
    else
        surrog(:,:,:,end+1) = oridat;
    end

    [tmp idx] = sort( surrog, myndims(surrog) );
    [tmp mx]  = max( idx,[], myndims(surrog));        
                
    len = size(surrog,  myndims(surrog) );
    pvals = 1-(mx-0.5)/len;
    if strcmpi(tail, 'both')
        pvals = min(pvals, 1-pvals);
        pvals = 2*pvals;
    end;    
    
function val = myndims(a)
    if ndims(a) > 2
        val = ndims(a);
    else
        if size(a,1) == 1,
            val = 2;
        elseif size(a,2) == 1,
            val = 1;
        else
            val = 2;
        end
    end; 


