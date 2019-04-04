% timefdetails() - details of the timef() function for time/frequency analysis 
%                  of multiple epochs of single-channel event-related data.
%
% Global Description:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% timef() performs normalized time/frequency averaging using either 
% FFT-, wavelet-, or multitaper DFT estimates. The wavelets are N-cycle 
% Hanning-windowed sinusoids. (Note: To substitute for hanning() windowing 
% gauss() or other windowing, replace the timef.m reference to hanning()).
%
% By default, the two image panels of the output plot show, respectively, the
% event-related spectral perturbation (ERSP) and inter-trial coherence (ITC)
% of the input data.
%
% The ERSP: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The ERSP (S. Makeig, Electroencephalogr Clin Neurophysiol 86:283-93, 1993) 
% shows mean log event-loced deviations from epoch-mean (or baseline-mean) power
% at each frequency. If bootstrap statistics are computed (via the 'alpha'
% probability input parameter), (time, freq) points with non-significant 
% differences from 0 are colored green in the image (but not in the 'ersp' 
% output variable - use the output 'erspboot' variable to re-mask the 'ersp'
% output if desired). The baseline mean spectrum removed from each epoch
% is available as output parameter "powbase". Note that log(power) differences
% are equivalent to log(power ratios) between baseline and observed spectra,
% meaning the implicit ERSP model is one of (multiplicative) amplitude modulation
% of the EEG/MEG spectrum by e.g. subcortical and/or intra-cortical influences.
%
% In the default view, the thin bottom panel below the upper (ERSP) image shows 
% the ERSP envelope (the most positive and most negative values at each output 
% time point). The thin left panel shows the mean (or baseline) log spectrum 
% (blue trace). When bootstrap statistics are computed (via the "alpha" argument), 
% the left panel (green trace) also shows the bootstrap significant levels (+/-) 
% at each frequency.
%
% The ITC: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (Inter-trial Coherence, cf. Tallon-Baudry et al., "Phase-locking factor")
% The lower panel shows the degree of tendency for spectral phase at each
% (time, freq) point to repeat across the data epochs. If bootstrap statistics
% are computed (as per the 'alpha' input parameter), non-significant points
% are colored green (but again not in the 'itc' output variable - use the 
% itcboot output to re-mask if desired in later plotting).
%
% The lower thin panel shows the time-domain average (ERP) of the input data 
% (blue) plus a zero-line (green). The average (ERP) is created principally by
% partial phase resetting of the EEG as measured by the ITC. (While phase resetting 
% dominates, event-related spectral power changes (as measured by the ERSP) 
% may also play a minor role). The thin left panel shows the frequency-mean ITC 
% (blue trace) and, if bootstrap statistics are computed, the ITC significance 
% limits at each frequency (green trace).
%
% ITC Math Derivation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By definition, linear coherence is
%       R=mean(Fxy)/sqrt(mean(abs(Fxx))*mean(abs(Fyy)));
% where Fxy is the cross-spectrum (FxFy*) and Fxx and Fyy the autospectra of  
% processes x and y.  We define the phase coherence to be
%       R=mean(Fxy/(abs(Fxx))*abs(Fyy));  % mean of individually normed Fxy
% To derive the ITC, we consider y to be a stimulus-locked process 
% such that, at each time and frequency, angle(Fyy) = 0 and abs(Fyy) = 1.
% Thus Pxy = Pxx, and the (complex) inter-trial phase coherence between x and
% the constant stimulus-locked process y is
%       ITC=mean(Pxx/abs(Pxx));
%
% USAGE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   >> [ersp,itc,powbase,times,freqs,erspboot,itcboot]  ...
%                = timef(data,frames,tlimits,srate,cycles,...
%                              'key1',value1,'key2',value2, ... );        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE:
% * Left-click on subplots to view and zoom in separate windows (uses axcopy()).
%
% Required Inputs:
%   data        = Single-channel data vector (1,frames*ntrials) (required)
%
%   frames      = Frames per trial                        {750}
%  Here, a frame is a data point (one channel at one time point) and the data
%  is assumed to be composed of concatenated epochs of the same length (frames).
%
%   tlimits     = [mintime maxtime] (ms) Epoch time limits {[-1000 2000]}
%  These should be the starting and ending times of the input epochs.
%
%   srate       = data sampling rate (Hz)                 {250}
%  This is sample rate per channel (i.e., data frames per second).
%
%   cycles      = >0 -> Number of cycles in each analysis wavelet 
%                 =0 -> Use FFTs (with constant window length) {0}
%                 If [wavecycles factor] -> wavelet cycles increase with frequency
%                 beginning at wavecyles (0<factor<1; factor=1 -> no increase,
%                 standard wavelets; factor=0 -> fixed epoch length, as in FFT.
%                 OR multitaper decomposition (with 'mtaper').
%  Here, the user chooses either to use the FFT method (fixed window size
%  for all frequencies) or the wavelet DFT method (the data window length
%  depends inversely on the frequency, e.g. '3' means that the data
%  windows will each be three cycles wide (at each frequency). A higher
%  number here (e.g., '5') will narrow the frequency band and widen the
%  time window.
%
% Optional Inter-Irial Coherence Type:
%   'type'      = ['coher'|'phasecoher'] Compute either linear coherence 
%                  ('coher') or phase coherence ('phasecoher') also known
%                  as the phase coupling factor           {'phasecoher'}.
%
% Optional Detrending:
%   'detret'    = ['on'|'off'], Detrend data in time.       {'off'}
%   'detrep'    = ['on'|'off'], Detrend data across trials (at each time point
%                 compute the linear trend across the set of ordered trials and 
%                 remove it; not that this also subtract the ERP) {'off'}
%
% Optional FFT/DFT Parameters:
%   'winsize'   = If cycles==0: data subwindow length (fastest is 2^n < frames);
%                 If cycles >0: *longest* window length to use. This determines 
%                 the lowest output frequency  {default: ~frames/8}
%  When cycles>0, winsize determines the lowest computed frequency. For example,
%  with srate=100 and cycles=3, a winsize of 100 means that 3 cycles must fit
%  within a 1-sec ( 100-sample) window. So, the lowest output frequency is 3 Hz.
%
%  When cycles=0, winsize is the length of data in each FFT window. This may be 
%  extended with zeroes (to give more output frequencies) using padratio (below).
%
%  'timesout'  = Number of output times (int<frames-winframes) {200}
%  The number of FFTs or wavelet DFTs computed and plotted.
%
%   'padratio'  = FFT-length/winframes (2^k)                    {2}
%                  Multiplies the number of output frequencies by
%                  dividing their spacing. When cycles==0, frequency
%                  spacing is (low_freq/padratio).
%  This factor multiplies the number of output frequencies. In the FFT method
%  (cycles=0), this is done by zero-padding each analysis window. In the wavelet
%  DFT method (cycles>0), this gives the number of frequencies per Hz.
%
%   'maxfreq'   = Maximum frequency (Hz) to plot (& to output, if cycles>0) 
%                  If cycles==0, all FFT frequencies are output. {50}
%   'baseline'  = Spectral baseline end-time (in ms). Use NaN for no baseline
%                 removal{0}
%   'powbase'   = Baseline spectrum to log-subtract. 'baseline' parameter is
%                 ignored if this parameter is used {def|NaN->from data}
%  This is useful only when you want to use a known baseline spectrum (e.g. from
%  another condition) instead of using the actual mean baseline spectrum of the data.
%  Otherwise, leave this out or specify as 'NaN' (not a number).
%
% Optional Multitaper Parameters:
%   'mtaper'    = If [N W], performs multitaper decomposition. 
%                  (N is the time resolution and W the frequency resolution; 
%                  maximum taper number is 2NW-1). Overwrites 'winsize' and 
%                  'padratio'. 
%                  If [N W K], forces the use of K Slepian tapers (if possible).
%                  Phase is calculated using standard methods.
%                  The use of mutitaper with wavelets (cycles>0) is not 
%                  recommended (as multiwavelets are not implemented). 
%                  Uses Matlab functions DPSS, PMTM.   {no multitaper}
%
% Optional Bootstrap Parameters:
%   'alpha'     = If non-0, compute two-tailed bootstrap significance prob. 
%                  level. Show non-signif. output values as green.   {0}
%  This optional parameter lengthens the computation time but gives bootstrap
%  estimates of which ERSP and ITC values are significantly different from 0,
%  by setting to 0 (green) all non-significant values in the ERSP and ITC images.
%  Normal values for alpha are 0 ([], or none) -> no bootstrap computation, or
%  0.01 (which should allow about 1% of random images to appear "significant" 
%
%   'naccu'     = Number of bootstrap replications to accumulate     {200}
%   'baseboot'  = Bootstrap baseline subtract (0 -> use 'baseline';
%                                                  1 -> use whole trial) {0}
% Optional Scalp Map Plotting Parameters:
%   'topovec'   = Scalp topography (map) to plot                     {none}
%   'elocs'     = Electrode location file for scalp map   {no default}
%                     File should be ascii in format of  >> topoplot example   
%  This is an optional map-plotting feature. Given an input map vector 
%  (one weight at each channel, and a electrode location file, timef() plots
%  a topoplot()-style 2-d scalp map on the left side of the figure. See
%  >> topoplot example % for the format of the electrode location file.
%
% Other Optional Plotting Parameters:
%   'vert'      = [vector of ms times] -> plot vertical dashed lines  {0 only}
%  Use this to add extra vertical dashed lines at significant epoch times.
%  Time 0 is marked by default.
%                     
%   'plotersp'  = ['on'|'off'] Plot power spectral perturbations    {'on'} 
%   'plotitc'   = ['on'|'off'] Plot inter trial coherence            {'on'}
%   'title'     = Optional figure title                              {none}
%
%   'pboot'     = Bootstrap power limits (e.g., from timef())   {from data}
%   'rboot'     = Bootstrap ITC limits (e.g., from timef())     {from data}
%  These are useful if you want to apply significance limits from another condition 
%  to new data. {default|NaN, compute from data}
%
%   'linewidth' = Line width for 'marktimes' traces (thick=2, thin=1) {2}
%   'axesfont'  = Axes text font size                                {10}
%   'titlefont' = Title text font size                               {8}
%
% Outputs: 
%        ersp   = Matrix (nfreqs,timesout) of log spectral diffs. from baseline (dB) 
%        itc    = Matrix (nfreqs,timesout) of inter-trial phase coherence (range: [0 1])
%   Note that when cycles=0, nfreqs is total number of FFT frequencies, which 
%   typically include frequencies higher than maxfreq. When cycles>0, *no* extra 
%   (higher) frequencies are computed.
%
%      powbase  = Baseline power spectrum (removed to compute the ERSP)
%        times  = Vector of output times (subwindow centers) (in ms).
%        freqs  = Vector of frequency bin centers (in Hz).
%
%     erspboot  = Matrix (2,nfreqs) of [lower;upper] ERSP significance diffs.
%      itcboot  = Matrix (2,nfreqs) of [lower;upper] ITC thresholds (abs., not diffs)
%  Note that the itcboot lower bound is practically meaningless.
%
%  Plot description:
%    Assuming both 'plotersp' and 'plotitc' options are 'on' (= default). The upper panel
%    presents the data ERSP (Event-Related Spectral Perturbation) in dB, with mean baseline
%    spectral activity (in dB) subtracted. Use "'baseline', NaN" to prevent timef() from
%    removing the baseline. The lower panel presents the data ITC (Inter-Trial Coherence).
%    Click on any plot axes to pop up a new window (using 'axcopy()')
%    -- Upper left marginal panel presents the mean spectrum during the baseline period
%       (blue), and when significance is set, the significance threshold at each frequency
%       (dotted green-black trace).
%    -- The marginal panel under the ERSP image shows the maximum (green) and minimum
%       (blue) ERSP values relative to baseline power at each frequency.
%    -- The lower left marginal panel shows mean ITC across the imaged time range (blue),
%       and when significance is set, the significance threshold (dotted green-black).
%    -- The marginal panel under the ITC image shows the ERP (which is produced by ITC
%       across the data spectral pass band).
%
% Author: Sigurd Enghoff, Arnaud Delorme & Scott Makeig
%          CNL / Salk Institute 1998- | SCCN/INC, UCSD 2002-
%
% See also: crossf() - event-related cross-spectral coherence between two input
%                      time series.
%
% History: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% timef() was coded by Sigurd Enghoff and Scott Makeig at The Salk 
% Institute, La Jolla CA in August, 1998, using methods developed 
% in Makeig, 1993. Arno Delorme added the multitaper option, recoded
% the function to use 'keyword','parameter' argument pairs, and added the
% 'type' argument with advice from Joern Anemueller at SCCN/Institute for
% Neural Computation, UCSD in early 2002.

% Copyright (C) 8/01/00 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 01-25-02 reformated help & license -ad 
