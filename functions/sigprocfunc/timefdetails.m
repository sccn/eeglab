% timefdetails() - details of the timef() function parameters for 
%                  time/frequency analysis of event-related 1-channel 
%                  multi-epoch data.
%
% Global description:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% timef() performs either FFT-based or wavelet DFT time/frequency
% averaging. The wavelets are N-cycle Hanning-windowed sinusoids. 
% (Note: To change the windowing to gauss(), etc., see line 258 of timef.m).
%
% The two image panels of the output plot show, respectively, the
% event-related spectral perturbation (ERSP) and inter-trial coherence (ITC)
% of the input data.
%
% THE ERSP: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (S. Makeig, Electroencephalogr Clin Neurophysiol 86:283-93, 1993) 
% shows log deviations from epoch-mean (or baseline period-mean) power at
% each frequency. If bootstrap statistics are computed (via the alpha
% input parameter), points with non-significant differences from 0
% are colored green (0, but not in the output ersp variable - use the
% output erspboot variable to re-mask if desired).
%
% The thin bottom panel shows the ERSP envelope (the most positive and most 
% negative values at each output time point).
%
% The thin left panel shows the mean (or baseline) log spectrum (blue trace).
% If bootstrap statistics are computed, the left panel (green trace) also shows 
% the bootstrap significant levels (+/-) at each frequency.
%
% THE ITC: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (Inter-trial Coherence, cf. Tallon-Buadry "Phase-locking factor")
% The lower panel shows the degree of tendency for the data phase at each
% frequency to be in the same phase in each trial. If bootstrap statistics
% are computed (as per the alpha input parameter), non-significant points
% are colored green (0 -- but not in the itc output variable - use the 
% itcboot output to re-mask if desired in later plotting).
%
% The lower thin panel shows the time-domain average of the input data (blue)
% against a zero-line (green). The output average is created principally by
% the resets to constant phase shown in the ITC. (While phase resets dominate,
% amplitude changes shown in the ERSP may also play a role).
%
% The left thin panel shows the frequency mean ITC (blue trace) and, if
% bootstrap statistics are computed, the ITC significance limits at each 
% frequency (green trace).
%
% ITC math: By definition, coherence R = mean(Pxy/(abs(Pxx)*abs(Pyy))) 
% where Pxy is the cross-spectrum and Pxx and Pyy the autospectra of  
% processes x and y. Here, we consider y to be a stimulus-locked process 
% such that, at each time and frequency, angle(Pyy) = 0 and abs(Pyy) = 1.
% Thus Pxy = Pxx, and the (complex) inter-trial coherence between x and
% the completely stimulus-locked process y is ITC = mean(Pxx/abs(Pxx)).
%
% USAGE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   >> [ersp,itc,powbase,times,freqs,erspboot,itcboot] =          ...
%                timef(data,frames,tlimits,titl,                  ...
%                         srate,cycles,winsize,timesout,          ...
%                            padratio,maxfreq,tvec,eloc_file,     ...
%                               alpha,marktimes,powbase,pboot,rboot);
%
% INPUTS:                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       data        = single-channel (1,frames*nepochs) data  {none}
%       frames      = frames per epoch                        {768}. 
%  Here, a frame is a data point (one channel at one time point) and the data
%  is assumed to be composed of concatenated epochs of the same length (frames).
%
%       tlimits     = epoch time limits (ms) [mintime maxtime]{[-1000 2000]}
%  These should be the starting and ending times of the input epochs.
%
%       titl        = figure title                            {none}
%       srate       = data sampling rate (Hz)                 {256}
%  This is sample rate per channel (i.e., data frames per second).
%
%       cycles      = >0 -> number of cycles in each analysis window 
%                     =0 -> use FFT (constant window length)  {0}
%  Here, the user chooses either to use the FFT method (fixed window size
%  for all frequencies) or the wavelet DFT method (the data window length
%  depends inversely on the frequency, e.g. '3' means that the data
%  windows will each be three cycles wide (at each frequency). A higher
%  number here (e.g., '5') will narrow the frequency band and widen the
%  time window.
%
%       winsize     = subepoch length: determines the lowest output frequency  
%                     cycles=0: data subwindow length (2^k<frames)
%                     cycles>0: *longest* subwindow length to use {~frames/8}
%  When cycles>0, winsize determines the lowest computed frequency. For example,
%   with srate=100 and cycles=3, a winsize of 100 means that 3 cycles must fit
%   within a 1-sec ( 100-sample) window. So, the lowest output frequency is 3 Hz.
%  When cycles=0, winsize is the length of data in each FFT window. This may be 
%   extended with zeroes (to give more output frequencies) using padratio (below).
%
%       timesout    = number of output times (int<frames-winsize){200}
%  The number of FFTs or wavelet DFTs computed and plotted.
%
%       padratio    = if cycles==0, FFT-length/winsize (2^k)  {2}
%                     if cycles >0, analysis frequencies per Hz. 
%  This factor multiplies the number of output frequencies. In the FFT method
%  (cycles=0), this is done by zero-padding each analysis window. In the wavelet
%  DFT method (cycles>0), this gives the number of frequencies per Hz.
%
%       maxfreq     = maximum frequency to plot (Hz)          {50}
%
%       tvec        = scalp topography (map) to plot          {[]}
%       eloc_file   = electrode location file for scalp map   {no default}
%                     ascii file in format of  >> topoplot example   
%  This is an optional map-plotting feature. Given an input map vector 
%  (one weight at each channel, and a electrode location file, timef() plots
%  a topoplot()-style 2-d scalp map on the left side of the figure.
%
%       alpha       = Two-tailed bootstrap significance prob. level {none}
%                     Sets n.s. plotted output values to green (0). 
%  This is optional, and lengthens the computation time. It gives bootstrap
%  estimates of which ERSP and ITC values are significantly different from 0,
%  by setting to 0 (green) all non-significant values in the ERSP and ITC images.
%  Normal values for alpha are 0 (or no argument) -> no bootstrap computation, or
%                              0.02 (which should allow about 2% of random images
%                                           to appear "significant" (non-green))
%
%       marktimes   = times to mark with a dotted vertical line{def|nan->none}
%  Use this to add extra vertical dashed lines at significant times in the epoch.
%  Time 0 is marked by default.
%
%       powbase     = baseline spectrum to log-subtract       {def|nan->from data}
%  This is useful only when you want to use a known baseline spectrum (e.g. from
%  another condition) instead of using the actual mean baseline spectrum of the data.
%  Otherwise, leave this out or as 'nan' (not a number).
%
%       pboot       = bootstrap power limits (cf. timef() out){def|nan -> from data}
%       rboot       = bootstrap ITC limits (cf. timef() out)  {def|nan -> from data}
%  This is useful if you want to apply significance limits from another condition to
%  new data. Else, leave these out or as 'nan'.
%
% OUTPUTS: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ersp        = log spectral differences from baseline (nfreqs,timesout)
%       itc         = inter-trial coherencies (nfreqs,timesout)
%   Note that when cycles=0, nfreqs is total number of FFT frequencies, which may 
%   include frequencies higher than maxfreq. When cycles>0, no extra (higher) 
%   frequencies are computed.
%
%          powbase  = baseline power spectrum (in whole epoch or given baseline)
%            times  = vector of output times (subwindow centers) in ms.
%            freqs  = vector of frequency bin centers in Hz.
%         erspboot  = [2,nfreqs] matrix of [lower;upper] ERSP significance diffs.
%          itcboot  = [2,nfreqs] matrix of [lower;upper] ITC thresholds (not diffs).
%  Note that the itcboot lower bound is practically meaningless.
%
% CONSTANTS SET NEAR THE TOP OF timef.m:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT_ERSP       = 1;      % Flag plot of spectral perturb. (ERSP) (1/0)
% PLOT_ITC        = 1;      % Flag plot of inter-trial coherence (1/0)
% LINEWIDTH       = 2;      % Plot thick (2) or thin (1) traces (can be non-int)
% BASE_BOOT       = 0;      % 0 = bootstrap ERSP data drawn from whole epoch (def)
%                           % 1 = bootstrap ERSP data drawn from baseline only
% BASELINE_END    = 200;    % Window center times < this are in the baseline.
% NACCU           = 200;    % Number of bootstrap repetitions to accumulate
% ERSP_CAXIS_LIMIT = 5;     % 0 -> use data limits; else positive value
%                           %         giving symmetric +/- caxis limits.
% ITC_CAXIS_LIMIT  = 0.5;   % 0 -> use data limits; else positive value
%                           %         giving symmetric +/- caxis limits.
% AXES_FONT       = 10;     % Axis text FontSize
% TITLE_FONT      = 8;      % Figure Title Font
%
% History: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% timef() was coded by Sigurd Enghoff and Scott Makeig at The Salk 
% Institute, La Jolla CA in August, 1998, using methods developed 
% in the first reference above.

% Copyright (C) 8/01/00 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 01-25-02 reformated help & license -ad 
