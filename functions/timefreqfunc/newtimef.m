% WARNING: this function is not part of the EEGLAB toolbox and should not be distributed
%          you must contact Arnaud Delorme (arno@salk.edu) for terms of use
%
% newtimef() - Returns estimates and plots of event-related (log) spectral
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
%                     To compare conditions 1 (data1) and condition 2 (data2)
%                     in place of data enter { data1 data2 }
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
%       'detrend'   = ['on'|'off'], Linearly detrend each data epoch   {'off'}
%       'rmerp'     = ['on'|'off'], Remove epoch mean from data epochs {'off'}
%
%    Optional FFT/DFT Parameters:
%       'winsize'   = If cycles==0: data subwindow length (fastest, 2^n<frames);
%                     If cycles >0: *longest* window length to use. This
%                      determines the lowest output frequency. Note that this
%                     parameter is overwritten if the minimum frequency requires
%                     a longer time window {~frames/8}
%       'timesout'  = Number of output times (int<frames-winframes). Enter a 
%                     negative value [-S] to subsample original time by S.
%                     Enter an array to obtain spectral decomposition at 
%                     specific time values (note: algorithm find closest time 
%                     point in data and this might result in an unevenly spaced
%                     time array. {def: 200}
%       'padratio'  = FFT-length/winframes (2^k)                    {2}
%                     Multiplies the number of output frequencies by dividing
%                     their spacing (standard FFT padding). When cycles~=0, 
%                     frequency spacing is divided by padratio.
%       'maxfreq'   = Maximum frequency (Hz) to plot (& to output, if cycles>0) 
%                     If cycles==0, all FFT frequencies are output. {50}
%                     DEPRECATED, use 'freqs' instead,
%       'freqs'     = [min max] frequency limits. Default [minfreq 50], 
%                     minfreq being determined by the number of data points, 
%                     cycles and sampling frequency.
%       'nfreqs'    = number of output frequencies. For FFT, closest computed
%                     frequency will be returned. Overwrite 'padratio' effects
%                     for wavelets. Default: use 'padratio'.
%       'freqscale' = ['log'|'linear'] frequency scale. Default is 'linear'.
%                     Note that for obtaining 'log' spaced freqs using FFT, 
%                     closest correspondant frequencies in the 'linear' space 
%                     are returned.
%       'baseline'  = Spectral baseline end-time (in ms). NaN imply that no
%                      baseline is used. A range [min max] may also be entered
%                     You may also enter one row per region for baseline
%                     e.g. [0 100; 300 400] considers the window 0 to 100 ms and
%                     300 to 400 ms. { default 0 }
%       'powbase'   = Baseline spectrum to log-subtract. {def|NaN->from data}
%       'lowmem'    = ['on'|'off'] compute frequency, by frequency to save
%                     memory. Default 'off'.
%       'verbose'   = ['on'|'off'] print text {'on'}
%       'subitc'    = ['on'|'off'] subtract stimulus locked Inter-Trial Coherence
%                     (ITC) from x and y. This computes the  'intrinsic' coherence
%                     x and y not arising from common synchronization to
%                     experimental events. See notes. {default: 'off'}
%
%    Optional Bootstrap Parameters:
%       'alpha'     = If non-0, compute two-tailed bootstrap significance prob. 
%                      level. Show non-signif. output values as green.   {0}
%       'naccu'     = Number of bootstrap replications to accumulate     {200}
%       'baseboot'  = Bootstrap baseline subtract (1 -> use 'baseline';
%                                                  0 -> use whole trial
%                                                  [min max] -> use time range) {1}
%                     You may also enter one row per region for baseline
%                     e.g. [0 100; 300 400] considers the window 0 to 100 ms and
%                     300 to 400 ms.
%       'boottype'  = ['shuffle'|'rand'|'randall'] shuffle time and trials or
%                     invert polarity of spectral data points in ERSP (or
%                     randomize phase in ITC) ('rand'). 'randall' compute 
%                     significance by accumulating radom polarity invertions 
%                     for each times/frequencies points (time consuming).
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
%       'plottype'  = ['image'|'curve'] plot time frequency images or
%                     curves (one curve per frequency). Default is 'image'.
%       'plotmean'  = ['on'|'off'] For 'curve' plots only. Average all
%                     frequencies given as input. Default: 'on'.
%       'highlightmode'  = ['background'|'bottom'] For 'curve' plots only,
%                     display significant time regions either in the plot background
%                     or underneatht the curve.
%       'plotersp'  = ['on'|'off'] Plot power spectral perturbations              {'on'} 
%       'plotitc'   = ['on'|'off'] Plot inter trial coherence                     {'on'}
%       'plotphasesign' = ['on'|'off'] Plot phase sign in the inter trial coherence {'on'}
%       'plotphase' = ['on'|'off'] Plot ITC phase instead ofITC amplitude         {'off'}
%       'erspmax'   = [real dB] set the ERSP max. for the color scale (min= -max) { auto }
%       'itcmax'    = [real] set the ITC image maximum for the color scale        { auto }
%       'erplim'    = [min max] ERP limits for ITC (below ITC image)              { auto }
%       'itcavglim' = [min max] average ITC limits for all freq. (left of ITC)    { auto }
%       'speclim'   = [min max] average spectrum limits (left of ERSP image)      { auto }
%       'erspmarglim' = [min max] average marginal ERSP limits (below ERSP image) { auto }
%       'title'     = Optional figure title                                       { none }
%       'marktimes' = Non-0 times to mark with a dotted vertical line (ms)        { none }
%       'linewidth' = Line width for 'marktimes' traces (thick=2, thin=1)         {2}
%       'axesfont'  = Axes text font size                                         {10}
%       'titlefont' = Title text font size                                        {8}
%       'vert'      = [times_vector] -> plot vertical dashed lines at specified times
%                     in ms.
%       'newfig'    = ['on'|'off'] Create new figure for difference plots {'on'}
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
% Revision 1.76  2005/07/27 18:23:34  arno
% keyword
%
% Revision 1.75  2005/07/27 18:21:24  arno
% same
%
% Revision 1.74  2005/07/27 18:17:54  arno
% implementing subitc
%
% Revision 1.73  2005/04/08 22:33:05  arno
% allowing to set limits; common limits for condition difference
%
% Revision 1.72  2005/04/07 02:02:22  arno
% plotting mbase
% ,
%
% Revision 1.71  2005/04/07 01:57:57  arno
% defining mbase
%
% Revision 1.70  2005/04/07 01:51:54  arno
% allowing to set powbase
%
% Revision 1.69  2005/03/07 21:25:47  arno
% chaninfo
%
% Revision 1.68  2005/01/28 00:52:01  arno
% same
%
% Revision 1.67  2005/01/28 00:47:41  arno
% fix timesout
%
% Revision 1.66  2004/12/10 23:07:00  arno
% add text message
%
% Revision 1.65  2004/12/07 22:29:58  arno
% new version with new bootstrap
%
% Revision 1.64  2004/10/25 18:37:02  hilit
% 'ploterps' -> 'plotersp'
%
% Revision 1.63  2004/08/10 01:13:33  arno
% cheking baseline
%
% Revision 1.62  2004/08/03 21:55:57  arno
% same scaling for both iamges
%
% Revision 1.61  2004/06/02 23:18:03  arno
% baseline, baseboot ...
%
% Revision 1.60  2004/06/02 22:50:19  arno
% do not remember
%
% Revision 1.59  2004/06/01 23:21:59  arno
% implementing baseboot and baseline windows
%
% Revision 1.58  2004/03/09 17:35:04  arno
% msg
%
% Revision 1.57  2004/03/01 16:18:20  arno
% default basebott 1
%
% Revision 1.56  2004/03/01 02:24:31  arno
% help message
%
% Revision 1.55  2004/02/27 18:59:22  arno
% reshape data
%
% Revision 1.54  2004/02/14 23:24:05  arno
% implement erspmax
%
% Revision 1.53  2004/01/06 17:01:06  arno
% header typo
%
% Revision 1.52  2003/12/09 23:23:31  arno
% nothing
%
% Revision 1.51  2003/12/09 23:13:03  arno
% fixing plotimef baseline problem
%
% Revision 1.50  2003/12/08 18:32:01  arno
% undoing baseline
%
% Revision 1.49  2003/12/08 18:10:44  arno
% problem in baseline subtraction
%
% Revision 1.48  2003/12/03 02:32:12  arno
% document verbose
%
% Revision 1.47  2003/12/03 02:31:34  arno
% verbose
%
% Revision 1.46  2003/10/15 18:46:25  arno
% *** empty log message ***
%
% Revision 1.45  2003/10/15 18:45:11  arno
% *** empty log message ***
%
% Revision 1.44  2003/08/08 22:19:10  arno
% *** empty log message ***
%
% Revision 1.43  2003/08/04 16:37:14  arno
% updating color for significance curves
%
% Revision 1.42  2003/08/04 14:43:09  arno
% updating ERP plot as for timef
%
% Revision 1.41  2003/08/01 21:05:12  arno
% commenting h(9)
%
% Revision 1.40  2003/07/09 21:33:17  arno
% timesout - > ntimesout for timefreq call
%
% Revision 1.39  2003/06/27 01:01:48  arno
% updating timesout help
%
% Revision 1.38  2003/06/19 01:22:09  arno
% debuging lowmem for multiple condition and 'coher'
%
% Revision 1.37  2003/06/19 00:16:06  arno
% nothing
%
% Revision 1.36  2003/06/18 00:37:02  arno
% debuging lowmem
%
% Revision 1.35  2003/06/17 23:19:12  arno
% j -> jj
%
% Revision 1.34  2003/05/29 15:07:03  arno
% debugging lowmem
%
% Revision 1.33  2003/05/24 19:09:16  arno
% debug lowmem
%
% Revision 1.32  2003/05/22 01:21:30  arno
% detrending param name change
%
% Revision 1.31  2003/05/21 02:26:16  arno
% debug condition
%
% Revision 1.30  2003/05/20 22:29:05  arno
% lowmem for 2condition debug
%
% Revision 1.29  2003/05/20 21:44:08  arno
% implemementing lowmem for 2 conditions
%
% Revision 1.28  2003/05/16 15:55:44  arno
% debuging ticks
%
% Revision 1.27  2003/05/02 21:56:27  arno
% adding low mem option
%
% Revision 1.26  2003/05/02 17:58:56  arno
% do not pass on maxfreq to timefreq
%
% Revision 1.25  2003/04/30 00:07:28  arno
% removing all ref to dispf
% /
%
% Revision 1.24  2003/04/29 18:41:06  arno
% implementing log vizualization
%
% Revision 1.23  2003/04/17 18:06:06  arno
% adding newfig option
%
% Revision 1.22  2003/01/20 05:30:42  cooper
% corrected ntrials in linear coher formula.
%
% Revision 1.21  2003/01/11 00:38:02  cooper
% fixed formula for linear coher.
%
% Revision 1.20  2003/01/10 19:54:29  arno
% debugging linear coherence bootstrap for 2 conditions
%
% Revision 1.19  2003/01/08 23:34:33  arno
% typo in bootstrap formula for coherence
%
% Revision 1.18  2003/01/02 03:15:27  cooper
% corrected text output msg.
%
% Revision 1.17  2002/11/20 01:33:39  arno
% crossf diff no plot
%
% Revision 1.16  2002/11/20 01:13:49  arno
% remove duplicate arguments
%
% Revision 1.15  2002/11/20 01:00:59  arno
% undo last
%
% Revision 1.14  2002/11/20 01:00:37  arno
% nothing
%
% Revision 1.13  2002/10/24 02:16:28  arno
% output of the right size
%
% Revision 1.12  2002/10/17 01:49:49  arno
% title check
%
% Revision 1.11  2002/10/15 23:02:28  arno
% time limits warning
%
% Revision 1.10  2002/10/15 20:54:16  arno
% title diff
%
% Revision 1.9  2002/10/15 19:17:25  arno
% cyclefact
%
% Revision 1.8  2002/10/15 19:07:11  arno
% for no signif diff
%
% Revision 1.7  2002/10/15 18:38:16  arno
% adding alpha to bootstat
%
% Revision 1.6  2002/10/15 18:27:08  arno
% missing argument
%
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

function [P,R,mbase,timesout,freqs,Pboot,Rboot,alltfX,PA] = timef( X, frame, tlimits, Fs, varwin, varargin);

% Note: PA is output of 'phsamp','on' 

%varwin,winsize,g.timesout,g.padratio,g.maxfreq,g.topovec,g.elocs,g.alpha,g.marktimes,g.powbase,g.pboot,g.rboot)

% ITC:   Normally, R = |Sum(Pxy)| / (Sum(|Pxx|)*Sum(|Pyy|)) is coherence.
%        But here, we consider    Phase(Pyy) = 0 and |Pyy| = 1 -> Pxy = Pxx
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
	help newtimef
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
    [tmp indices] = unique(varargin(1:2:end));
    varargin = varargin(sort(union(indices*2-1, indices*2))); % these 2 line remove duplicate arguments
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
try, g.boottype;   catch, g.boottype = 'shuffle'; end;
try, g.condboot;   catch, g.condboot = 'abs'; end;
try, g.title;      catch, g.title = DEFAULT_TITLE; end;
try, g.winsize;    catch, g.winsize = max(pow2(nextpow2(g.frame)-3),4); end;
try, g.pad;        catch, g.pad = max(pow2(nextpow2(g.winsize)),4); end;
try, g.timesout;   catch, g.timesout = DEFAULT_NWIN; end;
try, g.padratio;   catch, g.padratio = DEFAULT_OVERSMP; end;
try, g.topovec;    catch, g.topovec = []; end;
try, g.elocs;      catch, g.elocs = DEFAULT_ELOC; end;
try, g.alpha;      catch, g.alpha = DEFAULT_ALPHA; end;  
try, g.marktimes;  catch, g.marktimes = DEFAULT_MARKTIME; end;
try, g.powbase;    catch, g.powbase = NaN; end;
try, g.pboot;      catch, g.pboot = NaN; end;
try, g.rboot;      catch, g.rboot = NaN; end;
try, g.plotersp;   catch, g.plotersp = 'on'; end;
try, g.subitc;     catch, g.subitc  = 'off'; end;
try, g.plotitc;    catch, g.plotitc  = 'on'; end;
try, g.detrend;    catch, g.detrend = 'off'; end;
try, g.rmerp;      catch, g.rmerp = 'off'; end;
try, g.baseline;   catch, g.baseline = 0; end;
try, g.baseboot;   catch, g.baseboot = 1; end;
try, g.linewidth;  catch, g.linewidth = 2; end;
try, g.naccu;      catch, g.naccu = 200; end;
try, g.mtaper;     catch, g.mtaper = []; end;
try, g.maxfreq;    catch, g.maxfreq = DEFAULT_MAXFREQ; end;
try, g.freqs;      catch, g.freqs = [0 g.maxfreq]; end;
try, g.nfreqs;     catch, g.nfreqs = []; end;
try, g.freqscale;  catch, g.freqscale = 'linear'; end;
try, g.vert;       catch, g.vert = []; end;
try, g.newfig;     catch, g.newfig = 'on'; end;
try, g.type;       catch, g.type = 'phasecoher'; end;
try, g.phsamp;     catch, g.phsamp = 'off'; end;
try, g.plotphase;  catch, g.plotphase = 'off'; end;
try, g.plotphasesign;  catch, g.plotphasesign = 'on'; end;
try, g.outputformat;  catch, g.outputformat = 'new'; end;
try, g.itcmax;     catch, g.itcmax = []; end;
try, g.erspmax;    catch, g.erspmax = []; end;
try, g.lowmem;     catch, g.lowmem = 'off'; end;
try, g.verbose;    catch, g.verbose = 'on'; end;
try, g.plottype;   catch, g.plottype = 'image'; end;
try, g.plotmean;   catch, g.plotmean = 'on'; end;
try, g.highlightmode;   catch, g.highlightmode = 'background'; end;
try, g.chaninfo;        catch, g.chaninfo    = []; end; 
try, g.erspmarglim;     catch, g.erspmarglim = []; end;
try, g.itcavglim;       catch, g.itcavglim   = []; end;
try, g.erplim;          catch, g.erplim      = []; end;
try, g.speclim;         catch, g.speclim     = []; end;
g.AXES_FONT       = AXES_FONT;           % axes text FontSize
g.TITLE_FONT      = TITLE_FONT;
g.ERSP_CAXIS_LIMIT = ERSP_CAXIS_LIMIT;         
g.ITC_CAXIS_LIMIT  = ITC_CAXIS_LIMIT;        

if isfield(g, 'detret'), g.detrend = g.detret; end;
if isfield(g, 'detrep'), g.rmerp   = g.detrep; end;

% testing arguments consistency
% -----------------------------
switch lower(g.verbose)
    case { 'on', 'off' }, ;
    otherwise error('verbose must be either on or off');
end;
if ~ischar(g.title) & ~iscell(g.title)
	error('Title must be a string or a cell array.');
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

if ~all(isnumeric(g.timesout))
	error('Value of timesout must be a number.');
end
if length(g.timesout) == 1 & g.timesout > 0
    if g.timesout > g.frame-g.winsize
        g.timesout = g.frame-g.winsize;
        disp(['Value of timesout must be <= frame-winsize, timeout adjusted to ' int2str(g.timesout) ]);
    end
end;

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
	myprintf(g.verbose, ['Warning: value of maxfreq reduced to Nyquist rate' ...
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
if ~strcmpi(g.boottype, 'shuffle') & ~strcmpi(g.boottype, 'rand')  & ~strcmpi(g.boottype, 'randall')
	error('Boottype must be either ''times'' or ''timestrials''.');
end;	
if ~strcmpi(g.condboot, 'abs') & ~strcmpi(g.condboot, 'angle') ...
		& ~strcmpi(g.condboot, 'complex')
	error('Condboot must be either ''abs'', ''angle'' or ''complex''.');
end;
if (~isnumeric(g.alpha) | length(g.alpha)~=1)
	error('timef(): Value of g.alpha must be a number.\n');
elseif (round(g.naccu*g.alpha) < 2)
	myprintf(g.verbose, 'Value of g.alpha is out of the normal range [%g,0.5]\n',2/g.naccu);
    g.naccu = round(2/g.alpha);
	myprintf(g.verbose, '  Increasing the number of bootstrap iterations to %d\n',g.naccu);
end
if g.alpha>0.5 | g.alpha<=0
    error('Value of g.alpha is out of the allowed range (0.00,0.5).');
end
if ~isnan(g.alpha)
    if length(g.baseboot) == 2
     myprintf(g.verbose, 'Bootstrap analysis will use data in range %3.2g-%3.2g ms.\n', g.baseboot(1),  g.baseboot(2))
    elseif g.baseboot > 0
     myprintf(g.verbose, 'Bootstrap analysis will use data in baseline (pre-0) subwindows only.\n')
   else
     myprintf(g.verbose, 'Bootstrap analysis will use data in all subwindows.\n')
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

switch lower(g.plotphasesign)
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
switch lower(g.rmerp)
    case { 'on', 'off' }, ;
    otherwise error('rmerp must be either on or off');
end;
switch lower(g.detrend)
    case { 'on', 'off' }, ;
    otherwise error('detrend must be either on or off');
end;
switch lower(g.phsamp)
    case { 'on', 'off' }, ;
    otherwise error('phsamp must be either on or off');
end;
switch lower(g.newfig)
    case { 'on', 'off' }, ;
    otherwise error('newfig must be either on or off');
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
switch g.type
    case { 'coher', 'phasecoher', 'phasecoher2' },;
    otherwise error('Type must be either ''coher'' or ''phasecoher''');
end;    
if g.tlimits(2)-g.tlimits(1) < 30
    disp('Crossf WARNING: time range is very small (<30 ms). Times limits are in millisenconds not seconds.'); 
end;

% checking keywords
% -----------------
allfields = { 'tlimits' 'frame' 'srate' 'cycles' 'cyclesfact' 'boottype' 'condboot' 'title' 'winsize' 'pad' 'timesout' 'padratio' 'topovec' 'elocs' 'alpha' 'marktimes' 'powbase' 'pboot' 'rboot' 'plotersp' 'plotitc' 'detrend' 'rmerp' 'baseline' 'baseboot' 'linewidth' 'naccu' 'mtaper' 'maxfreq' 'freqs' ...
              'plotamp' 'subitc' 'nfreqs' 'freqscale' 'vert' 'newfig' 'type' 'phsamp' 'plotphase' 'plotphasesign' 'outputformat' 'itcmax' 'erspmax' 'lowmem' 'verbose' 'plottype' 'plotmean' 'highlightmode' 'chaninfo' 'erspmarglim' 'itcavglim' 'erplim' 'speclim' 'AXES_FONT' 'TITLE_FONT' 'ERSP_CAXIS_LIMIT' 'ITC_CAXIS_LIMIT' };
tmpfields = fieldnames(g);
for index = 1:length(tmpfields)
    if isempty(strmatch(tmpfields{index}, allfields))
        error([ 'Unknown keyword ''' tmpfields{index} '''']);
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute frequency by frequency if low memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(g.lowmem, 'on') & length(X) ~= g.frame & isempty(g.nfreqs) & ~iscell(X)
    
    % compute for first 2 trials to get freqsout
    XX = reshape(X, 1, frame, prod(size(X))/g.frame);    
    [P,R,mbase,timesout,freqsout] = newtimef(XX(1,:,1), frame, tlimits, Fs, varwin, 'plotitc', 'off', 'plotamp', 'off',varargin{:}, 'lowmem', 'off');
    
    % scan all frequencies
    for index = 1:length(freqsout)
        if nargout < 8
            [P(index,:),R(index,:),mbase(index),timesout,tmpfreqs(index),Pboot(index,:),Rboot(index,:)] = ...
                newtimef(X, frame, tlimits, Fs, varwin, 'freqs', [freqsout(index) freqsout(index)], 'nfreqs', 1, ...
                          'plotamp', 'off', 'plotphasesign', 'off',varargin{:}, 'lowmem', 'off', 'timesout', timesout);
        else
            [P(index,:),R(index,:),mbase(index),timesout,tmpfreqs(index),Pboot(index,:),Rboot(index,:), ...
            alltfX(index,:,:)] = ...
                newtimef(X, frame, tlimits, Fs, varwin, 'freqs', [freqsout(index) freqsout(index)], 'nfreqs', 1, ...
                          'plotamp', 'off', 'plotphasesign', 'off',varargin{:}, 'lowmem', 'off', 'timesout', timesout);
        end;
    end;
    
    % plot and return
    ERP = mean(X,2);
    plottimef(P, R, Pboot, Rboot, ERP, freqsout, timesout, mbase, g);
    return;
end;    


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare 2 conditions part
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iscell(X)
    vararginori = varargin;
	if length(X) ~= 2
		error('timef: to compare conditions, data must be 2-elements cell arrays');
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
		g.title = { 'Condition 1', 'Condition 2', 'Condition 1 - condition 2' };
	end;
	
	myprintf(g.verbose, 'Running newtimef on condition 1 *********************\n');
	myprintf(g.verbose, 'Note: if an out-of-memory error occurs, try reducing the\n');
	myprintf(g.verbose, '      number of time points or number of frequencies\n');
	myprintf(g.verbose, '      (the ''coher'' options takes 3 times more memory than other options)\n');
    [P1,R1,mbase1,timesout,freqs,Pboot1,Rboot1,alltfX1] = newtimef( X{1}, frame, tlimits, Fs, varwin, ...
                                                      'plotitc', 'off', 'plotersp', 'off', vararginori{:}, 'lowmem', 'off');
    
	myprintf(g.verbose, '\nRunning newtimef on condition 2 *********************\n');
    [P2,R2,mbase2,timesout,freqs,Pboot2,Rboot2,alltfX2] = newtimef( X{2}, frame, tlimits, Fs, varwin,  ...
                                                      'plotitc', 'off', 'plotersp', 'off', vararginori{:}, 'lowmem', 'off');
    
    % recompute baselines for power
    % -----------------------------
    if ~isnan( g.baseline(1) ) & ~isnan( mbase1 ) & isnan(g.powbase)
        disp('Recomputing baseline power by using the mean from both conditions');
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
        myprintf(g.verbose, '\nSubtracting common baseline\n');
    else
        mbase = NaN;
    end;

    % plotting
    % --------
    if strcmpi(g.plotersp, 'on') | strcmpi(g.plotitc, 'on')
        g.titleall = g.title;
        if strcmpi(g.newfig, 'on'), figure; end;
        
        % using same color scale
        % ----------------------
        if ~isfield(g, 'erspmax')
            g.erspmax = max( max(max(abs(Pboot1))), max(max(abs(Pboot2))) );
        end;
        if ~isfield(g, 'itcmax')
            g.itcmax  = max( max(max(abs(Rboot1))), max(max(abs(Rboot2))) );
        end;
            
        subplot(1,3,1); g.title = g.titleall{1}; 
        g = plottimef(P1, R1, Pboot1, Rboot1, mean(X{1},2), freqs, timesout, mbase, g);
        g.itcavglim = [];
        subplot(1,3,2); g.title = g.titleall{2}; 
        plottimef(P2, R2, Pboot2, Rboot2, mean(X{2},2), freqs, timesout, mbase, g);
        subplot(1,3,3); g.title = g.titleall{3};
    end;
    
    if isnan(g.alpha)
        switch(g.condboot)
            case 'abs',  Rdiff = abs(R1)-abs(R2);
            case 'angle',  Rdiff = angle(R1)-angle(R2);
            case 'complex',  Rdiff = R1-R2;
        end;
        if strcmpi(g.plotersp, 'on') | strcmpi(g.plotitc, 'on')
            plottimef(P1-P2, Rdiff, [], [], mean(X{1},2)-mean(X{2},2), freqs, timesout, mbase, g);
        end;
	else 		
		% preprocess data and run compstat
		% --------------------------------
		alltfX1power = alltfX1.*conj(alltfX1);
		alltfX2power = alltfX2.*conj(alltfX2);
		formula = {'log10(mean(arg1(:,:,X),3))'};
		switch g.type
		 case 'coher', % take the square of alltfx and alltfy first to speed up
		  formula = { formula{1} ['sum(arg2(:,:,X),3)./sqrt(sum(arg1(:,:,X),3)*length(X) )'] };
          if strcmpi(g.lowmem, 'on')
              for ind = 1:2:size(alltfX1power,1)
                  if ind == size(alltfX1,1), indarr = ind; else indarr = [ind:ind+1]; end;
                  [resdifftmp resimagestmp res1tmp res2tmp] = ...
                      condstat(formula, g.naccu, g.alpha, {'both' 'upper'}, { '' g.condboot}, ...
                               { alltfX1power(indarr,:,:) alltfX2power(indarr,:,:) }, {alltfX1(indarr,:,:) alltfX2(indarr,:,:)});
                  resdiff{1}(indarr,:)     = resdifftmp{1};   resdiff{2}(indarr,:)     = resdifftmp{2};
                  resimages{1}(indarr,:,:) = resimagestmp{1}; resimages{2}(indarr,:,:) = resimagestmp{2};
                  res1{1}(indarr,:)        = res1tmp{1};      res1{2}(indarr,:)        = res1tmp{2};
                  res2{1}(indarr,:)        = res2tmp{1};      res2{2}(indarr,:)        = res2tmp{2};
              end;     
          else
              [resdiff resimages res1 res2] = condstat(formula, g.naccu, g.alpha, {'both' 'upper'}, { '' g.condboot}, ...
                                                       { alltfX1power alltfX2power }, {alltfX1 alltfX2});
          end; 
		 case 'phasecoher2', % normalize first to speed up
		  formula = { formula{1} ['sum(arg2(:,:,X),3)./sum(arg3(:,:,X),3)'] }; 
		  alltfX1abs = sqrt(alltfX1power); % these 2 lines can be suppressed 
		  alltfX2abs = sqrt(alltfX2power); % by inserting sqrt(arg1(:,:,X)) instead of arg3(:,:,X))
          if strcmpi(g.lowmem, 'on')
              for ind = 1:2:size(alltfX1abs,1)
                  if ind == size(alltfX1,1), indarr = ind; else indarr = [ind:ind+1]; end;
                  [resdifftmp resimagestmp res1tmp res2tmp] = ...
                      condstat(formula, g.naccu, g.alpha, {'both' 'upper'}, { '' g.condboot}, ...
                               { alltfX1power(indarr,:,:) alltfX2power(indarr,:,:) }, {alltfX1(indarr,:,:) ...
                                      alltfX2(indarr,:,:)}, { alltfX1abs(indarr,:,:) alltfX2abs(indarr,:,:) });
                  resdiff{1}(indarr,:)     = resdifftmp{1};   resdiff{2}(indarr,:)     = resdifftmp{2};
                  resimages{1}(indarr,:,:) = resimagestmp{1}; resimages{2}(indarr,:,:) = resimagestmp{2};
                  res1{1}(indarr,:)        = res1tmp{1};      res1{2}(indarr,:)        = res1tmp{2};
                  res2{1}(indarr,:)        = res2tmp{1};      res2{2}(indarr,:)        = res2tmp{2};
              end;     
          else
              [resdiff resimages res1 res2] = condstat(formula, g.naccu, g.alpha, {'both' 'upper'}, { '' g.condboot}, ...
                                                       { alltfX1power alltfX2power }, {alltfX1 alltfX2}, { alltfX1abs alltfX2abs });
          end; 
		 case 'phasecoher',
		  formula = { formula{1} ['mean(arg2(:,:,X),3)'] }; 
          if strcmpi(g.lowmem, 'on')
              for ind = 1:2:size(alltfX1,1)
                  if ind == size(alltfX1,1), indarr = ind; else indarr = [ind:ind+1]; end;
                  alltfX1norm = alltfX1(indarr,:,:)./sqrt(alltfX1(indarr,:,:).*conj(alltfX1(indarr,:,:)));
                  alltfX2norm = alltfX2(indarr,:,:)./sqrt(alltfX2(indarr,:,:).*conj(alltfX2(indarr,:,:)));
                  [resdifftmp resimagestmp res1tmp res2tmp] = ...
                      condstat(formula, g.naccu, g.alpha, {'both' 'both'}, { '' g.condboot}, ...
                               { alltfX1power(indarr,:,:) alltfX2power(indarr,:,:) }, { alltfX1norm alltfX2norm });
                  resdiff{1}(indarr,:)     = resdifftmp{1};   resdiff{2}(indarr,:)     = resdifftmp{2};
                  resimages{1}(indarr,:,:) = resimagestmp{1}; resimages{2}(indarr,:,:) = resimagestmp{2};
                  res1{1}(indarr,:)        = res1tmp{1};      res1{2}(indarr,:)        = res1tmp{2};
                  res2{1}(indarr,:)        = res2tmp{1};      res2{2}(indarr,:)        = res2tmp{2};
              end;     
          else
              alltfX1norm = alltfX1./sqrt(alltfX1.*conj(alltfX1));
              alltfX2norm = alltfX2./sqrt(alltfX2.*conj(alltfX2)); % maybe have to suppress preprocessing -> lot of memory
              [resdiff resimages res1 res2] = condstat(formula, g.naccu, g.alpha, {'both' 'both'}, { '' g.condboot}, ...
                                                       { alltfX1power alltfX2power }, { alltfX1norm alltfX2norm });
          end; 
		end;
        
		% same as below: plottimef(P1-P2, R2-R1, 10*resimages{1}, resimages{2}, mean(X{1},2)-mean(X{2},2), freqs, times, mbase, g);
        if strcmpi(g.plotersp, 'on') | strcmpi(g.plotitc, 'on')
            g.erspmax = []; % auto scale
            g.itcmax  = []; % auto scale
            plottimef(10*resdiff{1}, resdiff{2}, 10*resimages{1}, resimages{2}, ...
                      mean(X{1},2)-mean(X{2},2), freqs, timesout, mbase, g);
		end;
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
myprintf(g.verbose, 'Computing Event-Related Spectral Perturbation (ERSP) and\n');
switch g.type
    case 'phasecoher',  myprintf(g.verbose, '  Inter-Trial Phase Coherence (ITC) images based on %d trials\n',trials);
    case 'phasecoher2', myprintf(g.verbose, '  Inter-Trial Phase Coherence 2 (ITC) images based on %d trials\n',trials);
    case 'coher',       myprintf(g.verbose, '  Linear Inter-Trial Coherence (ITC) images based on %d trials\n',trials);
end;
myprintf(g.verbose, '  of %d frames sampled at %g Hz.\n',g.frame,g.srate);
myprintf(g.verbose, 'Each trial contains samples from %1.0f ms before to\n',g.tlimits(1));
myprintf(g.verbose, '  %1.0 ms after the timelocking event.\n',g.tlimits(2));
if ~isnan(g.alpha)
  myprintf(g.verbose, 'Only significant values (bootstrap p<%g) will be colored;\n',g.alpha) 
  myprintf(g.verbose, '  non-significant values will be plotted in green\n');
end

% -----------------------------------------
% detrend over epochs (trials) if requested
% -----------------------------------------
X = reshape(X, g.frame, prod(size(X))/g.frame);
if strcmpi(g.rmerp, 'on')
    X = X - mean(X,2)*ones(1, length(X(:))/g.frame);
end;        

% ----------------------------------------------------
% compute time frequency decompositions, power and ITC
% ----------------------------------------------------
if length(g.timesout) > 1, tmioutopt = { 'timesout' , g.timesout };
else                       tmioutopt = { 'ntimesout', g.timesout };
end;
[alltfX freqs timesout R] = timefreq(X, g.srate, tmioutopt{:}, 'winsize', g.winsize, ...
                                     'tlimits', g.tlimits, 'detrend', g.detrend, 'itctype', ...
                                     g.type, 'subitc', g.subitc, 'wavelet', [g.cycles g.cyclesfact], ...
                                     'padratio', g.padratio, 'freqs', g.freqs, 'freqscale', g.freqscale, ...
                                     'nfreqs', g.nfreqs); 
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
    for jj=1:g.timesout
        PA(:,:,jj) = PA(:,:,jj) ./ repmat(PP(:,jj)', [size(PP,1) 1]);
    end
end

% --------
% baseline
% --------
if size(g.baseline,2) == 2
    baseln = [];
    for index = 1:size(g.baseline,1)
        tmptime   = find(timesout >= g.baseline(index,1) & timesout <= g.baseline(index,2));
        baseln = union(baseln, tmptime);
    end;
    if length(baseln)==0
        error('No point found in baseline');
    end;
else
    if ~isempty(find(timesout < g.baseline))
        baseln = find(timesout < g.baseline); % subtract means of pre-0 (centered) windows
    else
        baseln = 1:length(timesout); % use all times as baseline
    end
end;

if ~isnan(g.alpha) & length(baseln)==0
    myprintf(g.verbose, 'timef(): no window centers in baseline (times<%g) - shorten (max) window length.\n', g.baseline)
    return
end
if isnan(g.powbase)
  myprintf(g.verbose, 'Computing the mean baseline spectrum\n');
  mbase = mean(P(:,baseln),2)';
else
  myprintf(g.verbose, 'Using the input baseline spectrum\n');
  mbase = 10.^(g.powbase/10);
end
baselength = length(baseln);
if ~isnan( g.baseline(1) ) & ~isnan( mbase )
    P = 10 * (log10(P) - repmat(log10(mbase(1:size(P,1)))',[1 length(timesout)])); % convert to (10log10) dB
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
        if size(g.baseboot,2) == 1
            if g.baseboot == 0, baselntmp = [];
            elseif ~isnan(g.baseline(1))
                 baselntmp = baseln; 
            else baselntmp = find(timesout <= 0); % if it is empty use whole epoch
            end;
        else
            baselntmp = [];
            for index = 1:size(g.baseboot,1)
                tmptime   = find(timesout >= g.baseboot(index,1) & timesout <= g.baseboot(index,2));
                if isempty(tmptime),
                    fprintf('Warning: empty baseline bootstrap interval [%3.2f %3.2f]\n', g.baseboot(index,1), g.baseboot(index,2));
                end;
                baselntmp = union(baselntmp, tmptime);
            end;
		end;
        if prod(size(g.baseboot)) > 2
            fprintf('Bootstrap analysis will use data in multiple selected windows.\n');
        elseif size(g.baseboot,2) == 2
            fprintf('Bootstrap analysis will use data in range %3.2g-%3.2g ms.\n', g.baseboot(1),  g.baseboot(2));
        elseif g.baseboot
            fprintf('   %d bootstrap windows in baseline (times<%g).\n', length(baselntmp), g.baseboot)
        end;        

        % power significance
        % ------------------
        formula = 'mean(arg1,3);';
        inputdata = alltfX.*conj(alltfX);
        if ~isnan( g.baseline ) 
           inputdata = inputdata ./ repmat(mbase', [1 size(inputdata,2) size(inputdata,3)]);
        end;
        Pboot = bootstat(inputdata, formula, 'boottype', 'shuffle', ...
							   'label', 'ERSP', 'bootside', 'both', 'naccu', g.naccu, ...
                               'basevect', baselntmp, 'alpha', g.alpha, 'dimaccu', 2 );
    	Pboot = 10*log10(Pboot);
                     
        % ITC significance
        % ----------------
        inputdata = alltfX;
		switch g.type
		 	case 'coher',       formula = [ 'sum(arg1,3)./sqrt(sum(arg1.*conj(arg1),3))/ sqrt(' int2str(trials) ');' ];
		 	case 'phasecoher',  formula = [ 'mean(arg1,3);' ]; inputdata = alltfX./sqrt(alltfX.*conj(alltfX));
		 	case 'phasecoher2', formula = [ 'sum(arg1,3)./sum(sqrt(arg1.*conj(arg1)),3);' ];
		end;
        if strcmpi(g.boottype, 'randall'), dimaccu = []; g.boottype = 'rand';
        else										 dimaccu = 2;
        end;
		Rboot = bootstat(inputdata, formula, 'boottype', g.boottype, ...
						   'label', 'ITC', 'bootside', 'upper', 'naccu', g.naccu, ...
                           'basevect', baselntmp, 'alpha', g.alpha, 'dimaccu', 2 );
	end;
else 
    Pboot = []; Rboot = [];
end

% --------
% plotting
% --------
ERP = mean(X,2);
mbase = log10(mbase)*10;
if strcmpi(g.plotersp, 'on') | strcmpi(g.plotitc, 'on')
    if strcmpi(g.plottype, 'image')
        plottimef(P, R, Pboot, Rboot, ERP, freqs, timesout, mbase, g);
    else    
        plotallcurves(P, R, Pboot, Rboot, ERP, freqs, timesout, mbase, g);
    end;
end;
if strcmpi(g.outputformat, 'old')
    R = abs(R); % convert coherence vector to magnitude
    mbase = 10^(mbase/10);
end;
return;

% -----------------
% plotting function
% -----------------
function g = plottimef(P, R, Pboot, Rboot, ERP, freqs, times, mbase, g);
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
    
	if ~isreal(R)
        Rangle = angle(R);
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
        myprintf(g.verbose, '\nNow plotting...\n');
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
      set(h(1), 'tag', 'ersp');
      
      PP = P;
      if ~isnan(g.alpha) % zero out nonsignif. power differences
          if size(PP,1) == size(Pboot,1) & size(PP,2) == size(Pboot,2)
              PP(find(PP > Pboot(:,:,1) & (PP < Pboot(:,:,2)))) = 0;
              Pboot = squeeze(mean(Pboot,2));
          else
              PP(find((PP > repmat(Pboot(:,1),[1 length(times)])) ...
                  & (PP < repmat(Pboot(:,2),[1 length(times)])))) = 0;
          end;
      end

      if isempty(g.erspmax)
          if g.ERSP_CAXIS_LIMIT == 0
              g.erspmax = [-1 1]*1.1*max(max(abs(P(:,:))));
          else
              g.erspmax = g.ERSP_CAXIS_LIMIT*[-1 1];
          end
      elseif length(g.erspmax) == 1
          g.erspmax = [ -g.erspmax g.erspmax];
      end;

      if ~strcmpi(g.freqscale, 'log')
          if ~isnan( g.baseline ) 
              imagesc(times,freqs,PP(:,:),g.erspmax); 
          else
              imagesc(times,freqs,PP(:,:));
          end;
      else 
          if ~isnan( g.baseline ) 
              imagesclogy(times,freqs,PP(:,:),g.erspmax); 
          else
              imagesclogy(times,freqs,PP(:,:));
          end;
      end;
          
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
          end;
      end;

      h(2) = gca;
      h(3) = cbar('vert'); % ERSP colorbar axes
      set(h(2),'Position',[.1 ordinate1 .8 height].*s+q)
      set(h(3),'Position',[.95 ordinate1 .05 height].*s+q)
      title('ERSP (dB)')

      %
      %%%%% plot marginal ERSP mean below ERSP image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %

      h(4) = subplot('Position',[.1 ordinate1-0.1 .8 .1].*s+q);
      
      E = [min(P(:,:),[],1);max(P(:,:),[],1)];

      % plotting limits
      if isempty(g.erspmarglim)
          g.erspmarglim = [min(E(1,:))-max(max(abs(E)))/3 max(E(2,:))+max(max(abs(E)))/3];
      end;
      
      plot(times,E,[0 0],g.erspmarglim, '--m','LineWidth',g.linewidth)
      xlim([min(times) max(times)])
      ylim(g.erspmarglim)
      
      tick = get(h(4),'YTick');
      set(h(4),'YTick',[tick(1) ; tick(end)])
      set(h(4),'YAxisLocation','right')
      set(h(4),'TickLength',[0.020 0.025]);
      xlabel('Time (ms)')
      ylabel('dB')

      %
      %%%%% plot mean spectrum to left of ERSP image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %

      h(5) = subplot('Position',[0 ordinate1 .1 height].*s+q);
      
      E = mbase;

      if ~isnan(E)
          
          % ploting limits
          if isempty(g.speclim)
              g.speclim = [min(E)-max(abs(E))/3 max(E)+max(abs(E))/3];
          end;
          
          % plot curves
          if ~strcmpi(g.freqscale, 'log')
              plot(freqs,E,'LineWidth',g.linewidth); hold on;
              if ~isnan(g.alpha)
                  plot(freqs,Pboot(:,:)'+[E;E], 'g', 'LineWidth',g.linewidth)
                  plot(freqs,Pboot(:,:)'+[E;E], 'k:','LineWidth',g.linewidth)
              end
              if freqs(1) ~= freqs(end), xlim([freqs(1) freqs(end)]); end;
              ylim(g.speclim)
          else
              semilogx(freqs,E,'LineWidth',g.linewidth); hold on;
              if ~isnan(g.alpha)
                  semilogx(freqs,Pboot(:,:)'+[E;E],'g', 'LineWidth',g.linewidth)
                  semilogx(freqs,Pboot(:,:)'+[E;E],'k:','LineWidth',g.linewidth)
              end
              if freqs(1) ~= freqs(end), xlim([freqs(1) freqs(end)]); end;
              ylim(g.speclim)
              set(h(5),'View',[90 90])
              divs = linspace(log(freqs(1)), log(freqs(end)), 10);
              set(gca, 'xtickmode', 'manual');
              divs = ceil(exp(divs)); divs = unique(divs); % ceil is critical here, round might misalign
                                                           % out-of border label with within border ticks
              set(gca, 'xtick', divs);
          end;
          set(h(5),'TickLength',[0.020 0.025]);          
          set(h(5),'View',[90 90])
          tick = get(h(5),'YTick');
          if (length(tick)>2)
              set(h(5),'YTick',[tick(1) ; tick(end-1)])
          end
          xlabel('Frequency (Hz)')
          ylabel('dB')
      end;
    end;

    switch lower(g.plotitc)
     case 'on'
      %
      %%%%%%%%%%%% Image the ITC %%%%%%%%%%%%%%%%%%
      %
      h(6) = subplot('Position',[.1 ordinate2 .9 height].*s+q); % ITC image
      set(h(1), 'tag', 'itc');

      if abs(R(1,1)-1) < 0.0001, g.plotphase = 'on'; end;
      if strcmpi(g.plotphase, 'on')
          RR = Rangle/pi*180;
      else
          RR = R;
      end;
      if ~isnan(g.alpha)
          if size(RR,1) == size(Rboot,1) & size(RR,2) == size(Rboot,2)
              tmp = gcf;
			  if size(Rboot,3) == 2	 RR(find(RR > Rboot(:,:,1) & RR < Rboot(:,:,2))) = 0;			  
			  else                   RR(find(RR < Rboot)) = 0;				  
			  end;
			  Rboot = mean(Rboot(:,:,end),2);
		  else
              RR(find(RR < repmat(Rboot(:),[1 length(times)]))) = 0;
          end;
      end

      if g.ITC_CAXIS_LIMIT == 0
          coh_caxis = min(max(max(R(:,:))),1)*[-1 1]; % 1 WAS 0.4 !
      else
          coh_caxis = g.ITC_CAXIS_LIMIT*[-1 1];
      end

      if strcmpi(g.plotphase, 'on')
        if ~strcmpi(g.freqscale, 'log')
              imagesc(times,freqs,RR(:,:)); % <---
        else 
              imagesclogy(times,freqs,RR(:,:)); % <---
        end;
        g.itcmax = [-180 180];
        setylim = 0;
      else          
        if ~strcmpi(g.freqscale, 'log')
          if exist('Rsign') & strcmp(g.plotphasesign, 'on')
              imagesc(times,freqs,Rsign(:,:).*RR(:,:),coh_caxis); % <---
          else
              imagesc(times,freqs,RR(:,:),coh_caxis); % <---
          end
        else 
          if exist('Rsign') & strcmp(g.plotphasesign, 'on')
              imagesclogy(times,freqs,Rsign(:,:).*RR(:,:),coh_caxis); % <---
          else
              imagesclogy(times,freqs,RR(:,:),coh_caxis); % <---
          end
        end;
      end;

      if isempty(g.itcmax)
          g.itcmax = caxis;
      elseif length(g.itcmax) == 1
          g.itcmax = [ -g.itcmax g.itcmax ];
      end;
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
          end;
      end;

      h(7) = gca;
      h(8) = cbar('vert');
      %h(9) = get(h(8),'Children'); % make the function crash
      set(h(7),'Position',[.1 ordinate2 .8 height].*s+q)
      set(h(8),'Position',[.95 ordinate2 .05 height].*s+q)
	  if setylim
		  set(h(8),'YLim',[0 g.itcmax(2)]); 
      end;
      if strcmpi(g.plotphase, 'on')
	    title('ITC phase')
      else  
	    title('ITC')
      end;

      %
      %%%%% plot the ERP below the ITC image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %

      h(10) = subplot('Position',[.1 ordinate2-0.1 .8 .1].*s+q); % ERP

      if isempty(g.erplim)
          ERPmax = max(ERP);
          ERPmin = min(ERP);
          g.erplim = [ ERPmin - 0.1*(ERPmax-ERPmin) ERPmax + 0.1*(ERPmax-ERPmin) ];
      end;

      plot(ERPtimes,ERP, [0 0],g.erplim,'--m','LineWidth',g.linewidth);
      hold on;
      plot([times(1) times(length(times))],[0 0], 'k');
      xlim([min(ERPtimes) max(ERPtimes)]);
      ylim(g.erplim)

      tick = get(h(10),'YTick');
      set(h(10),'YTick',[tick(1) ; tick(end)])
      set(h(10),'TickLength',[0.02 0.025]);
      set(h(10),'YAxisLocation','right')
      xlabel('Time (ms)')
      ylabel('uV')

      E = mean(R(:,:)');
      
      %
      %%%%% plot the marginal mean ERP below the ITC image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %

      h(11) = subplot('Position',[0 ordinate2 .1 height].*s+q); % plot the marginal mean
                                                                % ITC left of the ITC image
      
      % set plotting limits
      if isempty(g.itcavglim)
          if ~isnan(g.alpha)
              g.itcavglim = [ min(E)-max(E)/3 max(Rboot)+max(Rboot)/3];
          else
              g.itcavglim = [ min(E)-max(E)/3 max(E)+max(E)/3];
          end;
      end;
      
      % plot marginal ITC
      if ~strcmpi(g.freqscale, 'log')
          plot(freqs,E,'LineWidth',g.linewidth); hold on;
          if ~isnan(g.alpha)
              plot(freqs,Rboot,'g', 'LineWidth',g.linewidth)
              plot(freqs,Rboot,'k:','LineWidth',g.linewidth)
          end
          if freqs(1) ~= freqs(end), xlim([freqs(1) freqs(end)]); end;
          ylim(g.itcavglim)
      else
          semilogx(freqs,E,'LineWidth',g.linewidth); hold on;
          if ~isnan(g.alpha)
              semilogx(freqs,Rboot(:),'g', 'LineWidth',g.linewidth)
              semilogx(freqs,Rboot(:),'k:','LineWidth',g.linewidth)
          end
          if freqs(1) ~= freqs(end), xlim([freqs(1) freqs(end)]); end;
          ylim(g.itcavglim)
          divs = linspace(log(freqs(1)), log(freqs(end)), 10);
          set(gca, 'xtickmode', 'manual');
          divs = ceil(exp(divs)); divs = unique(divs); % ceil is critical here, round might misalign
                                                       % out-of border label with within border ticks
          set(gca, 'xtick', divs);
      end;

      % plot details
      tick = get(h(11),'YTick');
      if length(tick) > 1
          set(h(11),'YTick',[tick(1) ; tick(length(tick))])
      end;
      set(h(11),'View',[90 90])
      %set(h(11),'TickLength',[0.020 0.025]);
      xlabel('Frequency (Hz)')
      ylabel('ERP')
      
      %
      %%%%%%%%%%%%%%% plot a topoplot() %%%%%%%%%%%%%%%%%%%%%%%
      %
      if (~isempty(g.topovec))
          h(12) = subplot('Position',[-.1 .43 .2 .14].*s+q);
          if length(g.topovec) == 1
              topoplot(g.topovec,g.elocs,'electrodes','off', ...
                       'style', 'blank', 'emarkersize1chan', 10, 'chaninfo', g.chaninfo);
          else
              topoplot(g.topovec,g.elocs,'electrodes','off', 'chaninfo', g.chaninfo);
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

% ---------------
% Plotting curves
% ---------------
function plotallcurves(P, R, Pboot, Rboot, ERP, freqs, times, mbase, g);

	if ~isreal(R)
		Rangle = angle(R);
		R = abs(R); % convert coherence vector to magnitude
		setylim = 1;
    else 
		Rsign = ones(size(R));
		setylim = 0;
	end;

    if strcmpi(g.plotitc, 'on') | strcmpi(g.plotersp, 'on')
        myprintf(g.verbose, '\nNow plotting...\n');
        pos = get(gca,'position');
        q = [pos(1) pos(2) 0 0];
        s = [pos(3) pos(4) pos(3) pos(4)];
    end;

    % time unit
    % ---------
    if times(end) > 10000
        times = times/1000;
        timeunit = 's';
    else
        timeunit = 'ms';
    end;
    
    if strcmpi(g.plotersp, 'on')
      %
      %%%%%%% image the ERSP %%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      if strcmpi(g.plotitc, 'on'), subplot(2,1,1); end; 
      set(gca, 'tag', 'ersp');
      alllegend = {};
      if strcmpi(g.plotmean, 'on') & freqs(1) ~= freqs(end)
          alllegend = { [ num2str(freqs(1)) '-' num2str(freqs(end)) ...
                      'Hz mean baseline ' num2str(mean(mbase)) ' dB' ] };
      else
          for index = 1:length(freqs)
             alllegend{index} = [ num2str(freqs(index)) 'Hz baseline ' num2str(mbase(index)) ' dB' ];
          end;
      end;
      plotcurve(times, P, 'maskarray', Pboot, 'title', 'ERSP', ...
                'xlabel', [ 'Time (' timeunit ')' ], 'ylabel', 'dB', 'ylim', [-g.erspmax g.erspmax], ...
                'vert', g.vert, 'marktimes', g.marktimes, 'legend', alllegend, ...
                'linewidth', g.linewidth, 'highlightmode', g.highlightmode, 'plotmean', g.plotmean);
    end;

    if strcmpi(g.plotitc, 'on')
      %
      %%%%%%%%%%%% Image the ITC %%%%%%%%%%%%%%%%%%
      %
      if strcmpi(g.plotersp, 'on'), subplot(2,1,2); end; 
      set(gca, 'tag', 'itc');
      if abs(R(1,1)-1) < 0.0001, g.plotphase = 'on'; end;
      if strcmpi(g.plotphase, 'on') % plot ITC phase instead of amplitude (e.g. for continuous data)
           RR = Rangle/pi*180;
      else RR = R;
      end; 
      
      % find regions of significance
      % ----------------------------
      alllegend = {};
      if strcmpi(g.plotmean, 'on') & freqs(1) ~= freqs(end)
          alllegend = { [ num2str(freqs(1)) '-' num2str(freqs(end)) 'Hz' ] };
      else 
          for index = 1:length(freqs)
              alllegend{index} = [ num2str(freqs(index)) 'Hz' ];
          end;
      end;
      plotcurve(times, RR, 'maskarray', Rboot, 'val2mask', R, 'title', 'ITC', ...
                'xlabel', [ 'Time (' timeunit ')' ], 'ylabel', 'dB', 'ylim', g.itcmax, ...
                'vert', g.vert, 'marktimes', g.marktimes, 'legend', alllegend, ...
                'linewidth', g.linewidth, 'highlightmode', g.highlightmode, 'plotmean', g.plotmean);
    end;

    if strcmpi(g.plotitc, 'on') | strcmpi(g.plotersp, 'on')
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
end;

if ~isempty(regions)
    axes(ax);
    in_a_region = 0;
    for index=1:length(regions)
        if regions(index) & ~in_a_region
            tmpreg(1) = times(index);
            in_a_region = 1;
        end;
        if ~regions(index) & in_a_region
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
            end;
        end;
    end;
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

function myprintf(verbose, varargin)
    if strcmpi(verbose, 'on')
        fprintf(varargin{:});
    end;
