% erpimage() - Image a collection of single-trial data epochs, optionally sorted on 
%              and/or aligned to an input sorting variable and smoothed across trials 
%              with a moving-average.  (To return event-aligned data without plotting, 
%              use eventlock()).  Optionally sort trials on value, amplitude or phase 
%              within a specified latency window. Optionally plot the ERP mean and std. dev.,
%              and moving-window spectral amplitude and inter-trial coherence at a
%              selected or peak frequency. Click on individual figures parts to examine 
%              them separately and zoom (using axcopy()).
% Usage:
%   >> [outdata,outvar,outtrials,limits,axhndls,erp, ...
%         amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx,erpsig] ...
%                  = erpimage(data,sortvar,times,'title',avewidth,decimate,...
%                             flag1,arg1,flag2,arg2,...);
% Necessary inputs:
%   data     - [vector or matrix] Single-channel input data to image. 
%               Formats (1,frames*trials) or (frames,trials)
%
% Optional ordered inputs {with defaults}:
%   sortvar  - [vector | []] Variable to sort epochs on (length(sortvar) = nepochs)
%              Example: sortvar may by subject response time in each epoch (in ms)
%              {default|[]: plot in input order}
%   times    - [vector | []] of latencies (ms) (length(times) = frames) 
%               ELSE [startms ntimes srate] Give start latency (ms), time points 
%               (i.e. frames) per epoch, sampling rate (Hz), {default|[]: 0:nframes-1}
%  'title'   - ['string'] Plot titla {default: none}
%   avewidth - Number of trials to smooth with a moving-average (may be non-integer) 
%               {default|0->1}
%   decimate - Factor to decimate ntrials out by (may be non-integer)  {default|0->1}
%               If this is large ( > sqrt(num. trials)), output this many trials.
%
% Unordered options ('keyword',argument pairs):
%
% Optionally realign data epochs: 
%   'align'  - [latency] Time-lock data to sortvar. Plot sortvar as at latency (ms)
%               If latency == Inf, plot at sortvar median {default: no align}
%   'renorm' - ['yes'|'no'|'formula(x)'] Normalize sorting variable to epoch 
%               latency range and plot. 'yes'= autoscale. Example of formula(x):
%               '3*x+2'. {default: 'no'}
%   'noplot' - Do not plot sortvar {default: Do plot sortvar if in times range}
%   'noshow' - ['yes'|'no'] Do not plot erpimage, simply return outputs {default: 'no'}
%
% Optionally sort input epochs: 
%   'nosort' - Do not sort data epochs.
%  {default} - Sort data epochs by sortvar (see Necessary inputs above).
%  'valsort' - [startms endms direction] Sort data on (mean) value 
%               between startms and (optional) endms. Direction is 1 or -1.
%              If -1, plot max-value epoch at bottom {Default: sort on sortvar}
% 'phasesort' - [ms_center prct freq maxfreq topphase] Sort epochs by phase in 
%                an n-cycle window centered at latency ms_center (ms). 
%                Percentile (prct) in range [0,100] gives percent of trials 
%                to reject for low amplitude. Else, if in range [-100,0], 
%                percent of trials to reject for high amplitude; freq (Hz) 
%                is the phase-sorting frequency. With optional maxfreq,
%                sort by phase at freq of max power in the data in range 
%                [freq,maxfreq] (Note: 'phasesort' arg freq overrides the 
%                frequency specified in 'coher'). With optional topphase, 
%                sort by phase, putting topphase (degrees, in range [-180,180]) 
%                at the top of the image. Note: 'phasesort' now uses circular 
%                smoothing. Use 'cycles' (below) for wavelet length. 
%                {Default: [0 25 8 13 180]}
%  'ampsort' - [center_ms prcnt freq maxfreq] Sort epochs by amplitude. 
%                See 'phasesort'. If ms_center is 'Inf', then sorting
%                is by mean power across the time window specified by 'winsort' below.
%                If arg 'freq' above is <0, sort by mean power in the range [freq maxfreq].
%  'sortwin' - [start_ms end_ms] With center_ms == Inf in 'ampsort' ars (above), sorts
%                by mean amplitude across window centers shifted from start_ms to end_ms by 10 ms.
%  'showwin' - Show sorting window behind ERP trace. {default: don't show sorting window}
%
% Plot time-varying spectral amplitude instead of potential:
% 'plotamps' - Image amplitudes at each trial and latency instead of potential values. 
%               Note: Currently requires 'coher' (below) with alpha signif. {default: no}
%
% Specify plot parameters:
%   'limits' - [lotime hitime minerp maxerp loamp hiamp locoher hicoher bamp]
%               Plot axes limits. Can use NaN (or nan, but not Nan) for missing items 
%               and omit late items. Use last input, bamp, to fix the baseline amplitude.
%   'signif' - [lo_amp, hi_amp, coher_signif_level] Use preassigned significance 
%               levels to save computation time. {default: none}
%   'caxis'  - [lo hi] Set color axis limits ELSE [fraction] Set caxis limits at 
%               (+/-)fraction*max(abs(data)) {default: symmetrical, based on data limits}
%
% Add epoch-mean ERP to plot:
%   'erp'    - Plot ERP time average of the trials below the image {default no ERP plotted}
%   'erpalpha' - [alpha] One-sided significance threshold (range: [.001 0.1]). 
%              Requires 'erp' {default no +/-alpha significance thresholds plotted}
%   'erpstd' - Plot ERP +/- stdev. Requires 'erp' {default no +/-stdev plotted}
%   'rmerp'  - Subtract the average ERP from each trial before processing {default no}
%
% Add time/frequency information:
%  'coher'   - [freq] Plot ERP average plus mean amplitude & coherence at freq (Hz)
%               ELSE [minfrq maxfrq] = same, but select frequency with max power in 
%               given range (Note: 'phasesort' freq (above) overwrites these parameters).
%               ELSE [minfrq maxfrq alpha] = plot coher. signif. level line at 
%               probability alpha (range: [0,0.1]) {default: no coher, no probabilities}
%   'srate'  - [freq] Specify the data sampling rate in Hz for amp/coher (if not 
%               implicit in third arg times) {default: as in icadefs.m}
%   'cycles' - Number of cycles in the wavelet time/frequency decomposition {default: 3}
%
% Add other features:
%   'cbar'   - Plot color bar to right of ERP-image {default no}
%   'topo'   - {map_vals,eloc_file} Plot a 2-D scalp map at upper left of image. 
%               See '>> topoplot example' for electrode location file structure.
%   'spec'   - [loHz,hiHz] Plot the mean data spectrum at upper right of image. 
%   'horz'   - [epochs_vector] Plot horizontal lines at specified epochs
%   'vert'   - [times_vector] Plot vertical dashed lines at specified latencies
%   'auxvar' - [(nvars,ntrials) matrix] Plot auxiliary variable(s) for each trial 
%               as separate traces. ELSE, 'auxvar',{[that_matrix],{colorstrings}} 
%               to specify N trace colors.  Ex: colorstrings = {'r','bo-','k:'} 
%               (See also: 'vert' above).
% Miscellaneous options:
% 'noxlabel' - Do not plot "Time (ms)" on the bottom x-axis
% 'yerplabel' - ['string'] ERP ordinate axis label (default is ERP). Get uV with '\muV'
%
% Optional outputs:
%    outdata  = (times,epochsout) data matrix (after smoothing)
%     outvar  = (1,epochsout) actual values trials are sorted on (after smoothing)
%   outtrials = (1,epochsout)  smoothed trial numbers
%     limits  = (1,10) array, 1-9 as in 'limits' above, then analysis frequency (Hz) 
%    axhndls  = vector of 1-7 plot axes handles (img,cbar,erp,amp,coh,topo,spec)
%        erp  = plotted ERP average
%       amps  = mean amplitude time course
%      coher  = mean inter-trial phase coherence time course
%     cohsig  = coherence significance level
%     ampsig  = amplitude significance levels [lo high]
%    outamps  = matrix of imaged amplitudes (from option 'allamps')
%   phsangls  = vector of sorted trial phases at the phase-sorting frequency
%     phsamp  = vector of sorted trial amplitudes at the phase-sorting frequency
%    sortidx  = indices of sorted data epochs plotted
%     erpsig  = trial average significance levels [2,frames]
%
% Example:  >> figure; erpimage(data,RTs,[-400 256 256],'Test',1,1,'erp','cbar','vert',-350);
%
% Plots an ERP-image of 1-s data epochs sampled at 256 Hz, sorted by RTs, title 'Test', 
% sorted epochs not smoothed or decimated. Also plots the epoch-mean ERP, a color bar, 
% and a dashed vertical line at -350 ms.

% Authors: Scott Makeig, Tzyy-Ping Jung & Arnaud Delorme, 
%          CNL/Salk Institute, La Jolla, 3-2-1998 -
%
% See also: erpimages(), phasecoher(), rmbase(), cbar(), movav()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Scott Makeig & Tzyy-Ping Jung, CNL / Salk Institute, La Jolla 3-2-98
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

% Uses external toolbox functions: phasecoher(), rmbase(), cbar(), movav()
% Uses included functions:         plot1trace(), phasedet()

% UNIMPLEMENTED - 'allcohers',[data2] -> image the coherences at each latency & epoch. 
%                 Requires arg 'coher' with alpha significance. 
%                 Shows projection on grand mean coherence vector at each latency 
%                 and trial. {default: no}
 
% $Log: not supported by cvs2svn $
% Revision 1.220  2004/09/02 00:04:00  scott
% added 'sortwin' and ampsort by mean freq in a range. -sm
%
% Revision 1.219  2004/09/01 21:56:28  scott
% adding 'sortwin' for 'ampsort' -sm
%
% Revision 1.218  2004/08/30 18:32:06  scott
% added figure(curfig) before plotting functions to keep
% version 7.0.0 from reverting plotting to the EEGLAB menu window
% when plotting from the gui (only).  -sm
%
% Revision 1.217  2004/08/30 17:08:49  scott
% debug last
%
% Revision 1.216  2004/08/30 17:05:17  scott
% reshape verttimes if necessary
%
% Revision 1.215  2004/08/20 00:39:04  arno
% debut phasearg if only using 1 frequency
%
% Revision 1.214  2004/08/13 20:27:12  scott
% help msg
%
% Revision 1.213  2004/07/29 23:18:23  arno
% same
%
% Revision 1.212  2004/07/29 23:16:50  arno
% convert to double for Matlab 7
%
% Revision 1.211  2004/06/09 01:48:43  arno
% make limits and image data consistent
%
% Revision 1.210  2004/06/05 01:52:45  arno
% aligntime as string to avoid really realigning
%
% Revision 1.209  2004/06/05 01:33:07  arno
% spetial option to preserve backward compatibility
%
% Revision 1.208  2004/05/07 04:47:17  scott
% made sotvar = [] work
%
% Revision 1.207  2004/03/26 00:21:06  arno
% same
%
% Revision 1.206  2004/03/26 00:19:48  arno
% minimum number of trials
%
% Revision 1.205  2004/03/26 00:18:11  arno
% plot aligntime
%
% Revision 1.204  2004/02/24 23:04:40  arno
% fixed vertical lines in ERP when RT-aligned
%
% Revision 1.203  2004/01/24 22:01:23  scott
% *** empty log message ***
%
% Revision 1.202  2004/01/24 21:58:33  scott
% same
%
% Revision 1.201  2004/01/24 21:53:38  scott
% same
%
% Revision 1.200  2004/01/24 21:43:59  scott
% *** empty log message ***
%
% Revision 1.199  2004/01/24 21:36:01  scott
% *** empty log message ***
%
% Revision 1.198  2004/01/24 21:26:50  scott
% same
%
% Revision 1.197  2004/01/24 21:25:04  scott
% same
%
% Revision 1.196  2004/01/24 21:22:59  scott
% same
%
% Revision 1.195  2004/01/24 21:17:00  scott
% same
%
% Revision 1.194  2004/01/24 21:16:32  scott
% same
%
% Revision 1.193  2004/01/24 21:10:31  scott
% same
%
% Revision 1.192  2004/01/24 21:08:29  scott
% same
%
% Revision 1.191  2004/01/24 21:04:09  scott
% same
%
% Revision 1.190  2004/01/24 21:01:18  scott
% same
%
% Revision 1.189  2004/01/24 20:59:53  scott
% same
%
% Revision 1.188  2004/01/24 20:51:57  scott
% same
%
% Revision 1.187  2004/01/24 20:45:11  scott
% same
%
% Revision 1.186  2004/01/24 20:42:05  scott
% same
%
% Revision 1.185  2004/01/24 20:40:23  scott
% plotting sorting window
%
% Revision 1.184  2003/12/17 21:42:13  scott
% adjust same
%
% Revision 1.183  2003/12/17 21:40:12  scott
% change y-labels on traces
%
% Revision 1.182  2003/12/04 17:53:52  arno
% debug spec()
%
% Revision 1.181  2003/12/03 02:18:48  arno
% use spec if psd is absent
%
% Revision 1.180  2003/11/26 18:16:23  scott
% help msg
%
% Revision 1.179  2003/11/19 01:06:22  arno
% including makehanning function
%
% Revision 1.178  2003/11/16 17:48:00  scott
% plot1erp() -> plot1trace(); printf "Done."
%
% Revision 1.177  2003/11/14 17:01:59  scott
% bootstrap msg
%
% Revision 1.176  2003/11/14 16:58:02  scott
% debug last
%
% Revision 1.175  2003/11/14 16:56:16  scott
% refining erpsig, sig fills
%
% Revision 1.174  2003/11/14 16:44:06  scott
% changed signif fill color
%
% Revision 1.173  2003/11/14 16:34:31  scott
% debug same
%
% Revision 1.172  2003/11/14 16:32:55  scott
% debug same
%
% Revision 1.171  2003/11/14 16:27:27  scott
% fill coher signif limits
%
% Revision 1.170  2003/11/13 02:32:42  scott
% fill the dB signif limits
%
% Revision 1.169  2003/11/13 02:29:51  scott
% debug
%
% Revision 1.168  2003/11/13 02:28:12  scott
% debug
%
% Revision 1.167  2003/11/13 02:26:25  scott
% debug
%
% Revision 1.166  2003/11/13 02:25:29  scott
% debug
%
% Revision 1.165  2003/11/13 02:24:19  scott
% fill amp signif
%
% Revision 1.164  2003/11/13 02:14:16  scott
% same
%
% Revision 1.163  2003/11/13 01:58:10  scott
% same
%
% Revision 1.162  2003/11/13 01:47:13  scott
% make erpalpha fill less saturated
%
% Revision 1.161  2003/11/10 23:40:32  scott
% fill the bootstrap limits behind the erp if erpalpha
%
% Revision 1.160  2003/11/10 16:57:46  arno
% msg typo
%
% Revision 1.159  2003/10/29 22:15:52  scott
% adjust erp alpha NBOOT to (low) alpha
%
% Revision 1.158  2003/10/29 22:07:10  arno
% nothing
%
% Revision 1.157  2003/09/24 19:36:19  scott
% fixed auxvar bug
%
% Revision 1.156  2003/09/24 18:56:22  scott
% debug
%
% Revision 1.155  2003/09/24 00:45:30  scott
% debug same
%
% Revision 1.154  2003/09/24 00:43:15  scott
% debug same
%
% Revision 1.153  2003/09/24 00:42:10  scott
% debug same
%
% Revision 1.152  2003/09/24 00:39:05  scott
% adding 'horz' -> horizontal line plotting
%
% Revision 1.151  2003/09/21 21:12:27  scott
% edited comments
%
% Revision 1.150  2003/09/11 22:23:25  scott
% debug same
%
% Revision 1.149  2003/09/11 22:19:49  scott
% made 'plotamps' and 'auxvar' work together
%
% Revision 1.148  2003/09/09 23:26:48  arno
% change && to &
%
% Revision 1.147  2003/09/07 00:41:51  arno
% fixing stdev, do not know why it crashed
%
% Revision 1.146  2003/09/06 22:45:04  scott
% same
%
% Revision 1.145  2003/09/06 22:43:57  scott
% add erpsig output
%
% Revision 1.144  2003/09/06 22:24:55  scott
% debug last, add \n before first printed line "Plotting input...
%
% Revision 1.143  2003/09/06 22:15:23  scott
% adjust auxvar if phase sort or amp sort
%
% Revision 1.142  2003/08/27 17:45:07  scott
% header
%
% Revision 1.141  2003/08/25 22:38:25  scott
% auxvars adjust tests -sm
%
% Revision 1.140  2003/08/24 04:49:01  scott
% help msg
%
% Revision 1.139  2003/08/24 04:38:18  scott
% same
%
% Revision 1.138  2003/08/24 04:37:47  scott
% debug last
%
% Revision 1.137  2003/08/24 04:36:29  scott
% fprintf
%
% Revision 1.136  2003/08/24 04:35:19  scott
% added help for 'erpalpha' -sm
%
% Revision 1.135  2003/08/24 04:27:41  scott
% fprintf adjust
%
% Revision 1.134  2003/08/24 04:25:05  scott
% adjust erpalpha line type
%
% Revision 1.133  2003/08/24 04:23:40  scott
% debug
%
% Revision 1.132  2003/08/24 04:20:26  scott
% same
%
% Revision 1.131  2003/08/24 04:19:58  scott
% added fprintf
%
% Revision 1.130  2003/08/24 04:17:42  scott
% same
%
% Revision 1.129  2003/08/24 04:17:17  scott
% same
%
% Revision 1.128  2003/08/24 04:16:39  scott
% debug same
%
% Revision 1.127  2003/08/24 04:12:33  scott
% added (erpalpha) bootstrap ERP estimation -sm
%
% Revision 1.126  2003/07/26 17:18:51  scott
% debug same
%
% Revision 1.125  2003/07/26 17:16:50  scott
% debug same
%
% Revision 1.124  2003/07/26 17:14:12  scott
% debug same
%
% Revision 1.123  2003/07/26 17:05:26  scott
% debug auxvars warning, hide vert matrix option
%
% Revision 1.122  2003/07/24 23:41:05  arno
% correcting output
%
% Revision 1.121  2003/07/24 18:17:33  scott
% *** empty log message ***
%
% Revision 1.120  2003/07/22 16:03:06  scott
% phasesort help message
%
% Revision 1.119  2003/07/22 15:40:26  scott
% adding circular smoothing to phase-sorted, allamps plots -sm
%
% Revision 1.118  2003/07/22 00:52:02  scott
% debug
%
% Revision 1.117  2003/07/22 00:48:36  scott
% debug
%
% Revision 1.116  2003/07/22 00:44:21  scott
% debug
%
% Revision 1.115  2003/07/21 21:38:32  scott
% debug
%
% Revision 1.114  2003/07/21 21:36:36  scott
% debug
%
% Revision 1.113  2003/07/21 21:33:48  scott
% debug
%
% Revision 1.112  2003/07/21 21:25:44  scott
% debug
%
% Revision 1.111  2003/07/21 21:19:15  scott
% debug
%
% Revision 1.110  2003/07/21 21:17:52  scott
% debug
%
% Revision 1.109  2003/07/21 21:13:09  scott
% debug
%
% Revision 1.108  2003/07/21 21:09:59  scott
% debug
%
% Revision 1.107  2003/07/21 20:46:19  scott
% debug
%
% Revision 1.106  2003/07/21 20:45:01  scott
% fg
%
% debug last
%
% Revision 1.105  2003/07/21 20:43:15  scott
% using wraparound smoothing for phase-sorted epochs
%
% Revision 1.104  2003/07/15 18:01:05  scott
% trial number -> trials throughout
%
% Revision 1.103  2003/05/06 17:31:29  arno
% debug Rmerp
%
% Revision 1.102  2003/05/06 15:58:51  scott
% header edits
%
% Revision 1.101  2003/05/06 15:54:15  scott
% edit header
%
% Revision 1.100  2003/05/06 00:49:52  arno
% debug last
%
% Revision 1.99  2003/05/06 00:45:43  arno
% implementing new option rmerp
%
% Revision 1.98  2003/04/26 01:06:38  arno
% debuging ampsort for phase sorting
%
% Revision 1.97  2003/04/25 22:32:26  arno
% doing the same for ampsort
%
% Revision 1.96  2003/04/25 22:30:40  arno
% interpolating phase value
%
% Revision 1.95  2003/04/25 18:04:31  arno
% nothing
%
% Revision 1.94  2003/04/24 21:54:24  arno
% removing debug message
%
% Revision 1.93  2003/04/23 23:59:48  arno
% remove debug msg
%
% Revision 1.92  2003/04/23 23:58:02  arno
% debuging phasedet in erpimage
%
% Revision 1.91  2003/04/23 23:41:19  arno
% restoring cycles to a default of 3, adding cycle paramete
% r
%
% Revision 1.90  2003/04/23 22:09:21  arno
% adding cycles to the phasedet script
%
% Revision 1.89  2003/04/23 01:24:15  arno
% chaning default to 3 cycles at 5 Hz
%
% Revision 1.88  2003/04/23 01:18:08  arno
% typo
%
% Revision 1.87  2003/04/23 01:15:51  arno
% debuging cycles
%
% Revision 1.86  2003/04/23 01:10:39  arno
% remove 3-cycle from help
%
% Revision 1.85  2003/04/23 01:05:55  arno
% changing default cycle for coherence
%
% Revision 1.84  2003/04/18 17:20:13  arno
% noshow option
%
% Revision 1.83  2003/03/23 22:46:39  scott
% *** empty log message ***
%
% Revision 1.82  2003/03/15 18:01:34  scott
% help msg
%
% Revision 1.81  2003/03/13 03:20:33  scott
% restoring
%
% Revision 1.80  2003/03/07 22:21:46  scott
% same
%
% Revision 1.79  2003/03/07 22:15:19  scott
% same
%
% Revision 1.78  2003/03/07 22:13:51  scott
% same
%
% Revision 1.77  2003/03/07 22:11:58  scott
% debugging axcopy call
%
% Revision 1.76  2003/03/07 22:09:47  scott
% typo
%
% Revision 1.75  2003/03/07 22:08:30  scott
% axcopy - adding ''''
%
% Revision 1.74  2003/03/07 21:18:13  scott
% testing axcopy
%
% Revision 1.73  2003/03/07 21:15:25  scott
% adding eegplot option if click on erpimage -sm
%
% Revision 1.72  2003/03/04 20:29:18  arno
% header typo
%
% Revision 1.71  2003/01/31 23:23:34  arno
% adding erpstd for ploting standard deviation of ERP
%
% Revision 1.70  2003/01/02 16:55:02  scott
% allowed sortvars to be negative -sm
%
% Revision 1.69  2002/11/19 23:27:55  arno
% debugging spectrum for very long epochs
%
% Revision 1.68  2002/11/09 21:19:57  scott
% added example
%
% Revision 1.67  2002/11/09 21:07:29  scott
% reorganized help message optional argument list
%
% Revision 1.66  2002/10/15 17:52:00  scott
% help msg - only two args required
%
% Revision 1.65  2002/10/14 16:12:38  scott
% fprintf ms
%
% Revision 1.64  2002/10/14 16:09:53  scott
% valsort fprintf
%
% Revision 1.63  2002/10/14 16:07:29  scott
% removed moreargs - consolidated help message -sm
%
% Revision 1.62  2002/10/14 14:56:45  scott
% working on ampsort
%
% Revision 1.61  2002/10/14 00:42:58  scott
% added valsort direction and help msg -sm
%
% Revision 1.60  2002/10/13 23:56:54  scott
% *** empty log message ***
%
% Revision 1.59  2002/10/13 23:55:47  scott
% debugging valsort
%
% Revision 1.58  2002/10/13 23:51:22  scott
% edit valsort fprint
%
% Revision 1.57  2002/10/13 23:49:43  scott
% *** empty log message ***
%
% Revision 1.56  2002/10/13 23:48:30  scott
% *** empty log message ***
%
% Revision 1.55  2002/10/13 23:48:04  scott
% *** empty log message ***
%
% Revision 1.54  2002/10/13 23:46:37  scott
% *** empty log message ***
%
% Revision 1.53  2002/10/13 23:43:47  scott
% debugging
%
% Revision 1.52  2002/10/13 23:41:57  scott
% debug ampargs
%
% Revision 1.51  2002/10/13 23:41:11  scott
% debug ampargs
%
% Revision 1.50  2002/10/13 23:37:01  scott
% debug valsort
%
% Revision 1.49  2002/10/13 23:35:51  scott
% debug
%
% Revision 1.48  2002/10/13 23:33:10  scott
% valsort debug
%
% Revision 1.47  2002/10/13 23:28:18  scott
% added Ampflag, Valflag defaults -sm
%
% Revision 1.46  2002/10/13 23:24:48  scott
% typo
%
% Revision 1.45  2002/10/13 23:23:32  scott
% added ampsort, valsort args -sm
%
% Revision 1.44  2002/09/03 21:35:27  arno
% removing oridata completelly
%
% Revision 1.43  2002/09/03 21:15:23  scott
% go back to averaging urdata instead of oridata -sm
%
% Revision 1.42  2002/08/31 17:00:50  arno
% add yerplabel option
%
% Revision 1.41  2002/08/30 18:18:09  arno
% same
%
% Revision 1.40  2002/08/30 18:14:08  arno
% same
%
% Revision 1.39  2002/08/30 18:12:12  arno
% same
%
% Revision 1.38  2002/08/30 18:09:15  arno
% same
%
% Revision 1.37  2002/08/30 18:07:20  arno
% same.
%
% Revision 1.36  2002/08/30 18:01:11  arno
% same
%
% Revision 1.35  2002/08/30 18:00:22  arno
% debug erp axis and average
%
% Revision 1.34  2002/08/21 18:36:13  arno
% adding error message
%
% Revision 1.33  2002/08/19 19:48:25  arno
% commenting crosscoher for Mac compatibility
%
% Revision 1.32  2002/08/12 01:37:00  arno
% color
%
% Revision 1.31  2002/08/11 22:34:04  arno
% color
%
% Revision 1.30  2002/08/09 16:28:07  arno
% debugging allamps
%
% Revision 1.29  2002/08/05 18:04:55  arno
% performing moving average on allamps amplitude (and not log)
%
% Revision 1.28  2002/07/27 01:25:45  arno
% debugging vert
%
% Revision 1.27  2002/07/26 16:18:32  arno
% removing debugging messages
%
% Revision 1.26  2002/07/26 16:14:03  arno
% removing trials with Nan values for sortvar, debugging 'vert'
%
% Revision 1.25  2002/07/15 02:00:48  arno
% same
%
% Revision 1.24  2002/07/15 01:55:45  arno
% debugging minamp maxamp
%
% Revision 1.23  2002/07/15 01:48:20  arno
% force yscale on amp ploterp
%
% Revision 1.22  2002/07/14 03:17:57  arno
% making allamps moving average of log power
%
% Revision 1.21  2002/07/14 02:39:16  arno
% same
%
% Revision 1.20  2002/07/14 02:23:16  arno
% same
%
% Revision 1.19  2002/07/14 01:58:00  arno
% same
%
% Revision 1.18  2002/07/14 01:47:23  arno
% testing amps limits
%
% Revision 1.17  2002/05/23 16:59:09  scott
% replaced nanmean with nan_mean() -sm
%
% Revision 1.16  2002/05/22 06:06:28  marissa
% changed line 1047 to remove nonexistent variable 'baseall'
%
% Revision 1.15  2002/05/20 17:58:53  scott
% adding fprintf info about ampsig plotting -sm
%
% Revision 1.14  2002/04/25 17:52:51  arno
% removing debug message
%
% Revision 1.13  2002/04/25 17:51:50  arno
% including renorm inside function
%
% Revision 1.12  2002/04/25 17:08:03  arno
% correcting error message
%
% Revision 1.11  2002/04/24 19:04:45  arno
% further check for coherfreq
%
% Revision 1.10  2002/04/24 18:29:21  arno
% two inputs for coher and time centering for phase
%
% Revision 1.9  2002/04/24 17:35:50  arno
% rechanged log
%
% 3/5/98 added nosort option -sm
% 3/22/98 added colorbar ylabel, sym. range finding -sm
% 5/08/98 added noplot option -sm
% 6/09/98 added align, erp, coher options -sm
% 6/10/98 added limits options -sm
% 6/26/98 made number of variables output 8, as above -sm 
% 9/16/98 plot out-of-bounds sortvars at nearest times boundary -sm
% 10/27/98 added cohsig, alpha -sm
% 10/28/98 adjust maxamp, maxcoh computation -sm
% 05/03/99 added horizontal ticks beneath coher trace, fixed vert. 
%          scale printing -t-pj
% 05/07/99 converted amps plot to log scaling -sm
% 05/14/99 added sort by alpha phase -se
% 07/23/99 made "Phase-sorted" axis label -sm
% 07/24/99 added 'allamps' option -sm
% 08/04/99 added new times spec., 'srate' arg, made 'phase' and 'allamps'
%          work together, plot re-aligned time zeros  -sm
% 06/26/99 debugged cbar; added vert lines at aligntime to plot1erp() axes -sm
% 09/29/99 fixed srate computation from times -sm & se
% 01/18/00 output outsort without clipping -sm
% 02/29/00 added 'vert' arg, fixed xticklabels, added ampsig -sm
% 03/03/00 added 'vert' arg lines to erp/amp/coher axes -sm
% 03/17/00 added axcopy -sm & tpj
% 03/29/00 added 'topo' option -sm 
% 05/05/00 fixed y-axis label bug when time limits given -sm
% 06/01/00 added topphase arg to 'phase' option for phasemovie.m -sm
% 07/12/00 adjusted prctle()
% 07/12/00 added 'spec' option -sm
% 08/22/00 added coherfreq to limits output -sm
% 09/13/00 added hard limit (1) to maxcoh -sm
% 09/14/00 made topoplot() and psd() plots relative to gca, not gcf -sm
% 10/10/00 added NoTimeflag -sm
% 11/03/00 changed method of rejecting small amplitude trials for phase sort -sm
% 11/03/00 added number_of_trials_out option for decfactor -sm
% 11/16/00 added ampoffset to center sig lines around baseline mean amp (0) -sm
% 01/06/01 edited help message; adjusted ampsig plot limits;initialized outputs -sm
%          rm'd 'allcoher' from help message - not fully implemented -sm
% 01/09/01 documented limits arg 'bamp' (baseline amplitude) -sm
% 02/13/01 debugged use of stored baseamp, ampsig parameters -sm
% 03/28/01 made erpimage(data) possible. Debugged ampsig change limits -sm
% 08/31/01 fixed allamps bug -sm
% 09/03/01 added 'auxvar' plotting -sm
% 01-25-02 reformated help & license, added links -ad
% 02-16-02 added matrix option to arg 'vert' -sm
% 04-05-02 corrected zero alignment problem (display only) -ad
%
% Known Bugs:
% 'limits', [lotime hitime] may not work with 'erp'
% 'limits', [... loerp hierp] (still??) may leave "ghost" grey numbers 
%       on the coher axis when printed (-djpeg or -depsc)
% 'allcohers' - not fully implemented, and has been omitted from the help msg

function [data,outsort,outtrials,limits,axhndls,erp,amps,cohers,cohsig,ampsig,allamps,phaseangles,phsamp,sortidx,erpsig] = erpimage(data,sortvar,times,titl,avewidth,decfactor,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24,arg25,arg26)

%
%%%%%%%%%%%%%%%%%%% Define defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize optional output variables:
erp = []; amps = []; cohers = []; cohsig = []; ampsig = []; 
allamps = []; phaseangles = []; phsamp = []; sortidx = [];
auxvar = []; erpsig = []; winloc = [];winlocs = [];
curfig = gcf;   % note current figure - to avoid v7.0.0 bug that draws
                % some elements on the EEGLAB window -sm 8-30-04

YES = 1;  % logical variables
NO  = 0;

DEFAULT_BASELINE_END = 0; % ms
TIMEX = 1;          % 1 -> plot time on x-axis; 
                    % 0 -> plot trials on x-axis

BACKCOLOR = [0.8 0.8 0.8]; % grey background
try, icadefs; catch, end;
                    % read BACKCOLOR for plot from defs file (edit this)
                    % read DEFAULT_SRATE for coher,phase,allamps, etc.

% Fix plotting text and line style parameters
SORTWIDTH = 2.5;    % Linewidth of plotted sortvar
ZEROWIDTH = 3.0;    % Linewidth of vertical 0 line
VERTWIDTH = 2.5;    % Linewidth of optional vertical lines
HORZWIDTH = 2.1;    % Linewidth of optional vertical lines
SIGNIFWIDTH = 1.9;  % Linewidth of red significance lines for amp, coher
DOTSTYLE   = 'k--'; % line style to use for vetical dotted/dashed lines
LINESTYLE = '-';    % solid line
LABELFONT = 14;     % font sizes for axis labels, tick labels
TICKFONT  = 11;

PLOT_HEIGHT = 0.2;  % fraction of y dim taken up by each time series axes
YGAP = 0.03;        % fraction gap between time axes
YEXPAND = 1.3;      % expansion factor for y-axis about erp, amp data limits

DEFAULT_AVEWIDTH  = 1; % smooth trials with this window size by default
DEFAULT_DECFACTOR = 1; % decimate by this factor by default
DEFAULT_CYCLES    = 3; % use this many cycles in amp,coher computation window
DEFAULT_CBAR      = NO;% do not plot color bar by default
DEFAULT_PHARGS = [0 25 8 13]; % Default arguments for phase sorting
DEFAULT_ALPHA     = 0.01;
alpha     = 0;      % default alpha level for coherence significance

MIN_ERPALPHA = 0.001; % significance bounds for ERP 
MAX_ERPALPHA = 0.1; 

Noshow    = NO;     % show sortvar by default
Nosort    = NO;     % sort on sortvar by default
Caxflag   = NO;     % use default caxis by default
Caxis     = [];
caxfraction = [];
Coherflag = NO;     % don't compute or show amp,coher by default
Cohsigflag= NO;     % default: do not compute coherence significance
Allampsflag=NO;     % don't image the amplitudes by default
Allcohersflag=NO;   % don't image the coherence amplitudes by default
Topoflag  = NO;     % don't plot a topoplot in upper left
Specflag  = NO;     % don't plot a spectrum in upper right
Erpflag   = NO;     % don't show erp average by default
Erpstdflag= NO; 
Erpalphaflag= NO; 
Alignflag = NO;     % don't align data to sortvar by default
Colorbar  = NO;     % if YES, plot a colorbar to right of erp image
Limitflag = NO;     % plot whole times range by default
Phaseflag = NO;     % don't sort by phase
Ampflag   = NO;     % don't sort by amplitude
Sortwinflag = NO;   % sort by amplitude over a window
Valflag   = NO;     % don't sort by value
Srateflag = NO;     % srate not given
Vertflag  = NO;
Horzflag  = NO;
Noshowflag  = NO;
Renormflag = NO;
Showwin = NO;
% yerplabel = '\muV';
yerplabel = 'ERP';
yerplabelflag = NO;
verttimes = [];
horzepochs = [];
NoTimeflag= NO;     % by default DO print "Time (ms)" below bottom axis
Signifflag= NO;     % compute significance instead of receiving it
Auxvarflag= NO;
Cycleflag = NO;
signifs   = NaN;
coherfreq = nan;    % amp/coher-calculating frequency
freq = 0;           % phase-sorting frequency
srate = DEFAULT_SRATE; % from icadefs.m
aligntime = nan;
timelimits= nan;
topomap   = [];     % topo map vector
lospecHz  = [];     % spec lo frequency
topphase = 180;     % default top phase for 'phase' option
renorm    = 'no';
noshow    = 'no';
Rmerp     = 'no';

minerp = NaN; % default limits
maxerp = NaN;
minamp = NaN;
maxamp = NaN;
mincoh = NaN;
maxcoh = NaN;
baseamp =NaN;
allamps = []; % default return 

ax1    = NaN; % default axes handles
axcb   = NaN;
ax2    = NaN;
ax3    = NaN;
ax4    = NaN;

%
%%%%%%%%%%%%%%%%%%% Test, fill in commandline args %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargin < 1
  help erpimage
  return
end

if nargin < 3 | isempty(times)
  if size(data,1)==1 | size(data,2)==1
   fprintf('erpimage(): either input a times vector or make data size = (frames,trials).\n')
   return
  end
  times = 1:size(data,1);
  NoTimesPassed= 1;
end

if nargin < 2 | isempty(sortvar)
  sortvar = 1:size(data,2);
  Noshow = 1; % don't plot the dummy sortvar
end

framestot = size(data,1)*size(data,2);
ntrials = length(sortvar);
if ntrials < 2
  help erpimage
  fprintf('\nerpimage(): too few trials.\n');
  return
end

frames = floor(framestot/ntrials);
if frames*ntrials ~= framestot
  help erpimage
  fprintf(...
    '\nerpimage(); length of sortvar doesn''t divide number of data elements??\n')
  return
end

if nargin < 6
  decfactor = 0;
end
if nargin < 5
  avewidth = 0;
end
if nargin<4
  titl = ''; % default no title
end
if nargin<3
  times = NO;
end
if length(times) == 1 | times == NO,  % make default times
   times = 0:frames-1;
   srate = 1000*(length(times)-1)/(times(length(times))-times(1));
   fprintf('Using sampling rate %g Hz.\n',srate);
elseif length(times) == 3
   mintime = times(1);
   frames = times(2);
   srate = times(3);
   times = mintime:1000/srate:mintime+(frames-1)*1000/srate;
   fprintf('Using sampling rate %g Hz.\n',srate);
else
   % Note: might use default srate read from icadefs here...
   srate = 1000*(length(times)-1)/(times(end)-times(1));
end
if length(times) ~= frames
   fprintf(...
'erpimage(): length(data)(%d) ~= length(sortvar)(%d) * length(times)(%d).\n\n',...
                  framestot,              length(sortvar),   length(times));
   return
end
if avewidth == 0,
  avewidth = DEFAULT_AVEWIDTH;
end
if decfactor == 0,
  decfactor = DEFAULT_DECFACTOR;
end
if avewidth < 1
  help erpimage
  fprintf('\nerpimage(): Variable avewidth cannot be < 1.\n')
  return
end
if avewidth > ntrials
  fprintf('Setting variable avewidth to max %d.\n',ntrials)
  avewidth = ntrials;  
end
if decfactor > ntrials
  fprintf('Setting variable decfactor to max %d.\n',ntrials)
  decfactor = ntrials;  
end
%
%%%%%%%%%%%%%%%%% Collect optional args %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargin > 6
  flagargs = [];

  for a=7:nargin % for each remaining Arg

	  Arg = eval(['arg' int2str(a-6)]);
	  if Caxflag == YES
		  if size(Arg,1) ~= 1 | size(Arg,2) > 2
			  help erpimage
			  fprintf('\nerpimage(): caxis arg must be a scalar or (1,2) vector.\n');
			  return
		  end
		  if size(Arg,2) == 2
			  Caxis = Arg;
		  else
			  caxfraction = Arg;
		  end
		  Caxflag = NO;
	  elseif Coherflag == YES
		  if length(Arg) > 3 | length(Arg) < 1
			  help erpimage
			  fprintf('\nerpimage(): coher arg must be size <= 3.\n');
			  return
		  end
		  coherfreq = Arg(1);
		  if size(Arg,2) == 1
			  coherfreq = Arg(1);
		  else
			  coherfreq = Arg(1:2);
		  end;
		  if size(Arg,2) == 3
			  Cohsigflag = YES;
			  alpha  = Arg(3);
			  if alpha < 0 | alpha > 0.1
				  fprintf('erpimage(): alpha value %g out of bounds.\n',alpha); 
				  return
			  end
		  end
		  Coherflag = NO;
		  Erpflag = YES;  % plot amp, coher below erp time series
	  elseif Topoflag == YES;
		  if length(Arg) ~= 2
			  help erpimage
			  fprintf('\nerpimage(): topo arg must be a list of length 2.\n');
			  return
		  end
		  topomap = Arg{1};
		  eloc_file = Arg{2};
		  Topoflag = NO;
	  elseif Specflag == YES;
		  if length(Arg) ~= 2
			  help erpimage
			  fprintf('\nerpimage(): spec arg must be a list of length 2.\n');
			  return
		  end
		  lospecHz = Arg(1);
		  hispecHz = Arg(2);
		  Specflag = NO;
	  elseif Renormflag == YES
		  renorm = Arg;
		  Renormflag = NO;
	  elseif Noshowflag == YES
		  noshow = Arg;
		  Noshowflag = NO;
	  elseif Alignflag == YES
		  aligntime = Arg;
		  Alignflag = NO;
	  elseif Limitflag == YES
		  %  [lotime hitime loerp hierp loamp hiamp locoher hicoher]
		  if size(Arg,1) ~= 1 | size(Arg,2) < 2 ...
				  | size(Arg,2) > 9 ...
				  help erpimage
			  fprintf('\nerpimage(): limits arg must be a vector sized (1,2<->9).\n');
			  return
		  end
		  if  ~isnan(Arg(1)) & (Arg(2) <= Arg(1))
			  help erpimage
			  fprintf('\nerpimage(): time limits out of order or out of range.\n');
			  return
		  end
		  if Arg(1) < min(times)
			  Arg(1) = min(times);
			  fprintf('Adjusting mintime limit to first data value %g\n',min(times));
		  end
		  if Arg(2) > max(times)
			  Arg(2) = max(times);
			  fprintf('Adjusting maxtime limit to last data value %g\n',max(times));
		  end
		  timelimits = Arg(1:2);
		  if length(Arg)> 2
			  minerp = Arg(3);
		  end
		  if length(Arg)> 3
			  maxerp = Arg(4);
		  end
		  if ~isnan(maxerp) & maxerp <= minerp
			  help erpimage
			  fprintf('\nerpimage(): erp limits args out of order.\n');
			  return
		  end
		  if length(Arg)> 4
			  minamp = Arg(5);
		  end
		  if length(Arg)> 5
			  maxamp = Arg(6);
		  end
		  if maxamp <= minamp
			  help erpimage
			  fprintf('\nerpimage(): amp limits args out of order.\n');
			  return
		  end
		  if length(Arg)> 6
			  mincoh = Arg(7);
		  end
		  if length(Arg)> 7
			  maxcoh = Arg(8);
		  end
		  if maxcoh <= mincoh
			  help erpimage
			  fprintf('\nerpimage(): coh limits args out of order.\n');
			  return
		  end
		  if length(Arg)>8
			  baseamp = Arg(9);    % for 'allamps'
		  end
		  Limitflag = NO;
		  
	  elseif Srateflag == YES
          srate = Arg(1);
          Srateflag = NO;
	  elseif Cycleflag == YES
          DEFAULT_CYCLES = Arg;
          Cycleflag = NO;
	  elseif Auxvarflag == YES;
          if isa(Arg,'cell')==YES & length(Arg)==2
			  auxvar = Arg{1};
			  auxcolors = Arg{2};
          elseif isa(Arg,'cell')==YES
			  fprintf('erpimage(): auxvars argument must be a matrix or length-2 cell array.\n');
			  return
          else
			  auxvar = Arg; % no auxcolors specified
          end
          Auxvarflag = NO;
	  elseif Vertflag == YES
            verttimes = Arg;
            Vertflag = NO;
	  elseif Horzflag == YES
            horzepochs = Arg;
            Horzflag = NO;
	  elseif yerplabelflag == YES
            yerplabel = Arg;
            yerplabelflag = NO;
	  elseif Signifflag == YES
		  signifs = Arg; % [low_amp hi_amp coher]
		  if length(signifs) ~= 3
			  fprintf('\nerpimage(): signif arg [%g] must have 3 values\n',Arg);
			  return
		  end
		  Signifflag = NO;
	  elseif Allcohersflag == YES
            data2=Arg;
            if size(data2) ~= size(data)
			  fprintf('erpimage(): allcohers data matrix must be the same size as data.\n');
              return
            end
            Allcohersflag = NO;
	  elseif Phaseflag == YES
            n = length(Arg);
            if n > 5
			  error('erpimage(): Too many arguments for keyword ''phasesort''');
            end
            phargs = Arg;
		  
            if phargs(3) < 0
              error('erpimage(): Invalid negative frequency argument for keyword ''phasesort''');
            end
            if n>=4
			  if phargs(4) < 0
				  error('erpimage(): Invalid negative argument for keyword ''phasesort''');
			  end
            end
            if min(phargs(1)) < times(1) | max(phargs(1)) > times(end)
			  error('erpimage(): time for phase sorting filter out of bound.');
            end
            if phargs(2) >= 100 | phargs(2) < -100
			  error('%-argument for keyword ''phasesort'' must be (-100;100)');
            end
            if length(phargs) >= 4 & phargs(3) > phargs(4)
			  error('erpimage(): Phase sorting frequency range must be increasing.');
            end
            if length(phargs) == 5
			  topphase = phargs(5);
            end
            Phaseflag = NO;
	  elseif Sortwinflag == YES % 'ampsort' mean amplitude over a time window
          n = length(Arg);
          sortwinarg = Arg;
          if n > 2
			  error('erpimage(): Too many arguments for keyword ''sortwin''');
          end
          if min(sortwinarg(1)) < times(1) | max(sortwinarg(1)) > times(end)
			  		error('erpimage(): start time for value sorting out of bounds.');
          end
          if n > 1
           if min(sortwinarg(2)) < times(1) | max(sortwinarg(2)) > times(end)
			  		error('erpimage(): end time for value sorting out of bounds.');
           end
          end
          if n > 1 & sortwinarg(1) > sortwinarg(2)
			  		error('erpimage(): Value sorting time range must be increasing.');
          end
          Sortwinflag = NO;
	  elseif Ampflag == YES % 'ampsort',[center_time,prcnt_reject,minfreq,maxfreq]
          n = length(Arg);
          if n > 4
			  error('erpimage(): Too many arguments for keyword ''ampsort''');
          end
          ampargs = Arg;
		  
          % if ampargs(3) < 0
           %    error('erpimage(): Invalid negative argument for keyword ''ampsort''');
          % end
          if n>=4
			  		if ampargs(4) < 0
				  		error('erpimage(): Invalid negative argument for keyword ''ampsort''');
			  		end
          end
          
          if ~isinf(ampargs(1))
             if min(ampargs(1)) < times(1) | max(ampargs(1)) > times(end)
			  		error('erpimage(): time for amplitude sorting filter out of bounds.');
             end
          end

          if ampargs(2) >= 100 | ampargs(2) < -100
			  		error('percentile argument for keyword ''ampsort'' must be (-100;100)');
          end
          
          if length(ampargs) == 4 & abs(ampargs(3)) > abs(ampargs(4))
			  		error('erpimage(): Amplitude sorting frequency range must be increasing.');
          end
          Ampflag = NO;

	  elseif Valflag == YES % sort by potential value in a given window
          % Usage: 'valsort',[mintime,maxtime,direction]
          n = length(Arg);
          if n > 3
			  error('erpimage(): Too many arguments for keyword ''valsort''');
          end
          valargs = Arg;
		  
          if min(valargs(1)) < times(1) | max(valargs(1)) > times(end)
			  		error('erpimage(): start time for value sorting out of bounds.');
          end
          if n > 1
           if min(valargs(2)) < times(1) | max(valargs(2)) > times(end)
			  		error('erpimage(): end time for value sorting out of bounds.');
           end
          end
          if n > 1 & valargs(1) > valargs(2)
			  		error('erpimage(): Value sorting time range must be increasing.');
          end
          if n==3 & (~isnumeric(valargs(3)) | valargs(3)==0)
			  		error('erpimage(): Value sorting direction must be +1 or -1.');
          end
          Valflag = NO;
      elseif Erpalphaflag == YES
          erpalpha = Arg(1);
          if erpalpha < MIN_ERPALPHA | erpalpha > MAX_ERPALPHA
             fprintf('erpimage(): erpalpha value is out of bounds [%g, %g]\n',...
                               MIN_ERPALPHA,MAX_ERPALPHA);
             return
          end
          Erpalphaflag = NO;
	  elseif strcmp(Arg,'nosort')
		  Nosort = YES;
	  elseif strcmp(Arg,'showwin')
		  Showwin = YES;
	  elseif strcmp(Arg,'renorm')
		  Renormflag = YES;
	  elseif strcmp(Arg,'noshow')
		  Noshowflag = YES;
	  elseif strcmp(Arg,'noplot')|strcmp(Arg,'noshow')
		  Noshow = YES;
	  elseif strcmp(Arg,'caxis')
		  Caxflag = YES;
	  elseif strcmp(Arg,'coher')
		  Coherflag = YES;
	  elseif (strcmp(Arg,'allamps') | strcmp(Arg,'plotamps'))
		  Allampsflag = YES;
	  elseif strcmp(Arg,'allcohers')
		  Allcohersflag = YES;
	  elseif strcmp(Arg,'topo') | strcmp(Arg,'topoplot')
		  Topoflag = YES;
	  elseif strcmp(Arg,'spec') | strcmp(Arg,'spectrum')
		  Specflag = YES;
	  elseif strcmp(Arg,'erp')| strcmp(Arg,'ERP')
		  Erpflag = YES;
	  elseif strcmpi(Arg,'erpstd')
		  Erpstdflag = YES;
	  elseif strcmpi(Arg,'erpalpha')
		  Erpalphaflag = YES;
	  elseif strcmpi(Arg,'rmerp')
		  Rmerp = 'yes';
	  elseif strcmp(Arg,'align')
		  Alignflag = YES;
	  elseif strcmp(Arg,'cbar') | strcmp(Arg,'colorbar')
		  Colorbar = YES;
	  elseif strcmp(Arg,'limits')
		  Limitflag = YES;
	  elseif (strcmp(Arg,'phase') | strcmp(Arg,'phasesort'))
		  Phaseflag = YES;
	  elseif strcmp(Arg,'ampsort') 
		  Ampflag = YES;
	  elseif strcmp(Arg,'sortwin') 
		  Sortwinflag = YES;
	  elseif strcmp(Arg,'valsort')
		  Valflag = YES;
	  elseif strcmp(Arg,'auxvar')
		  Auxvarflag = YES;
	  elseif strcmp(Arg,'cycles')
		  Cycleflag = YES;
	  elseif strcmp(Arg,'yerplabel')
		  yerplabelflag = YES;
	  elseif strcmp(Arg,'srate')
		  Srateflag = YES;
	  elseif strcmp(Arg,'vert') |  strcmp(Arg,'verttimes')
		  Vertflag = YES;
	  elseif strcmp(Arg,'horz') |  strcmp(Arg,'horiz') | strcmp(Arg,'horizontal')
		  Horzflag = YES;
	  elseif strcmp(Arg,'signif')|strcmp(Arg,'signifs')|strcmp(Arg,'sig')|strcmp(Arg,'sigs')
		  Signifflag = YES;
	  elseif strcmp(Arg,'noxlabel') | strcmp(Arg,'noxlabels') | strcmp(Arg,'nox')
		  NoTimeflag = YES;
	  else
		  help erpimage
		  if isstr(Arg)
			  fprintf('\nerpimage(): unknown arg %s\n',Arg);
		  else
			  fprintf('\nerpimage(): unknown arg %d, size(%d,%d)\n',a,size(Arg,1),size(Arg,2));
		  end
		  return
	  end
  end % Arg
end

if   Caxflag == YES ...
  |Coherflag == YES ...
  |Alignflag == YES ...
  |Limitflag == YES
    help erpimage
    fprintf('\nerpimage(): missing option arg.\n')
    return
end
if (Allampsflag | exist('data2')) & ( any(isnan(coherfreq)) | ~Cohsigflag )
	fprintf('\nerpimage(): allamps and allcohers flags require coher freq, srate, and cohsig.\n');
	return
end
if Allampsflag & exist('data2')
	fprintf('\nerpimage(): cannot image both allamps and allcohers.\n');
	return
end
if ~exist('srate') | srate <= 0 
	fprintf('\nerpimage(): Data srate must be specified and > 0.\n');
	return
end
if ~isempty(auxvar)
	% whos auxvar
    if size(auxvar,1) ~= ntrials & size(auxvar,2) ~= ntrials
		fprintf('erpimage(): auxvar size should be (N,ntrials), e.g., (N,%d)\n',...
                           ntrials);
        return
    end
	if size(auxvar,1) == ntrials & size(auxvar,2) ~= ntrials  % make (N,frames)
		auxvar = auxvar';               
	end
	if size(auxvar,2) ~= ntrials
		fprintf('erpimage(): auxvar size should be (N,ntrials), e.g., (N,%d)\n',...
                           ntrials);
		return
	end
	if exist('auxcolors')==YES % if specified
		if isa(auxcolors,'cell')==NO % if auxcolors is not a cell array
			fprintf(...
                 'erpimage(): auxcolors argument to auxvar flag must be a cell array.\n');
			return
		end
	end
end
if exist('phargs')
	if phargs(3) > srate/2
	   fprintf(...
            'erpimage(): Phase-sorting frequency (%g Hz) must be less than Nyquist rate (%g Hz).',...
                phargs(3),srate/2);
	end
                                    % DEFAULT_CYCLES = 9*phargs(3)/(phargs(3)+10); % 3 cycles at 5 Hz
	if frames < DEFAULT_CYCLES*srate/phargs(3)
		fprintf('\nerpimage(): phase-sorting freq. (%g) too low: epoch length < %d cycles.\n',...
				phargs(3),DEFAULT_CYCLES);
		return
	end
	if length(phargs)==4 & phargs(4) > srate/2
		phargs(4) = srate/2;
	end
	if length(phargs)==5 & (phargs(5)>180 | phargs(5) < -180)
		fprintf('\nerpimage(): coher topphase (%g) out of range.\n',topphase);
		return
	end
end
if exist('ampargs')
	if abs(ampargs(3)) > srate/2
    	  fprintf(...
           'erpimage(): amplitude-sorting frequency (%g Hz) must be less than Nyquist rate (%g Hz).',...
              abs(ampargs(3)),srate/2);
	end
    % DEFAULT_CYCLES = 9*abs(ampargs(3))/(abs(ampargs(3))+10); % 3 cycles at 5 Hz
	if frames < DEFAULT_CYCLES*srate/abs(ampargs(3))
		fprintf('\nerpimage(): amplitude-sorting freq. (%g) too low: epoch length < %d cycles.\n',...
				abs(ampargs(3)),DEFAULT_CYCLES);
		return
	end
	if length(ampargs)==4 & abs(ampargs(4)) > srate/2
		ampargs(4) = srate/2;
  fprintf('> Reducing max ''ampsort'' frequency to Nyquist rate (%g Hz)\n',srate/2)
	end
end
if ~any(isnan(coherfreq))
	if coherfreq(1) <= 0 | srate <= 0
		fprintf('\nerpimage(): coher frequency (%g) out of range.\n',coherfreq(1));
		return
	end
	if coherfreq(end) > srate/2 | srate <= 0
		fprintf('\nerpimage(): coher frequency (%g) out of range.\n',coherfreq(end));
		return
	end
    %DEFAULT_CYCLES = 9*coherfreq(1)/(coherfreq(1)+10); % 3 cycles at 5 Hz
	if frames < DEFAULT_CYCLES*srate/coherfreq(1)
		fprintf('\nerpimage(): coher freq. (%g) too low:  epoch length < %d cycles.\n',...
				coherfreq(1),DEFAULT_CYCLES);
		return
	end
end
          
if isnan(timelimits)
   timelimits = [min(times) max(times)];
end
if ~isstr(aligntime) & ~isnan(aligntime)
	if ~isinf(aligntime) ...
			& (aligntime < timelimits(1) | aligntime > timelimits(2))
		help erpimage
		fprintf('\nerpimage(): requested align time outside of time limits.\n');
		return
	end
end
% 
%%%%%%%%%%%%%%%%  Replace nan's with 0s %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
nans= find(isnan(data));
if length(nans)
  fprintf('Replaced %d nan in data with 0s.\n');
  data(nans) = 0;
end
%
%%%%%%%%%%%%%% Reshape data to (frames,ntrials) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if size(data,2) ~= ntrials
	if size(data,1)>1
		% fprintf('frames %d, ntrials %d length(data) %d\n',frames,ntrials,length(data));
		data=reshape(data,1,frames*ntrials);
   end
   data=reshape(data,frames,ntrials);
end
fprintf('\nPlotting input data as %d epochs of %d frames sampled at %3.1f Hz.\n',...
                             ntrials,frames,srate);
%
%%%%%%%%%%%%%% Reshape data2 to (frames,ntrials) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if exist('data2') == 1
  if size(data2,2) ~= ntrials
   if size(data2,1)>1
     data2=reshape(data2,1,frames*ntrials);
   end
   data2=reshape(data2,frames,ntrials);
  end
end
%
%%%%%%%%%%%%%%% if sortvar=NaN, remove lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -ad
%
if any(isnan(sortvar))
	nanlocs = find(isnan(sortvar));
	fprintf('Removing %d trials with NaN sortvar values.\n', length(nanlocs));
	data(:,nanlocs) = [];
	sortvar(nanlocs) = [];
    if length(sortvar) < 4 & avewidth > 1
        error('Not enough trials');
    end;
	if exist('data2') == 1
		data2(:,nanlocs) = [];
	end;
	if ~isempty(auxvar)
		auxvar(:,nanlocs) = [];
	end
	if ~isempty(verttimes)
		if size(verttimes,1) == ntrials
			verttimes(nanlocs,:) = [];
		end;
	end;
	ntrials = size(data,2);
	if ntrials <= 1, close(gcf); error('Too few trials'); end;
end;

%
%%%%%%%%%%%%%%%%%%% Renormalize sortvar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
switch lower(renorm)
 case 'yes',
  disp('erpimage warning: *** sorting variable renormalized ***');  
  sortvar = (sortvar-min(sortvar)) / (max(sortvar) - min(sortvar)) * ...
		   0.5 * (max(times) - min(times)) + min(times) + 0.4*(max(times) - min(times));
 case 'no',;
 otherwise,
  if ~isempty(renorm)
      locx = findstr('x', lower(renorm));
      if length(locx) ~= 1, error('erpimage: unrecognized renormalizing formula'); end;
      eval( [ 'sortvar =' renorm(1:locx-1) 'sortvar' renorm(locx+1:end) ';'] );
  end;
end;
%
%%%%%%%%%%%%%%%%%%% Align data to sortvar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if isstr(aligntime) | ~isnan(aligntime)
    if ~isstr(aligntime) & isinf(aligntime)
        aligntime= median(sortvar);
        fprintf('Aligning data to median sortvar.\n'); 
        % Alternative below: trimmed median - ignore top/bottom 5%
        %   ssv = sort(sortvar); % ssv = 'sorted sortvar'
        %   aligntime= median(ssv(ceil(ntrials/20)):floor(19*ntrials/20)); 
    end
    
    if ~isstr(aligntime)
        fprintf('Realigned sortvar plotted at %g ms.\n',aligntime);
        aligndata=zeros(frames,ntrials); % begin with matrix of zeros()
        shifts = zeros(1,ntrials);
        for t=1:ntrials, %%%%%%%%% foreach trial %%%%%%%%%
            shft = round((aligntime-sortvar(t))*srate/1000);
            shifts(t) = shft;
            if shft>0, % shift right
                if frames-shft > 0
                    aligndata(shft+1:frames,t)=data(1:frames-shft,t);
                else
                    fprintf('No aligned data for epoch %d - shift (%d frames) too large.\n',t,shft);
                end
            elseif shft < 0 % shift left
                if frames+shft > 0
                    aligndata(1:frames+shft,t)=data(1-shft:frames,t);
                else
                    fprintf('No aligned data for epoch %d - shift (%d frames) too large.\n',t,shft);
                end
            else % shft == 0
                aligndata(:,t) = data(:,t);
            end 
        end % end trial
        fprintf('Shifted epochs by %d to %d frames.\n',min(shifts),max(shifts));
        data = aligndata;                       % now data is aligned to sortvar
    else
        aligntime = str2num(aligntime);
        if isinf(aligntime),  aligntime= median(sortvar); end;
    end;
end 

%
%%%%%%%%%%%%%%%%%%%%%% Remove the ERP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmpi(Rmerp, 'yes')
    data = data - repmat(nan_mean(data')', [1 size(data,2)]);
end;

%
%%%%%%%%%%%%%%% Sort the data trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if exist('phargs') == 1 % if phase-sort
	if length(phargs) >= 4 && phargs(3) ~= phargs(4) % find max frequency in specified band
        if exist('psd') == 2 % requires Singla Processing Toolbox
            [pxx,freqs] = psd(data(:),max(1024, pow2(ceil(log2(frames)))),srate,frames,0);
        else
            [pxx,freqs] = spec(data(:),max(1024, pow2(ceil(log2(frames)))),srate,frames,0);
        end;
	%gf = gcf; % figure;plot(freqs,pxx); %xx=axis; %axis([phargs(3) phargs(4) xx(3) xx(4)]); %figure(gf);
		
		pxx = 10*log10(pxx);
		n = find(freqs >= phargs(3) & freqs <= phargs(4));
		if ~length(n)
			freq = phargs(3);
		end
		[dummy maxx] = max(pxx(n));
		freq = freqs(n(maxx));
	else
		freq = phargs(3); % else use specified frequency
	end
        
       phwin = phargs(1);
       [dummy minx] = min(abs(times-phwin)); % closest time to requested
       winlen = floor(DEFAULT_CYCLES*srate/freq);
       winloc = minx-linspace(floor(winlen/2), floor(-winlen/2), winlen+1); 
       tmprange = find(winloc>0 & winloc<=frames);
       winloc = winloc(tmprange); % sorting window times
       [phaseangles phsamp] = phasedet(data,frames,srate,winloc,freq);
       winlocs = winloc;
	
    if length(tmprange) ~=  winlen+1
        filtersize = DEFAULT_CYCLES * length(tmprange) / (winlen+1);
        timecenter = median(winloc)/srate*1000+times(1); % center of window in ms
        phaseangles = phaseangles + 2*pi*(timecenter-phargs(1))*freq;
        fprintf('Sorting data epochs by phase at frequency %2.1f Hz: \n', freq);
        fprintf('    Data time limits reached -> now uses a %1.1f cycles (%1.0f ms) window centered at %1.0f ms\n', ...
                filtersize, 1000/freq*filtersize, timecenter);
        fprintf(...
  '    Filter length is %d; Phase has been linearly interpolated to latency at %1.0f ms.\n', ...
                        length(winloc), phargs(1));
    else
        fprintf(...
  'Sorting data epochs by phase at %2.1f Hz in a %1.1f-cycle (%1.0f ms) window centered at %1.0f ms.\n',...  
			freq,DEFAULT_CYCLES,1000/freq*DEFAULT_CYCLES,times(minx));
        fprintf('Phase is computed using a wavelet of %d frames.\n',length(winloc));
    end;
	%
	% Reject small (or large) phsamp trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	phargs(2) = phargs(2)/100; % convert rejection rate from % to fraction
	[tmp n] = sort(phsamp); % sort amplitudes
	if phargs(2)>=0
	  n = n(ceil(phargs(2)*length(n))+1:end); % if rej 0, select all trials
	  fprintf('Retaining %d epochs (%g percent) with largest power at the analysis frequency,\n',...
	     length(n),100*(1-phargs(2)));
	else % phargs(2) < 0
	   phargs(2) = 1+phargs(2); % subtract from end
	   n = n(1:floor(phargs(2)*length(n)));
	   fprintf('Retaining %d epochs (%g percent) with smallest power at the analysis frequency,\n',...
                      length(n),phargs(2)*100);
	end
	%
	% Remove low|high-amplitude trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	data = data(:,n); % amp-sort the data, removing rejected-amp trials
	phsamp = phsamp(n);           % amp-sort the amps
	phaseangles = phaseangles(n); % amp-sort the phaseangles
	sortvar = sortvar(n);         % amp-sort the trial indices
	ntrials = length(n);          % number of trials retained
	if ~isempty(auxvar)
	   auxvar = auxvar(:,n);
        end
	%
	% Sort remaining data by phase angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	topphase = topphase/360*2*pi; % convert from degrees to radians
	phaseangles = -phaseangles;
	ip = find(phaseangles>topphase);
	phaseangles(ip) = phaseangles(ip)-2*pi; % rotate so topphase at top of plot
	
	[phaseangles sortidx] = sort(phaseangles); % sort trials on (rotated) phase
	data    =  data(:,sortidx);                % sort data by phase
	phsamp  =  phsamp(sortidx);                % sort amps by phase
	sortvar = sortvar(sortidx);                % sort input sortvar by phase
	phaseangles = -phaseangles; % Note: phsangles now descend from pi 
	if ~isempty(auxvar)
	   auxvar = auxvar(:,sortidx);
	end
	
	fprintf('Size of data = [%d,%d]\n',size(data,1),size(data,2));
	sortidx = n(sortidx); % return original trial indices in final sorted order
%
% %%%%%%%%%%%%%%% Sort data by amplitude %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
elseif exist('ampargs') == 1 % if amplitude-sort
	if length(ampargs) == 4 % find max frequency in specified band
          if exist('psd') == 2
            [pxx,freqs] = psd(data(:),max(1024, pow2(ceil(log2(frames)))),srate,frames,0);
          else
            [pxx,freqs] = spec(data(:),max(1024, pow2(ceil(log2(frames)))),srate,frames,0);
          end;
	  pxx = 10*log10(pxx);
	  n = find(freqs >= abs(ampargs(3)) & freqs <= abs(ampargs(4)));
	  if ~length(n)
		  freq = abs(ampargs(3));
	  end
          if ampargs(3)>=0
	     [dummy maxx] = max(pxx(n));
	     freq = freqs(n(maxx)); % use the highest-power frequency
          else
             freq = freqs(n);  % use all frequencies in the specified range
          end
	else
	  freq = abs(ampargs(3)); % else use specified frequency
	end
        if length(freq) == 1
           fprintf('Sorting data epochs by amplitude at frequency %2.1f Hz \n', freq);
        else
           fprintf('Sorting data epochs by amplitude at %d frequencies (%2.1f Hz to %.1f Hz) \n',...
                                  length(freq),freq(1),freq(end));
        end
        SPECWININCR = 10;	% make spectral sorting time windows increment by 10 ms
        if isinf(ampargs(1))
            ampwins = sortwinarg(1):SPECWININCR:sortwinarg(2);
        else
            ampwins = ampargs(1);
        end
        if ~isinf(ampargs(1)) % single time given
            if length(freq) == 1
        	fprintf(...
'   in a %1.1f-cycle (%1.0f ms) time window centered at %1.0f ms.\n',...  
			DEFAULT_CYCLES,1000/freq(1)*DEFAULT_CYCLES,times(minx));
           else
        	fprintf(...
'   in %1.1f-cycle (%1.0f-%1.0f ms) time windows centered at %1.0f ms.\n',...  
	  DEFAULT_CYCLES,1000/freq(1)*DEFAULT_CYCLES,1000/freq(end)*DEFAULT_CYCLES,times(minx));
           end
        else % range of times
            [dummy sortwin_st ] = min(abs(times-ampwins(1)));
            [dummy sortwin_end] = min(abs(times-ampwins(end)));
            if length(freq) == 1
        	fprintf(...
'   in %d %1.1f-cycle (%1.0f ms) time windows centered from %1.0f to  %1.0f ms.\n',...  
			length(ampwins),DEFAULT_CYCLES,1000/freq(1)*DEFAULT_CYCLES,times(sortwin_st),times(sortwin_end));
           else
        	fprintf(...
'   in %d %1.1f-cycle (%1.0f-%1.0f ms) time windows centered from %1.0f to %1.0f ms.\n',...  
			length(ampwins),DEFAULT_CYCLES,1000/freq(1)*DEFAULT_CYCLES,1000/freq(end)*DEFAULT_CYCLES,times(sortwin_st),times(sortwin_end));
           end
        end
	
        phsamps = 0; %%%%%%%%%%%%%%%%%%%%%%%%%% sort by (mean) amplitude %%%%%%%%%%%%%%%%%%%%%%%%%%
        minxs = [];
        for f = 1:length(freq)  % use one or range of frequencies
         frq = freq(f);
         for ampwin = ampwins
	       [dummy minx] = min(abs(times-ampwin)); % find nearest time point to requested
           minxs = [minxs minx];
	       winlen = floor(DEFAULT_CYCLES*srate/frq);
	       % winloc = minx-[winlen:-1:0]; % ending time version
	       winloc = minx-linspace(floor(winlen/2), floor(-winlen/2), winlen+1);
           tmprange = find(winloc>0 & winloc<=frames);
           winloc = winloc(tmprange); % sorting window frames
           if f==1
              winlocs = [winlocs;winloc];  % store tme windows
           end
	       [phaseangles phsamp] = phasedet(data,frames,srate,winloc,frq);
           phsamps = phsamps+phsamp;  % accumulate amplitudes across 'sortwin'
	     end
        end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     	if length(tmprange) ~=  winlen+1 % ????????? 
       	 	filtersize = DEFAULT_CYCLES * length(tmprange) / (winlen+1);
        	timecenter = median(winloc)/srate*1000+times(1); % center of window in ms
        	phaseangles = phaseangles + 2*pi*(timecenter-ampargs(1))*freq(end);
        	fprintf(...
'    Data time limits reached -> now uses a %1.1f cycles (%1.0f ms) window centered at %1.0f ms\n', ...
                	filtersize, 1000/freq(1)*filtersize, timecenter);
        	fprintf(...
'    Wavelet length is %d; Phase has been linearly interpolated to latency et %1.0f ms.\n', ...
           		length(winloc(1,:)), ampargs(1));
        end 
        if length(freq) == 1
           fprintf('Amplitudes are computed using a wavelet of %d frames.\n',length(winloc(1,:)));
        end

	%
	% Reject small (or large) phsamp trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	ampargs(2) = ampargs(2)/100; % convert rejection rate from % to fraction
	[tmp n] = sort(phsamps); % sort amplitudes
	if ampargs(2)>=0
	   n = n(ceil(ampargs(2)*length(n))+1:end); % if rej 0, select all trials
	   fprintf('Retaining %d epochs (%g percent) with largest power at the analysis frequency,\n',...
	      length(n),100*(1-ampargs(2)));
	else % ampargs(2) < 0
		ampargs(2) = 1+ampargs(2); % subtract from end
		n = n(1:floor(ampargs(2)*length(n)));
		fprintf(...
            'Retaining %d epochs (%g percent) with smallest power at the analysis frequency,\n',...
                   length(n),ampargs(2)*100);
	end
	%
	% Remove low|high-amplitude trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	data = data(:,n); % amp-sort the data, removing rejected-amp trials
	phsamps = phsamps(n);           % amp-sort the amps
	phaseangles = phaseangles(n); % amp-sort the phaseangles
	sortvar = sortvar(n);         % amp-sort the trial indices
	ntrials = length(n);          % number of trials retained
	if ~isempty(auxvar)
           auxvar = auxvar(:,n);
        end
	%
	% Sort remaining data by amplitude %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	[phsamps sortidx] = sort(phsamps); % sort trials on amplitude
	data    =  data(:,sortidx);                % sort data by amp
	phaseangles  =  phaseangles(sortidx);      % sort angles by amp
	sortvar = sortvar(sortidx);                % sort input sortvar by amp
	if ~isempty(auxvar)
	   auxvar = auxvar(:,sortidx);
	end
	fprintf('Size of data = [%d,%d]\n',size(data,1),size(data,2));

	sortidx = n(sortidx); % return original trial indices in final sorted order
%
%%%%%%%%%%%%%%%%%%%%%% Don't Sort trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
elseif Nosort == YES
  fprintf('Not sorting data on input sortvar.\n');
  sortidx = 1:ntrials;	
%
%%%%%%%%%%%%%%%%%%%%%% Sort trials on (mean) value %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
elseif exist('valargs')
    [sttime stframe] = min(abs(times-valargs(1)));
    sttime = times(stframe);
    if length(valargs)>1
        [endtime endframe] = min(abs(times-valargs(2)));
        endtime = times(endframe);
    else
        endframe = stframe;
        endtime = times(endframe);
    end
    if length(valargs)==1 | sttime == endtime
        fprintf('Sorting data on value at time %4.0f ms.\n',sttime);
    elseif length(valargs)>1
        fprintf('Sorting data on mean value between %4.0f and %4.0f ms.\n',...
                sttime,endtime);
    end
    if endframe>stframe
        sortval = mean(data(stframe:endframe,:));
    else
        sortval = data(stframe,:);
    end
    [sortval,sortidx] = sort(sortval);
    if length(valargs)>2
        if valargs(3) <0
            sortidx = sortidx(end:-1:1); % plot largest values on top
                                         % if direction < 0   
        end
    end
    data = data(:,sortidx);
    sortvar = sortvar(sortidx);                % sort input sortvar by amp
    if ~isempty(auxvar)
        auxvar = auxvar(:,sortidx);
    end
    if ~isempty(phaseangles)
	phaseangles  =  phaseangles(sortidx);      % sort angles by amp
    end
    winloc = [stframe,endframe];
    winlocs = winloc;
%
%%%%%%%%%%%%%%%%%%%%%% Sort trials on sortvar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
else 
    fprintf('Sorting data on input sortvar.\n');
    [sortvar,sortidx] = sort(sortvar);
    data = data(:,sortidx);
    if ~isempty(auxvar)
        auxvar = auxvar(:,sortidx);
    end
end

%if max(sortvar)<0
%   fprintf('Changing the sign of sortvar: making it positive.\n');
%   sortvar = -sortvar;
%end
%
%%%%%%%%%%%%%%%%%%% Adjust decfactor if phargs or ampargs %%%%%%%%%%%%%%%%%%%%%
%
if decfactor < 0
    decfactor = -decfactor;
    invdec = 1;
else
    invdec = 0;
end;
if decfactor > sqrt(ntrials) % if large, output this many trials
    n = 1:ntrials;
    if exist('phargs') & length(phargs)>1
        if phargs(2)>0
            n = n(ceil(phargs(2)*ntrials)+1:end); % trials after rejection
        elseif phargs(2)<0
            n = n(1:floor(phargs(2)*length(n)));  % trials after rejection
        end
    elseif exist('ampargs') & length(ampargs)>1
        if ampargs(2)>0
            n = n(ceil(ampargs(2)*ntrials)+1:end); % trials after rejection
        elseif ampargs(2)<0
            n = n(1:floor(ampargs(2)*length(n)));  % trials after rejection
        end
    end
    if invdec
        decfactor = (length(n)-avewidth)/decfactor;
    else
        decfactor = length(n)/decfactor;
    end;
end

% 
%%%%%%%%%%%%%%%%%% Smooth data using moving average %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
urdata = data; % save data to compute amp, coher on unsmoothed data
if ~Allampsflag & ~exist('data2') % if imaging potential,
    if avewidth > 1 | decfactor > 1
        if Nosort == YES
            fprintf('Smoothing the data using a window width of %g epochs ',avewidth);
        else
            fprintf('Smoothing the sorted epochs with a %g-epoch moving window.',...
                    avewidth);
        end
        fprintf('\n');
        fprintf('  and a decimation factor of %g\n',decfactor);
        if ~exist('phargs') % if not phase-sorted trials
           [data,outtrials] = movav(data,1:ntrials,avewidth,decfactor); 
           % Note: movav() here sorts using square window
           [outsort,outtrials] = movav(sortvar,1:ntrials,avewidth,decfactor); 
        else % if phase-sorted trials, use circular / wrap-around smoothing
           backhalf  = floor(avewidth/2);
           fronthalf = floor((avewidth-1)/2);
           if avewidth > 2
            [data,outtrials] = movav([data(:,[(end-backhalf+1):end]),...
                                      data,...
                                      data(:,[1:fronthalf])],...
                                      [1:(ntrials+backhalf+fronthalf)],avewidth,decfactor); 
            [outsort,outtrials] = movav([sortvar((end-backhalf+1):end),...
                                        sortvar,...
                                        sortvar(1:fronthalf)],...
                                        1:(ntrials+backhalf+fronthalf),avewidth,decfactor); 
            % outtrials = 1:ntrials;
           else % avewidth==2
            [data,outtrials] = movav([data(:,end),data],...
                                       [1:(ntrials+1)],avewidth,decfactor); 
            % Note: movav here sorts using square window
            [outsort,outtrials] = movav([sortvar(end) sortvar],...
                                        1:(ntrials+1),avewidth,decfactor); 
            % outtrials = 1:ntrials;
           end
        end
        if ~isempty(auxvar)
          if ~exist('phargs') % if not phase-sorted trials
            [auxvar,tmp] = movav(auxvar,1:ntrials,avewidth,decfactor); 
          else % if phase-sorted trials
           if avewidth>2 
            [auxvar,tmp] = movav([auxvar(:,[(end-backhalf+1):end]),...
                                  auxvar,...
                                  auxvar(:,[1:fronthalf])],...
                                 [1:(ntrials+backhalf+fronthalf)],avewidth,decfactor); 
           else % avewidth==2
            [auxvar,tmp] = movav([auxvar(:,end),auxvar],[1:(ntrials+1)],avewidth,decfactor); 
           end
          end
        end
        fprintf('Output data will be %d frames by %d smoothed trials.\n',...
                frames,length(outtrials));
        fprintf('Outtrials: %3.2f to %4.2f\n',min(outtrials),max(outtrials));
    else % don't smooth
        outtrials = 1:ntrials;
        outsort = sortvar;
    end

    %
    %%%%%%%%%%%%%%%%%%%%%%%%% Find color axis limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if ~isempty(Caxis) 
        mindat = Caxis(1);
        maxdat = Caxis(2);
        fprintf('Using the specified caxis range of [%g,%g].\n', mindat, maxdat);
    else
        mindat = min(min(data));
        maxdat = max(max(data));
        maxdat =  max(abs([mindat maxdat])); % make symmetrical about 0
        mindat = -maxdat;
        if ~isempty(caxfraction)
            adjmax = (1-caxfraction)/2*(maxdat-mindat);
            mindat = mindat+adjmax;
            maxdat = maxdat-adjmax;
            fprintf(...
                'The caxis range will be %g times the sym. abs. data range -> [%g,%g].\n',...
                caxfraction,mindat,maxdat);
        else
            fprintf(...
                'The caxis range will be the sym. abs. data range -> [%g,%g].\n',...
                mindat,maxdat);
        end
    end
end % if ~Allampsflag & ~exist('data2')

%
%%%%%%%%%%%%%%%%%%%%%%%%%% Set time limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if isnan(timelimits(1))
    timelimits = [min(times) max(times)];
end
fprintf('Data will be plotted between %g and %g ms.\n',timelimits(1),timelimits(2));

%
%%%%%%%%%%%%% Image the aligned/sorted/smoothed data %%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmpi(noshow, 'no')  
    if ~any(isnan(coherfreq))       % if plot three time axes
        image_loy = 3*PLOT_HEIGHT;
    elseif Erpflag == YES   % elseif if plot only one time axes
        image_loy = 1*PLOT_HEIGHT;
    else                    % else plot erp-image only
        image_loy = 0*PLOT_HEIGHT;
    end
    gcapos=get(gca,'Position');
    delete(gca)
    if isempty(topomap)
        image_top = 1;
    else
        image_top = 0.9;
    end
    ax1=axes('Position',...
             [gcapos(1) gcapos(2)+image_loy*gcapos(4) ...
              gcapos(3) (image_top-image_loy)*gcapos(4)]);
end;    
ind = isnan(data);    % find nan's in data
[i j]=find(ind==1);
if ~isempty(i)
    data(i,j) = 0;      % plot shifted nan data as 0 (=green)
end

%
%%%%%%%%%%%%% Determine coherence freqeuncy %%%%%%%%%%%%%%%%%%%%%%%%%%
%
if length(coherfreq) == 2 & coherfreq(1) ~= coherfreq(2) & freq <= 0 
	% find max frequency in specified band
    if exist('psd') == 2
        [pxx,tmpfreq] = psd(data(:),max(1024,pow2(ceil(log2(frames)))),srate,frames,0);
    else
        [pxx,tmpfreq] = spec(data(:),max(1024,pow2(ceil(log2(frames)))),srate,frames,0);
    end;
    
	pxx = 10*log10(pxx);
	n = find(tmpfreq >= coherfreq(1) & tmpfreq <= coherfreq(2));
	if ~length(n)
		coherfreq = coherfreq(1);
	end
	[dummy maxx] = max(pxx(n));
	coherfreq = tmpfreq(n(maxx));	
else 
	coherfreq = coherfreq(1);
end

if ~Allampsflag & ~exist('data2') %%%%%%%% Plot ERP image %%%%%%%%%%

    if strcmpi(noshow, 'no')
        if TIMEX
            imagesc(times,outtrials,data',[mindat,maxdat]);% plot time on x-axis
            set(gca,'Ydir','normal');
            axis([timelimits(1) timelimits(2) ...
                  min(outtrials) max(outtrials)]);
        else
            imagesc(outtrials,times,data,[mindat,maxdat]); % plot trials on x-axis
            axis([min(outtrials) max(outtrials)...
                  timelimits(1) timelimits(2)]);
        end
        hold on
        drawnow
    end;
    
elseif Allampsflag %%%%%%%%%%%%%%%% Plot allamps instead of data %%%%%%%%%%%%%%

    if freq > 0 
        coherfreq = mean(freq); % use phase-sort frequency
    end
        
    if ~isnan(signifs) % plot received significance levels
        fprintf('Computing and plotting received amp and ITC signif. levels...\n');
        [amps,cohers,cohsig,ampsig,allamps] = ...
            phasecoher(urdata,length(times),srate,coherfreq,DEFAULT_CYCLES,0);
        % need to receive cohsig and ampsig to get allamps
        ampsig = signifs([1 2]); % assume these already in dB
        cohsig = signifs(3);
        
    elseif alpha>0 % compute significance levels
        fprintf('Computing and plotting %g amp and ITC signif. level...\n',alpha);
        [amps,cohers,cohsig,ampsig,allamps] = ...
            phasecoher(urdata,length(times),srate,coherfreq,DEFAULT_CYCLES,alpha);
        % need to receive cohsig and ampsig to get allamps
        ampsig = 20*log10(ampsig); % convert to dB
        fprintf('Coherence significance level: %g\n',cohsig);

    else % no plotting of significance
        [amps,cohers,cohsig,ampsig,allamps] = ...
            phasecoher(urdata,length(times),srate,coherfreq,DEFAULT_CYCLES,0);
        % need to receive cohsig and ampsig to get allamps
    end
    
    % fprintf('#1 Size of allamps = [%d %d]\n',size(allamps,1),size(allamps,2));
    
    base = find(times<=DEFAULT_BASELINE_END);
    if length(base)<2
        base = 1:floor(length(times)/4); % default first quarter-epoch
        fprintf('Using %g to %g ms as amplitude baseline.\n',...
                times(1),times(base(end)));
    end
    amps = 20*log10(amps); % convert to dB
    
    if alpha>0
        fprintf('Amplitude significance levels: [%g %g] dB\n',ampsig(1),ampsig(2));
    end

    % baseall = mean(mean(allamps(base,:)));
    % fprintf('#2 Size of allamps = [%d %d]\n',size(allamps,1),size(allamps,2));

    fprintf('Subtracting the mean baseline log amplitude \n');

    %fprintf('Subtracting the mean baseline log amplitude %g\n',baseall);
    % allamps = allamps./baseall;
    % fprintf('#3 Size of allamps = [%d %d]\n',size(allamps,1),size(allamps,2));
    
    if avewidth > 1 | decfactor > 1
        % Note: using square window
        if Nosort == YES
            fprintf(...
                'Smoothing the amplitude epochs using a window width of %g epochs ',...
                avewidth);
        else % sort trials
            fprintf(...
                'Smoothing the sorted amplitude epochs with a %g-epoch moving window.',...
                avewidth);
        end
        fprintf('\n');
        fprintf('  and a decimation factor of %g\n',decfactor);

        %fprintf('4 Size of allamps = [%d %d]\n',size(allamps,1),size(allamps,2));

        if exist('phargs') % if phase-sorted trials, use circular/wrap-around smoothing
           backhalf  = floor(avewidth/2);
           fronthalf = floor((avewidth-1)/2);
           if avewidth > 2
            [allamps,outtrials] = movav([allamps(:,[(end-backhalf+1):end]),...
                                      allamps,...
                                      allamps(:,[1:fronthalf])],...
                                      [1:(ntrials+backhalf+fronthalf)],avewidth,decfactor); 
            				% Note: sort using square window
            [outsort,outtrials] = movav([sortvar((end-backhalf+1):end),...
                                      sortvar,...
                                      sortvar(1:fronthalf)],...
                                      1:(ntrials+backhalf+fronthalf),avewidth,decfactor); 
                                        % outtrials = 1:ntrials;
            if ~isempty(auxvar)
               [auxvar,tmp] = movav([auxvar(:,[(end-backhalf+1):end]),...
                                      auxvar,...
                                      auxvar(:,[1:fronthalf])],...
                                      [1:(ntrials+backhalf+fronthalf)],avewidth,decfactor); 
            end
           else % avewidth==2
            [allamps,outtrials] = movav([allamps(:,end),allamps],...
                                      [1:(ntrials+1)],avewidth,decfactor); 
            				% Note: sort using square window
            [outsort,outtrials] = movav([sortvar(end) sortvar],...
                                      1:(ntrials+1),avewidth,decfactor); 
                                        % outtrials = 1:ntrials;
            [auxvar,tmp] = movav([auxvar(:,end),auxvar],[1:(ntrials+1)],avewidth,decfactor); 
           end
        else % if trials not phase sorted, no wrap-around
            [allamps,outtrials] = movav(allamps,1:ntrials,avewidth,decfactor); 
                                        % Note: using square window
            %fprintf('5 Size of allamps = [%d %d]\n',size(allamps,1),size(allamps,2));
            [outsort,outtrials] = movav(sortvar,1:ntrials,avewidth,decfactor); 
            if ~isempty(auxvar)
              [auxvar,tmp] = movav(auxvar,1:ntrials,avewidth,decfactor); 
            end
        end
        fprintf('Output allamps data will be %d frames by %d smoothed trials.\n',...
                                      frames,length(outtrials));

    else % if no smoothing
        outtrials = 1:ntrials;
        outsort = sortvar;
    end

    allamps = 20*log10(allamps);
    if isnan(baseamp) % if not specified in 'limits'
        [amps,baseamp] = rmbase(amps,length(times),base); % remove (log) baseline
        allamps = allamps - baseamp; % divide by (non-log) baseline amplitude
    else
        amps = amps-baseamp; % use specified (log) baseamp
        allamps = allamps - baseamp; % = divide by (non-log) baseline amplitude
        if isnan(signifs);
            ampsig = ampsig-baseamp;
        end
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%% Find color axis limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if ~isempty(Caxis) 
        mindat = Caxis(1);
        maxdat = Caxis(2);
        fprintf('Using the specified caxis range of [%g,%g].\n',...
                mindat,maxdat);
    else
        mindat = min(min(allamps));
        maxdat = max(max(allamps));
        maxdat =  max(abs([mindat maxdat])); % make symmetrical about 0
        mindat = -maxdat;
        if ~isempty(caxfraction)
            adjmax = (1-caxfraction)/2*(maxdat-mindat);
            mindat = mindat+adjmax;
            maxdat = maxdat-adjmax;
            fprintf(...
                'The caxis range will be %g times the sym. abs. data range -> [%g,%g].\n',...
                caxfraction,mindat,maxdat);
        else
            fprintf(...
                'The caxis range will be the sym. abs. data range -> [%g,%g].\n',...
                mindat,maxdat);
        end
    end
    %
    %%%%%%%%%%%%%%%%%%%%% Image amplitudes at coherfreq %%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if strcmpi(noshow, 'no')
        fprintf('Plotting amplitudes at freq %g Hz instead of potentials.\n',coherfreq);
        if TIMEX
            imagesc(times,outtrials,allamps',[mindat,maxdat]);% plot time on x-axis
            set(gca,'Ydir','normal');
            axis([timelimits(1) timelimits(2) ...
                  min(outtrials) max(outtrials)]);
        else
            imagesc(outtrials,times,allamps,[mindat,maxdat]); % plot trials on x-axis
            axis([min(outtrials) max(outtrials)...
                  timelimits(1) timelimits(2)]);
        end
        drawnow
        hold on
    end;
    
elseif exist('data2') %%%%%% Plot allcohers instead of data %%%%%%%%%%%%%%%%%%%
                      %%%%%%%%% UNDOCUMENTED AND DEPRECATED OPTION %%%%%%%%%%%%
    if freq > 0 
        coherfreq = mean(freq); % use phase-sort frequency
    end
    if alpha>0
        fprintf('Computing and plotting %g coherence significance level...\n',alpha);
        %[amps,cohers,cohsig,ampsig,allcohers] = ...
        %crosscoher(urdata,data2,length(times),srate,coherfreq,DEFAULT_CYCLES,alpha);
        fprintf('Inter-Trial Coherence significance level: %g\n',cohsig);
        fprintf('Amplitude significance levels: [%g %g]\n',ampsig(1),ampsig(2));
    else
        %[amps,cohers,cohsig,ampsig,allcohers] = ...
        % crosscoher(urdata,data2,length(times),srate,coherfreq,DEFAULT_CYCLES,0);
    end
    if ~exist('allcohers')
        fprintf('erpimage(): allcohers not returned....\n')
        return
    end
    allamps = allcohers; % output variable
                         % fprintf('Size allcohers = (%d, %d)\n',size(allcohers,1),size(allcohers,2));
                         % fprintf('#1 Size of allcohers = [%d %d]\n',size(allcohers,1),size(allcohers,2));
    base = find(times<=0);
    if length(base)<2
        base = 1:floor(length(times)/4); % default first quarter-epoch
    end
    amps = 20*log10(amps); % convert to dB
    ampsig = 20*log10(ampsig); % convert to dB
    if isnan(baseamp)
        [amps,baseamp] = rmbase(amps,length(times),base); % remove baseline
    else
        amps = amps - baseamp;
    end
    % fprintf('#2 Size of allcohers = [%d %d]\n',size(allcohers,1),size(allcohers,2));
    
    if avewidth > 1 | decfactor > 1
        % Note: using square window
        if Nosort == YES
            fprintf(...
                'Smoothing the amplitude epochs using a window width of %g epochs '...
                ,avewidth);
        else
            fprintf(...
                'Smoothing the sorted amplitude epochs with a %g-epoch moving window.'...
                ,avewidth);
        end
        fprintf('\n');
        fprintf('  and a decimation factor of %g\n',decfactor);
        % fprintf('4 Size of allcohers = [%d %d]\n',size(allcohers,1),size(allcohers,2));
        
        [allcohers,outtrials] = movav(allcohers,1:ntrials,avewidth,decfactor); 
        % Note: using square window
        % fprintf('5 Size of allcohers = [%d %d]\n',size(allcohers,1),size(allcohers,2));
        [outsort,outtrials] = movav(sortvar,1:ntrials,avewidth,decfactor); 
        fprintf('Output data will be %d frames by %d smoothed trials.\n',...
                frames,length(outtrials));
    else
        outtrials = 1:ntrials;
        outsort = sortvar;
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%% Find color axis limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if ~isempty(Caxis) 
        mindat = Caxis(1);
        maxdat = Caxis(2);
        fprintf('Using the specified caxis range of [%g,%g].\n',...
                mindat,maxdat);
    else
        mindat = -1;
        maxdat = 1
        fprintf(...
            'The caxis range will be the sym. abs. data range [%g,%g].\n',...
            mindat,maxdat);
    end
    %
    %%%%%%%%%%%%%%%%%%%%% Image coherences at coherfreq %%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if strcmpi(noshow, 'no')
        fprintf('Plotting coherences at freq %g Hz instead of potentials.\n',coherfreq);
        if TIMEX
            imagesc(times,outtrials,allcohers',[mindat,maxdat]);% plot time on x-axis
            set(gca,'Ydir','normal');
            axis([timelimits(1) timelimits(2) ...
                  min(outtrials) max(outtrials)]);
        else
            imagesc(outtrials,times,allcohers,[mindat,maxdat]); % plot trials on x-axis
            axis([min(outtrials) max(outtrials)...
                  timelimits(1) timelimits(2)]);
        end
        drawnow
        hold on
    end;
    
end %%%%%%%%%%%%%%%%%%%%%%%%%%% End plot image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% plot vert lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(verttimes)
 if size(verttimes,1) ~= 1 & size(verttimes,2) == 1 & size(verttimes,1) ~= ntrials
        verttimes = verttimes';
 end
 if size(verttimes,1) ~= 1 & size(verttimes,1) ~= ntrials
    fprintf('\nerpimage(): vert arg matrix must have 1 or %d rows\n',ntrials);
    return
 end;
 if strcmpi(noshow, 'no')
     if size(verttimes,1) == 1
         fprintf('Plotting %d lines at times: ',size(verttimes,2));
     else
         fprintf('Plotting %d traces starting at times: ',size(verttimes,2));
     end
     for vt = verttimes % for each column
         fprintf('%g ',vt(1));
         if isnan(aligntime) % if nor re-aligned data
             if TIMEX          % overplot vt on image
                 if length(vt)==1
                     figure(curfig);plot([vt vt],[0 max(outtrials)],DOTSTYLE,'Linewidth',VERTWIDTH);
                 elseif length(vt)==ntrials
                     [outvt,ix] = movav(vt,1:ntrials,avewidth,decfactor); 
                     figure(curfig);plot(outvt,outtrials,DOTSTYLE,'Linewidth',VERTWIDTH);
                 end
             else
                 if length(vt)==1
                     figure(curfig);plot([0 max(outtrials)],[vt vt],DOTSTYLE,'Linewidth',VERTWIDTH);
                 elseif length(vt)==ntrials
                     [outvt,ix] = movav(vt,1:ntrials,avewidth,decfactor); 
                     figure(curfig);plot(outtrials,outvt,DOTSTYLE,'Linewidth',VERTWIDTH);
                 end
             end
         else                % re-aligned data
             if TIMEX          % overplot vt on image
                 if length(vt)==ntrials
                     [outvt,ix] = movav(vt,1:ntrials,avewidth,decfactor); 
                     figure(curfig);plot(aligntime+outvt-outsort,outtrials,DOTSTYLE,'LineWidth',VERTWIDTH); 
                 elseif length(vt)==1
                     figure(curfig);plot(aligntime+vt-outsort,outtrials,DOTSTYLE,'LineWidth',VERTWIDTH); 
                 end
             else
                 if length(vt)==ntrials
                     [outvt,ix] = movav(vt,1:ntrials,avewidth,decfactor); 
                     figure(curfig);plot(outtrials,aligntime+outvt-outsort,DOTSTYLE,'LineWidth',VERTWIDTH); 
                 elseif length(vt)==1
                     figure(curfig);plot(outtrials,aligntime+vt-outsort,DOTSTYLE,'LineWidth',VERTWIDTH); 
                 end
             end
         end
     end
     %end
     fprintf('\n');
 end;
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% plot horizontal ('horz') lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~isempty(horzepochs)
 if size(horzepochs,1) > 1 & size(horzepochs,1) > 1
    fprintf('\nerpimage(): horz arg must be a vector\n');
    return
 end;
 if strcmpi(noshow, 'no')
     fprintf('Plotting %d lines at epochs: ',length(horzepochs));
     for he = horzepochs % for each horizontal line
         fprintf('%g ',he);
             if TIMEX          % overplot he on image
                 figure(curfig);plot([timelimits(1) timelimits(2)],[he he],LINESTYLE,'Linewidth',HORZWIDTH);
             else
                 figure(curfig);plot([he he], [timelimits(1) timelimits(2)],LINESTYLE,'Linewidth',HORZWIDTH);
             end
     end
     %end
     fprintf('\n');
 end;
end
if strcmpi(noshow, 'no')
    set(gca,'FontSize',TICKFONT)
    hold on;
end;
%
%%%%%%%%%%% plot vertical line at 0 or align time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmpi(noshow, 'no')
    if ~isnan(aligntime) % if trials time-aligned 
        if times(1) <= aligntime & times(frames) >= aligntime
            figure(curfig);plot([aligntime aligntime],[min(outtrials) max(outtrials)],...
                               'k','Linewidth',ZEROWIDTH); % plot vertical line at time 0
            % plot vertical line at aligntime
        end
    else % trials not time-aligned 
        if times(1) <= 0 & times(frames) >= 0
            figure(curfig);plot([0 0],[min(outtrials) max(outtrials)],...
                               'k','Linewidth',ZEROWIDTH); % plot smoothed sortwvar
        end
    end
end;

if Noshow == NO & ( min(outsort) < timelimits(1) ...
                   |max(outsort) > timelimits(2))
  ur_outsort = outsort; % store the pre-adjusted values
  fprintf('Not all sortvar values within time vector limits: \n')
  fprintf('        outliers will be shown at nearest limit.\n');
  i = find(outsort< timelimits(1));
  outsort(i) = timelimits(1);
  i = find(outsort> timelimits(2));
  outsort(i) = timelimits(2);
end

if strcmpi(noshow, 'no')
    if TIMEX
        if Nosort == YES
            figure(curfig);l=ylabel('Trials');
        else
            if exist('phargs')
                figure(curfig);l=ylabel('Phase-sorted Trials');
            elseif exist('ampargs')
                figure(curfig);l=ylabel('Amplitude-sorted Trials');
            else
                l=ylabel('Sorted Trials');
            end
        end
    else % if switch x<->y axes
        if Nosort == YES & NoTimeflag==NO
            figure(curfig);l=xlabel('Trials');
        else
            if exist('phargs')
                figure(curfig);l=ylabel('Phase-sorted Trials');
            elseif NoTimeflag == NO
                figure(curfig);l=xlabel('Sorted Trials');
            end
        end
    end
    set(l,'FontSize',LABELFONT);
    
    t=title(titl);
    set(t,'FontSize',LABELFONT);
    
    set(gca,'Box','off');
    set(gca,'Fontsize',TICKFONT);
    set(gca,'color',BACKCOLOR);
    if Erpflag == NO & NoTimeflag == NO
        if exist('NoTimesPassed')~=1
            figure(curfig);l=xlabel('Time (ms)');
        else
            figure(curfig);l=xlabel('Frames');
        end
        set(l,'Fontsize',LABELFONT);
    end
end;

%
%%%%%%%%%%%%%%%%%%%% Overplot sortvar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmpi(noshow, 'no')

    if Noshow == YES
        fprintf('Not overplotting sorted sortvar on data.\n');
        
    elseif isnan(aligntime) % plot sortvar on un-aligned data
        
        if Nosort == NO;
            fprintf('Overplotting sorted sortvar on data.\n');
        end
        hold on; 
        if TIMEX      % overplot sortvar
            figure(curfig);plot(outsort,outtrials,'k','LineWidth',SORTWIDTH); 
        else
            figure(curfig);plot(outtrials,outsort,'k','LineWidth',SORTWIDTH);
        end                                                 
        drawnow
    else % plot re-aligned zeros on sortvar-aligned data
        if Nosort == NO;
            fprintf('Overplotting sorted sortvar on data.\n');
        end
        hold on; 
        if TIMEX      % overplot re-aligned 0 time on image
            figure(curfig);plot([aligntime aligntime],[min(outtrials) max(outtrials)],...
                                            'k','LineWidth',SORTWIDTH);
        else
            figure(curfig);plot([[min(outtrials) max(outtrials)],aligntime aligntime],...
                                            'k','LineWidth',SORTWIDTH);
        end
        fprintf('Overplotting realigned times-zero on data.\n');
        hold on; 
        
        if TIMEX      % overplot realigned sortvar on image
            figure(curfig);plot(0+aligntime-outsort,outtrials,'k','LineWidth',ZEROWIDTH); 
        else
            figure(curfig);plot(0+outtrials,aligntime-outsort,'k','LineWidth',ZEROWIDTH); 
        end                                                 
        drawnow
    end
end;

%
%%%%%%%%%%%%%%%%%%%% Overplot auxvar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmpi(noshow, 'no')
    if ~isempty(auxvar)
        fprintf('Overplotting auxvar(s) on data.\n');
        hold on; 
        auxtrials = outtrials(:)' ; % make row vector
        
        if exist('auxcolors')~=1 % If no auxcolors specified
            auxcolors = cell(1,size(auxvar,1));
            for c=1:size(auxvar,1)
                auxcolors(c) = {'k'};       % plot auxvars as black trace(s)
            end
        end
        for c=1:size(auxvar,1)
            if isnan(aligntime) % plot auxvar on un-aligned data
                auxcolor = auxcolors{c};
                if TIMEX      % overplot auxvar
                    figure(curfig);plot(auxvar(c,:)',auxtrials',auxcolor,'LineWidth',SORTWIDTH); 
                else
                    figure(curfig);plot(auxtrials',auxvar(c,:)',auxcolor,'LineWidth',SORTWIDTH);
                end                                                 
                drawnow
            else % plot re-aligned zeros on sortvar-aligned data
                auxcolor = auxcolors{c};
                if TIMEX      % overplot realigned 0-time on image
                    figure(curfig);plot(0+aligntime-auxvar(c,:)',auxtrials',auxcolor,'LineWidth',ZEROWIDTH); 
                else
                    figure(curfig);plot(0+auxtrials',aligntime-auxvar(c,:)',auxcolor,'LineWidth',ZEROWIDTH); 
                end                                                 
                drawnow
            end % aligntime
        end % c
    end % auxvar
end;
%
%%%%%%%%%%%%%%%%%%%%%%%% Plot colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmpi(noshow, 'no')
    if Colorbar == YES
        pos=get(ax1,'Position');
        axcb=axes('Position',...
                  [pos(1)+pos(3)+0.02 pos(2) ...
                   0.03 pos(4)]);
        figure(curfig);cbar(axcb,0,[mindat,maxdat]); % plot colorbar to right of image
        set(axcb,'fontsize',TICKFONT);
        % drawnow
        axes(ax1); % reset current axes to the erpimage
    end
end;

%
%%%%%%%%%%%%%%%%%%%%%%% Compute ERP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
erp = [];
if Erpflag == YES 
  if exist('erpalpha')
    [erp erpsig] = nan_mean(urdata',erpalpha);   
    fprintf('   Mean ERP (p<%g) significance threshold: +/-%g\n',erpalpha,mean(erpsig));
  else
    [erp] = nan_mean(urdata');   % compute erp average, ignoring nan's
  end
end;          
if Erpflag == YES & strcmpi(noshow, 'no')
    axes(ax1); % reset current axes to the erpimage
    xtick = get(ax1,'Xtick');     % remember x-axis tick locations
    xticklabel = get(ax1,'Xticklabel');     % remember x-axis tick locations
    set(ax1, 'xticklabel', []);
    widthxticklabel = size(xticklabel,2);
    xticklabel = cellstr(xticklabel);
    for tmpindex = 1:length(xticklabel)
        if length(xticklabel{tmpindex}) < widthxticklabel
            spaces = char(ones(1,ceil((widthxticklabel-length(xticklabel{tmpindex}))/2) )*32);
            xticklabel{tmpindex} = [spaces xticklabel{tmpindex}];
        end;
    end;
    xticklabel = strvcat(xticklabel);
    if Erpstdflag == YES
        stdev = nan_std(urdata');
    end;
    %
    %%%%%% Plot ERP time series below image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if isnan(maxerp)
        fac = 10;
        maxerp = 0;
        while maxerp == 0
            maxerp = round(fac*YEXPAND*max(erp))/fac; % minimal decimal places
            fac = 10*fac;
        end
        if Erpstdflag == YES
            fac = fac/10;
            maxerp = max(maxerp, round(fac*YEXPAND*max(erp+stdev))/fac);
        end;
    end
    if isnan(minerp)
        fac = 1;
        minerp = 0;
        while minerp == 0
            minerp = round(fac*YEXPAND*min(erp))/fac; % minimal decimal places
            fac = 10*fac;
        end
        if Erpstdflag == YES
            fac = fac/10;
            minerp = min(minerp, round(fac*YEXPAND*min(erp-stdev))/fac);
        end;
    end
    limit = [timelimits(1:2) -max(abs([minerp maxerp])) max(abs([minerp maxerp]))];
          
    if ~isnan(coherfreq)
        set(ax1,'Xticklabel',[]);     % remove tick labels from bottom of image
        ax2=axes('Position',...
                 [gcapos(1) gcapos(2)+2/3*image_loy*gcapos(4) ...
                  gcapos(3) (image_loy/3-YGAP)*gcapos(4)]);
    else
        ax2=axes('Position',...
                 [gcapos(1) gcapos(2) ...
                  gcapos(3) image_loy*gcapos(4)]);
    end
    fprintf('Plotting the ERP trace below the ERP image\n');
    if Erpstdflag == YES
        if Showwin
          tmph = plot1trace(ax2,times,erp,limit, [], stdev,[],times(winlocs)); % plot ERP +/-stdev
        else
          tmph = plot1trace(ax2,times,erp,limit, [], stdev,[],[]); % plot ERP +/-stdev
        end
    elseif ~isempty('erpsig')
        erpsig = [erpsig;-1*erpsig];
        if Showwin
          tmph = plot1trace(ax2,times,erp,limit,erpsig,[],times(winlocs)); % plot ERP and 0+/-alpha threshold
        else
          tmph = plot1trace(ax2,times,erp,limit,erpsig,[],[]); % plot ERP and 0+/-alpha threshold
        end
    else
        if Showwin
          tmph = plot1trace(ax2,times,erp,limit,[],[],times(winlocs)); % plot ERP alone
        else
          tmph = plot1trace(ax2,times,erp,limit,[],[],[]); % plot ERP alone
        end
    end;
    
    if ~isnan(aligntime)
        line([aligntime aligntime],[limit(3:4)*1.1],'Color','k','LineWidth',ZEROWIDTH); % x=median sort value
        line([0 0],[limit(3:4)*1.1],'Color','k','LineWidth',ZEROWIDTH); % x=median sort value
        % remove y axis
        if length(tmph) > 1
            delete(tmph(end));
        end;
    end
    
    set(ax2,'Xtick',xtick);        % use same Xticks as erpimage above
    if ~isnan(coherfreq)
        set(ax2,'Xticklabel',[]);    % remove tick labels from ERP x-axis
    else % bottom axis
        set(ax2,'Xticklabel',xticklabel); % add ticklabels to ERP x-axis
    end
    
    set(ax2,'Yticklabel',[]);      % remove tick labels from left of image
    set(ax2,'YColor',BACKCOLOR);
    if isnan(coherfreq)            % if no amp and coher plots below . . .
        if TIMEX & NoTimeflag == NO
            if exist('NoTimesPassed')~=1
                figure(curfig);l=xlabel('Time (ms)');
            else
                figure(curfig);l=xlabel('Frames');
            end
            set(l,'FontSize',LABELFONT);
        else
            if exist('NoTimesPassed')~=1
                figure(curfig);l=ylabel('Time (ms)');
            else
                figure(curfig);l=ylabel('Frames');
            end
            set(l,'FontSize',LABELFONT);
        end
    end
    
    if ~isempty(verttimes) 
        if size(verttimes,1) == ntrials
            vts=sort(verttimes); 
            vts = vts(ceil(ntrials/2),:); % plot median verttimes values if a matrix
        else 
            vts = verttimes(:)';  % make verttimes a row vector
        end
        for vt = vts
            if isnan(aligntime)
                if TIMEX      % overplot vt on ERP
                    figure(curfig);plot([vt vt],[limit(3:4)],DOTSTYLE,'Linewidth',VERTWIDTH);
                else
                    figure(curfig);plot([min(outtrials) max(outtrials)],[limit(3:4)],DOTSTYLE,'Linewidth',VERTWIDTH);
                end
            else
                if TIMEX      % overplot realigned vt on ERP
                    figure(curfig);plot(repmat(median(aligntime+vt-outsort),1,2),[limit(3),limit(4)],...
                         DOTSTYLE,'LineWidth',VERTWIDTH); 
                else
                    figure(curfig);plot([limit(3),limit(4)],repmat(median(aligntime+vt-outsort),1,2),...
                         DOTSTYLE,'LineWidth',VERTWIDTH); 
                end                                                 
            end
        end
    end
    
    limit = double(limit);
    ydelta = double(1/10*(limit(2)-limit(1))); 
    ytextoffset = double(limit(1)-1.1*ydelta);
    ynumoffset  = double(limit(1)-0.3*ydelta); % double for Matlab 7
    
    t=text(ynumoffset,0.7*limit(3), num2str(limit(3)));
    set(t,'HorizontalAlignment','right','FontSize',TICKFONT)
    
    t=text(ynumoffset,0.7*limit(4), num2str(limit(4)));
    set(t,'HorizontalAlignment','right','FontSize',TICKFONT)
    
    ynum = 0.7*(limit(3)+limit(4))/2;
    t=text(ytextoffset,ynum,yerplabel,'Rotation',90);
    set(t,'HorizontalAlignment','center','FontSize',LABELFONT)
    
    set(ax2,'Fontsize',TICKFONT);
    set(ax2,'Box','off','color',BACKCOLOR);
    drawnow
end

%
%%%%%%%%%%%%%%%%%%%%% Plot amp, coher time series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~isnan(coherfreq) 
    if freq > 0 
        coherfreq = mean(freq); % use phase-sort frequency
    end
    %
    %%%%%% Plot amp axis below ERP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if ~Allampsflag %%%% don't repeat computation if already done for 'allamps'
        
        fprintf('Computing and plotting amplitude at %g Hz.\n',coherfreq);
        
        if ~isnan(signifs) | Cohsigflag==NO % don't compute or plot signif. levels
            [amps,cohers] = phasecoher(urdata,size(times,2),srate,coherfreq,DEFAULT_CYCLES);
            if ~isnan(signifs) 
                ampsig = signifs([1 2]);
                fprintf('Using supplied amplitude significance levels: [%g,%g]\n',...
                        ampsig(1),ampsig(2));
                cohsig = signifs(3);
                fprintf('Using supplied coherence significance level: %g\n',cohsig);
            end
        else % compute amps, cohers with significance
            fprintf(...
                'Computing and plotting %g coherence significance level at %g Hz...\n',...
                alpha,                          coherfreq);
            [amps,cohers,cohsig,ampsig] = ...
                phasecoher(urdata,size(times,2),srate,coherfreq,DEFAULT_CYCLES,alpha);
            fprintf('Coherence significance level: %g\n',cohsig);
            ampsig = 20*log10(ampsig); % convert to dB
        end
        
        amps   = 20*log10(amps);      % convert to dB
        fprintf('Data amplitude levels: [%g %g] dB\n',min(amps),max(amps));
        if alpha>0 % if computed significance levels
            fprintf('Data amplitude significance levels: [%g %g] dB\n',ampsig(1),ampsig(2));
        end
        
        if isnan(baseamp) % if baseamp not specified in 'limits'
            base = find(times<=DEFAULT_BASELINE_END); % use default baseline end point (ms)
            if length(base)<2
                base = 1:floor(length(times)/4); % default first quarter-epoch
                fprintf('Using %g to %g ms as amplitude baseline.\n',...
                        times(1),times(base(end)));
            end
            [amps,baseamp] = rmbase(amps,length(times),base); % remove (log) baseline
            fprintf('Removed baseline amplitude of %d dB for plotting.\n',baseamp);
        else
            fprintf('Removing specified baseline amplitude of %d dB for plotting.\n',...
                    baseamp);
            amps = amps-baseamp;
        end
        if Cohsigflag
            ampsig = ampsig - baseamp;
        end;
    end % ~Allampsflag
    
    if strcmpi(noshow, 'no')
        axis('off') % rm ERP axes axis and labels
        ax3=axes('Position',...
                 [gcapos(1) gcapos(2)+1/3*image_loy*gcapos(4) ...
                  gcapos(3) (image_loy/3-YGAP)*gcapos(4)]);
        
        if isnan(maxamp) % if not specified
            fac = 1;
            maxamp = 0;
            while maxamp == 0
                maxamp = floor(YEXPAND*fac*max(amps))/fac; % minimal decimal place
                fac = 10*fac;
            end
            maxamp = maxamp + 10/fac;
            if Cohsigflag
                if ampsig(2)>maxamp
                    if ampsig(2)>0
                        maxamp = 1.01*(ampsig(2));
                    else
                        maxamp = 0.99*(ampsig(2));
                    end
                end
            end
        end
        
        if isnan(minamp) % if not specified
            fac = 1;
            minamp = 0;
            while minamp == 0
                minamp = floor(YEXPAND*fac*max(-amps))/fac; % minimal decimal place
                fac = 10*fac;
            end
            minamp = minamp + 10/fac;
            minamp = -minamp;
            if Cohsigflag
                if ampsig(1)< minamp
                    if ampsig(1)<0
                        minamp = 1.01*(ampsig(1));
                    else
                        minamp = 0.99*(ampsig(1));
                    end
                end
            end
        end
        
        fprintf('Plotting the ERSP amplitude trace below the ERP\n');
        fprintf('Min, max plotting amplitudes: [%g, %g] dB\n',minamp,maxamp);
        fprintf('     relative to baseamp: %g dB\n',baseamp);
        if Cohsigflag
                ampsiglims = [repmat(ampsig(1)-mean(ampsig),1,length(times))];
                ampsiglims = [ampsiglims;-1*ampsiglims];
        	plot1trace(ax3,times,amps,[timelimits minamp(1) maxamp(1)],ampsiglims,[],[]); % plot AMP
        else
        	plot1trace(ax3,times,amps,[timelimits minamp(1) maxamp(1)],[],[],[]); % plot AMP
        end
        
        if ~isnan(aligntime)
            line([aligntime aligntime],[minamp(1) maxamp(1)]*1.1,'Color','k'); 
            % x=median sort value
        end
        set(ax3,'Xtick',xtick);
        set(ax3,'Xticklabel',[]);   % remove tick labels from bottom of image
        set(ax3,'Yticklabel',[]);   % remove tick labels from left of image
        set(ax3,'YColor',BACKCOLOR);
        axis('off');
        set(ax3,'Box','off','color',BACKCOLOR);
        
        if ~isempty(verttimes)
            if size(verttimes,1) == ntrials
                vts=sort(verttimes); 
                vts = vts(ceil(ntrials/2),:); % plot median values if a matrix
            else 
                vts=verttimes(:)';	
            end
            for vt = vts
                if isnan(aligntime)
                    if TIMEX      % overplot vt on amp
                        figure(curfig);plot([vt vt],[minamp(1) maxamp(1)],DOTSTYLE,...
                             'Linewidth',VERTWIDTH);
                    else
                        figure(curfig);plot([min(outtrials) max(outtrials)],[minamp(1) maxamp(1)],DOTSTYLE,...
                             'Linewidth',VERTWIDTH);
                    end
                else
                    if TIMEX      % overplot realigned vt on amp
                        figure(curfig);plot(repmat(median(aligntime+vt-outsort),1,2),[minamp(1),maxamp(1)],DOTSTYLE,...
                             'LineWidth',VERTWIDTH); 
                    else
                        figure(curfig);plot([minamp,maxamp],repmat(median(aligntime+vt-outsort),1,2),DOTSTYLE,...
                             'LineWidth',VERTWIDTH); 
                    end                                                 
                end
            end
        end
        
        if 0 % Cohsigflag % plot amplitude significance levels
            hold on
            figure(curfig);plot([timelimits(1) timelimits(2)],[ampsig(1) ampsig(1)] - mean(ampsig),'r',...
                 'linewidth',SIGNIFWIDTH);
            figure(curfig);plot([timelimits(1) timelimits(2)],[ampsig(2) ampsig(2)] - mean(ampsig),'r',...
                 'linewidth',SIGNIFWIDTH);
        end
        
        t=text(ynumoffset,maxamp, num2str(maxamp,3));
        set(t,'HorizontalAlignment','right','FontSize',TICKFONT);
        
        t=text(ynumoffset,0, num2str(0));
        set(t,'HorizontalAlignment','right','FontSize',TICKFONT);
        
        t=text(ynumoffset,minamp, num2str(minamp,3));
        set(t,'HorizontalAlignment','right','FontSize',TICKFONT);
        
        t=text(ytextoffset,(maxamp+minamp)/2,'ERSP','Rotation',90);
        set(t,'HorizontalAlignment','center','FontSize',LABELFONT);
        
        axtmp = axis;
        dbtxt= text(1/13*(axtmp(2)-axtmp(1))+axtmp(1), ...
                    11/13*(axtmp(4)-axtmp(3))+axtmp(3), ...
                    [num2str(baseamp,4) ' dB']);
        set(dbtxt,'fontsize',TICKFONT);
        drawnow;
        set(ax3, 'xlim', timelimits);
        set(ax3, 'ylim', [minamp(1) maxamp(1)]);
   
        %
        %%%%%% Make coher axis below amp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        ax4=axes('Position',...
                 [gcapos(1) gcapos(2) ...
                  gcapos(3) (image_loy/3-YGAP)*gcapos(4)]);
        if isnan(maxcoh)
            fac = 1;
            maxcoh = 0;
            while maxcoh == 0
                maxcoh = floor(YEXPAND*fac*max(cohers))/fac; % minimal decimal place
                fac = 10*fac;
            end
            maxcoh = maxcoh + 10/fac;
            if maxcoh>1
                maxcoh=1; % absolute limit
            end
        end
        if isnan(mincoh)
            mincoh = 0;
        end
        fprintf('Plotting the ITC trace below the ERSP\n');
        if Cohsigflag % plot coherence significance level
            cohsiglims = [repmat(cohsig,1,length(times));zeros(1,length(times))];
            coh_handle = plot1trace(ax4,times,cohers,[timelimits mincoh maxcoh],cohsiglims,[],[]); 
                                                                           % plot COHER, fill sorting window
        else
            coh_handle = plot1trace(ax4,times,cohers,[timelimits mincoh maxcoh],[],[],[]); % plot COHER
        end
        if ~isnan(aligntime)
            line([aligntime aligntime],[[mincoh maxcoh]*1.1],'Color','k'); 
            % x=median sort value
        end
        % set(ax4,'Xticklabel',[]);    % remove tick labels from bottom of image
        set(ax4,'Xtick',xtick);
        set(ax4,'Xticklabel',xticklabel);
        set(ax4,'Ytick',[]);
        set(ax4,'Yticklabel',[]);      % remove tick labels from left of image
        set(ax4,'YColor',BACKCOLOR);
        
        if ~isempty(verttimes)
            if size(verttimes,1) == ntrials
                vts=sort(verttimes); 
                vts = vts(ceil(ntrials/2),:); % plot median values if a matrix
            else 
                vts = verttimes(:)';  % make verttimes a row vector
            end
            for vt = vts
                if isnan(aligntime)
                    if TIMEX      % overplot vt on coher
                        figure(curfig);plot([vt vt],[mincoh maxcoh],DOTSTYLE,'Linewidth',VERTWIDTH);
                    else
                        figure(curfig);plot([min(outtrials) max(outtrials)],[mincoh maxcoh],DOTSTYLE,'Linewidth',VERTWIDTH);
                    end
                else
                    if TIMEX      % overplot realigned vt on coher
                        figure(curfig);plot(repmat(median(aligntime+vt-outsort),1,2),[mincoh,maxcoh],DOTSTYLE,'LineWidth',VERTWIDTH); 
                    else
                        figure(curfig);plot([mincoh,maxcoh],repmat(median(aligntime+vt-outsort),1,2),DOTSTYLE,'LineWidth',VERTWIDTH); 
                    end                                                 
                end
            end
        end
        
        t=text(ynumoffset,0, num2str(0));
        set(t,'HorizontalAlignment','right','FontSize',TICKFONT);
        
        t=text(ynumoffset,maxcoh, num2str(maxcoh));
        set(t,'HorizontalAlignment','right','FontSize',TICKFONT);
        
        t=text(ytextoffset,maxcoh/2,'ITC','Rotation',90);
        set(t,'HorizontalAlignment','center','FontSize',LABELFONT);
        drawnow
        
        %if Cohsigflag % plot coherence significance level
            %hold on
            %plot([timelimits(1) timelimits(2)],[cohsig cohsig],'r',...
                 %'linewidth',SIGNIFWIDTH);
        %end
        
        set(ax4,'Box','off','color',BACKCOLOR);
        set(ax4,'Fontsize',TICKFONT);
        if NoTimeflag==NO
            if exist('NoTimesPassed')~=1
                figure(curfig);l=xlabel('Time (ms)');
            else
                figure(curfig);l=xlabel('Frames');
            end
            set(l,'Fontsize',LABELFONT);
        end
        axtmp = axis;
        hztxt=text(10/13*(axtmp(2)-axtmp(1))+axtmp(1), ...
                   8/13*(axtmp(4)-axtmp(3))+axtmp(3), ...
                   [num2str(coherfreq,4) ' Hz']);
        set(hztxt,'fontsize',TICKFONT);
    end; % noshow
else
    amps   = [];    % null outputs unless coherfreq specified
    cohers = [];
end
axhndls = [ax1 axcb ax2 ax3 ax4];
if exist('ur_outsort')
    outsort = ur_outsort; % restore outsort clipped values, if any
end
if nargout<1
    data = []; % don't spew out data if no args out and no ;
end

%   
%%%%%%%%%%%%%%% plot a topoplot() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
if (~isempty(topomap)) & strcmpi(noshow, 'no') 
    h(12)=axes('Position',...
               [gcapos(1)+0.10*gcapos(3) gcapos(2)+0.92*gcapos(4),...
                0.20*gcapos(3) 0.14*gcapos(4)]);
    % h(12) = subplot('Position',[.10 .86 .20 .14]); 
    fprintf('Plotting a topo map in upper left.\n');
	if length(topomap) == 1
		topoplot(topomap,eloc_file,'electrodes','off', ...
				 'style', 'blank', 'emarkersize1chan', 10);
	else
		topoplot(topomap,eloc_file,'electrodes','off');
	end;
    axis('square')
    axhndls = [axhndls h(12)];
end 

%   
%%%%%%%%%%%%%%% plot a spectrum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
SPECFONT = 10;
if (~isempty(lospecHz)) & strcmpi(noshow, 'no')  
    h(13)=axes('Position',...
               [gcapos(1)+0.82*gcapos(3) ...
                gcapos(2)+0.96*gcapos(4),...
                0.15*gcapos(3)*(0.8/gcapos(3))^0.5 ...
                0.10*gcapos(4)*(0.8/gcapos(4))^0.5]);
    
    % h(13) = subplot('Position',[.75 .88 .15 .10]); 
    fprintf('Plotting the data spectrum in upper right.\n');
    winlength = frames;
    if winlength > 512
        for k=2:5
            if rem(winlength,k) == 0
                break
            end
        end
        winlength = winlength/k;
    end
    
    % [Pxx, Pxxc, F] = PSD(X,NFFT,Fs,WINDOW,NOVERLAP,P)
    if exist('psd') == 2
        [Pxx,F] = psd(reshape(urdata,1,size(urdata,1)*size(urdata,2)),...
                           512,srate,winlength,0,0.05);
    else
        [Pxx,F] = spec(reshape(urdata,1,size(urdata,1)*size(urdata,2)),...
                           512,srate,winlength,0);
    end;
    figure(curfig);plot(F,10*log10(Pxx));
    goodfs = find(F>= lospecHz & F <= hispecHz);
    maxgfs = max(10*log10(Pxx(goodfs)));
    mingfs = min(10*log10(Pxx(goodfs)));
    axis('square')
    axis([lospecHz hispecHz mingfs-1 maxgfs+1]);
    set(h(13),'Box','off','color',BACKCOLOR);
    set(h(13),'Fontsize',SPECFONT);
    figure(curfig);l=ylabel('dB');
    set(l,'Fontsize',SPECFONT);
    if ~isnan(coherfreq)
        hold on; figure(curfig);plot([coherfreq,coherfreq],[mingfs maxgfs],'r');
    end
    axhndls = [axhndls h(13)];
end 

%   
%%%%%%%%%%%%%%% save plotting limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
limits = [min(times) max(times) minerp maxerp minamp maxamp mincoh maxcoh];
limits = [limits baseamp coherfreq];  % add coherfreq to output limits array

%   
%%%%%%%%%%%%%%% turn on axcopy() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmpi(noshow, 'no')  
    axcopy(gcf);
    % eegstr = 'img=get(gca,''''children''''); if (strcmp(img(end),''''type''''),''''image''''), img=get(img(end),''''CData''''); times=get(img(end),''''Xdata''''); clf; args = [''''limits'''' '''','''' times(1) '''','''' times(end)]; if exist(''''EEG'''')==1, args = [args '''','''' ''''srate'''' '''','''' EEG.srate]; end eegplot(img,args); end';
    % axcopy(gcf,eegstr);
end;

fprintf('Done.\n');

%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End erpimage() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
if strcmpi(noshow, 'no')  
    axes('position',gcapos);
    axis off
end;

return
%
%%%%%%%%%%%%%%%%%%% function plot1trace() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [plot_handle] = plot1trace(ax,times,erp,axlimits,signif,stdev,winlocs)
%                           If signif is [], plot erp +/- stdev
%                           Else if signif, plot erp and signif(1,:)&signif(2,:) fill
%                           Else, plot erp alone
%                           If winlocs not [], plot grey back image(s) in sort window
%                                       winlocs(1,1)-> winlocs(1,end) (ms)
%                                        ...
%                                       winlocs(end,1)-> winlocs(end,end) (ms)
  FILLCOLOR    = [.66 .76 1];
  WINFILLCOLOR    = [.88 .92 1];
  ERPDATAWIDTH = 2;
  ERPZEROWIDTH = 2;
  axes(ax);
  if ~isempty(winlocs)
   for k=1:size(winlocs,1)
    winloc = winlocs(k,:);
    fillwinx = [winloc winloc(end:-1:1)];
    hannwin = makehanning(length(winloc));
    hannwin = hannwin./max(hannwin); % make max = 1
    hannwin = hannwin(:)'; % make row vector
    if ~isempty(axlimits) & sum(isnan(axlimits))==0
       % fillwiny = [repmat(axlimits(3),1,length(winloc)) repmat(axlimits(4),1,length(winloc))];
       fillwiny = [hannwin*axlimits(3) hannwin*axlimits(4)];
    else
       % fillwiny = [repmat(min(erp)*1.1,1,length(winloc)) repmat(max(erp)*1.1,1,length(winloc))];
       fillwiny = [hannwin*2*min(erp) hannwin*2*max(erp)];
    end
    fillwh = fill(fillwinx,fillwiny, WINFILLCOLOR); hold on    % plot 0+alpha
    set(fillwh,'edgecolor',WINFILLCOLOR-[.00 .00 0]); % make edges NOT highlighted
   end
  end
  if ~isempty(signif);% (2,times) array giving upper and lower signif limits
      filltimes = [times times(end:-1:1)];
      if size(signif,1) ~=2 | size(signif,2) ~= length(times)
         fprintf('plot1trace(): signif array must be size (2,frames)\n')
         return
      end
      fillsignif = [signif(1,:) signif(2,end:-1:1)];
      fillh = fill(filltimes,fillsignif, FILLCOLOR); hold on    % plot 0+alpha
      set(fillh,'edgecolor',FILLCOLOR-[.02 .02 0]); % make edges slightly highlighted
      % [plot_handle] = plot(times,signif, 'r','LineWidth',1); hold on    % plot 0+alpha
      % [plot_handle] = plot(times,-1*signif, 'r','LineWidth',1); hold on % plot 0-alpha
  end
  if ~isempty(stdev)
      [plot_handle] = plot(times,erp+stdev, 'r--','LineWidth',1); hold on % plot erp+stdev
      [plot_handle] = plot(times,erp-stdev, 'r--','LineWidth',1); hold on % plot erp-stdev
  end
  [plot_handle] = plot(times,erp,'LineWidth',ERPDATAWIDTH); hold on
  if ~isempty(axlimits) & sum(isnan(axlimits))==0
    if axlimits(2)>axlimits(1) & axlimits(4)>axlimits(3)
      axis([axlimits(1:2) 1.1*axlimits(3:4)])
    end
    l1=line([axlimits(1:2)],[0 0],    'Color','k',...
                 'linewidth',ERPZEROWIDTH); % y=zero-line
    l2=line([0 0],[axlimits(3:4)*1.1],'Color','k',...
                 'linewidth',ERPZEROWIDTH); % x=zero-line
    plot_handle = [plot_handle l1 l2];
  end
%
%%%%%%%%%%%%%%%%%%% function phasedet() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% phasedet() - function used in erpimage.m
%              Constructs a complex filter at frequency freq
%

function [ang,amp,win] = phasedet(data,frames,srate,nwin,freq)
% Typical values:
%   frames = 768;
%   srate = 256;
%   nwin = [200:300];
%   freq = 10;
    
data = reshape(data,[frames prod(size(data))/frames]);
% number of cycles depend on window size 
% number of cycles automatically reduced if smaller window
% note: as the number of cycle changes, the frequency shifts a little
%       this has to be fixed
win = exp(2i*pi*freq(:)*[1:length(nwin)]/srate);
win = win .* repmat(makehanning(length(nwin))',length(freq),1);
%tmp =gcf; figure; plot(real(win)); figure(tmp);
%fprintf('ANY NAN ************************* %d\n', any(any(isnan( data(nwin,:)))));

tmpdata = data(nwin,:) - repmat(mean(data(nwin,:), 2), [1 size(data,2)]);
resp = win * tmpdata;
ang = angle(resp);
amp = abs(resp);

%
%%%%%%%%%%%%%% function prctle() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function prctl = prctle(data,pc); % return percentile of a distribution
[prows pcols] = size(pc);
if prows ~= 1 & pcols ~= 1
    error('pc must be a scalar or a vector.');
end
if any(pc > 100) | any(pc < 0)
    error('pc must be between 0 and 100');
end
[i,j] = size(data);
sortdata = sort(data);
if i==1 | j==1 % if data is a vector
  i = max(i,j); j = 1;
  if i == 1,
    fprintf('  prctle() note: input data is a single scalar!\n')
    y = data*ones(length(pc),1); % if data is scalar, return it
    return;
  end
  sortdata = sortdata(:);
end
pt = [0 100*((1:i)-0.5)./i 100];
sortdata = [min(data); sortdata; max(data)];
prctl = interp1(pt,sortdata,pc);
%
%%%%%%%%%%%%%%%%%%%%%%% function nan_mean() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% nan_mean() - Take the column means of a matrix, ignoring NaN values.
%              Return significance bounds if alpha (0 < alpha< <1) is given.
% 
function [out, outalpha]  = nan_mean(in,alpha)
   NBOOT = 500;

   intrials = size(in,1);
   inframes = size(in,2);
   nans = find(isnan(in));
   in(nans) = 0;
   sums = sum(in);
   nonnans = ones(size(in));
   nonnans(nans) = 0;
   nonnans = sum(nonnans);
   nononnans = find(nonnans==0);
   nonnans(nononnans) = 1;
   out = sum(in)./nonnans;
   outalpha = [];
   if nargin>1 
     if NBOOT < round(3/alpha)
        NBOOT = round(3/alpha);
     end
     fprintf('Computing %d bootstrap ERP values... ',NBOOT);
     booterps = zeros(NBOOT,inframes);
     for n=1:NBOOT
         signs = sign(randn(1,intrials)'-0.5);
         booterps(n,:) = sum(repmat(signs,1,inframes).*in)./nonnans;
         if ~rem(n,50)
             fprintf('%d ',n);
         end
     end
     fprintf('\n');
     booterps = sort(abs(booterps));
     alpha = 1+floor(2*alpha*NBOOT); % one-sided probability threshold
     outalpha = booterps(end+1-alpha,:);
   end
   out(nononnans) = NaN;

function out = nan_std(in)
    
    nans = find(isnan(in));
    in(nans) = 0;
   
    nonnans = ones(size(in));
    nonnans(nans) = 0;
    nonnans = sum(nonnans);
    nononnans = find(nonnans==0);
    nonnans(nononnans) = 1;
   
    out = sqrt((sum(in.^2)-sum(in).^2./nonnans)./(nonnans-1));
    out(nononnans) = NaN;

% symmetric hanning function
function w = makehanning(n)
if ~rem(n,2)
   w = 0.5*(1 - cos(2*pi*(1:n/2)'/(n+1)));
   w = [w; w(end:-1:1)];
else
   w = 0.5*(1 - cos(2*pi*(1:(n+1)/2)'/(n+1)));
   w = [w; w(end-1:-1:1)];
end

