% erpimage() - Image single-trial ERPs optionally sorted on and/or aligned to 
%              an input variable and smoothed by moving-average (Note: to
%              return event-aligned data without plotting, use eventlock()).
%              Click on axes to examine separately and zoom.
% Usage:
%   >> [outdata,outvar,outtrials,limits,axhndls,erp, ...
%         amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx] ...
%            ...  
%             = erpimage(data,sortvar,times,'title',avewidth,decimate,...
%                             flag1,arg1,flag2,arg2,...);
% Inputs:
%   data     - single-channel data: format (1,frames*trials) or (frames,trials)
%   sortvar  - vector variable to sort trials on (ntrials = length(sortvar))
%              for example, >> sortvar = rts (in ms)
%   times    - vector of times in ms (frames=length(times)){def|0->[0:frames-1]}
%              else  [startms ntimes srate] = start time (ms), time points per epoch,
%                sampling rate (Hz),
%  'title'   - title string {default none}
%   avewidth - ntrials in moving average window (may be non-int) {def|0->1}
%   decimate - factor to decimate ntrials out by (may be non-int) {def|0->1}
%                If this is large (>sqrt(num. trials)), this many trials output.
% Options:
%   'renorm' - ['yes'|'no'|'formula(x)'] normalize sorting variable. 'yes'=auto.
%              Ex for formula(x): '3*x+2'. Default is 'no'.
%   'align'  - [time] -> time lock data to sortvar aligned to time in msec
%              (time=Inf -> align to median sortvar) {default: no align}
%   'nosort' - don't sort data on sortvar {default: sort}
%   'noplot' - don't plot sortvar {default: plot if in times range}
%   'limits' - [lotime hitime minerp maxerp loamp hiamp locoher hicoher bamp]
%              Can use NaN for missing items and omit late items; use
%              bamp to fix baseline amplitude.
%   'caxis'  - [lo hi] -> set color axis limits {default: data bounds}
%                else [fraction] = set caxis limits at (+/-)fraction*max(abs(data))
%   'cbar'   - plot color bar to right of erp-image {default no}
%   'erp'    - plot erp time average of the trials below the image
%   'auxvar' - [matrix] -> plot auxiliary variable(s) for each trial as separate
%              traces. To plot N traces, the auxvar matrix should be size (N,frames) 
%              ELSE, 'auxvar',{[matrix],{colorstrings}} specifies the N trace colors. 
%              e.g. colorstrings = {'r','bo-','k:'} 
%
% Note:
%       FOR MORE INPUT ARGS: phase,coher,allamps,topo,spec,srate,signif,vert,noxlabel
%                       SEE: >> erpimage moreargs 
% Outputs:
%   outdata   = (times,epochsout) data matrix (after smoothing)
%   outvar    = (1,epochsout)  sortvar vector (after smoothing)
%   outtrials = (1,epochsout)  smoothed trial numbers 
%   limits    = (1,10) array, 1-9 as in 'limits' above, then analysis frequency (Hz)
%   axhndls   = vector of 1-7 plot axes handles (img,cbar,erp,amp,coh,topo,spec)
%   erp       = plotted ERP average
%
% Note: 
%      FOR MORE OUTPUT ARGS: amps,coher,cohsig,ampsig,outamps,phsangls,phsamp,sortidx
%                  SEE:  >> erpimage moreargs
%
% Authors: Scott Makeig, Tzyy-Ping Jung & Arnaud Delorme, 
%          CNL/Salk Institute, La Jolla, 3-2-1998 
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
% Uses included functions:         plot1erp(), phasedet()

% UNIMPLEMENTED - 'allcohers',[data2] -> image the coherences at each time & trial. 
%                   Requires arg 'coher' with alpha significance. 
%                   Shows projection on grand mean coherence vector at each time 
%                   and trial. {default: no}
 
% $Log: not supported by cvs2svn $
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
% 'limits', [lotime hitime] does not work with 'erp'
% 'limits', [... loerp hierp] (still??) may leave "ghost" grey numbers 
%       on the coher axis when printed (-djpeg or -depsc)
% 'allcohers' - not fully implemented, and has been dropped from the help msg

function [data,outsort,outtrials,limits,axhndls,erp,amps,cohers,cohsig,ampsig,allamps,phaseangles,phsamp,sortidx] = erpimage(data,sortvar,times,titl,avewidth,decfactor,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24,arg25,arg26)

%
%%%%%%%%%%%%%%%%%%% Define defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize optional output variables:
erp = []; amps = []; cohers = []; cohsig = []; ampsig = []; 
allamps = []; phaseangles = []; phsamp = []; sortidx = [];
auxvar = [];

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
SIGNIFWIDTH = 1.9;  % Linewidth of red significance lines for amp, coher
DOTSTYLE   = 'k--'; % line style to use for vetical dotted/dashed lines
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
Alignflag = NO;     % don't align data to sortvar by default
Colorbar  = NO;     % if YES, plot a colorbar to right of erp image
Limitflag = NO;     % plot whole times range by default
Phaseflag = NO;     % don't sort by phase
Ampflag   = NO;     % don't sort by amplitude
Valflag   = NO;     % don't sort by value
Srateflag = NO;     % srate not given
Vertflag  = NO;
Renormflag = NO;
yerplabel = '\muV';
yerplabelflag = NO;
verttimes = [];
NoTimeflag= NO;     % by default DO print "Time (ms)" below bottom axis
Signifflag= NO;     % compute significance instead of receiving it
Auxvarflag= NO;
signifs   = NaN;
coherfreq = nan;    % amp/coher-calculating frequency
freq = 0;           % phase-sorting frequency
srate = DEFAULT_SRATE; % from icadefs.m
aligntime = nan;
timelimits= nan;
topomap   = [];     % topo map vector
lospecHz  = [];     % spec lo frequency
topphase = 180;     % default top phase for 'phase' option
renorm    ='no';

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

if isstr(data) 
   ans=strcmp(data,'moreargs');
   if ans==1
     more on
     help erpimopt
     more off
     return
   else
     more on
     help erpimage
     more off
     return
   end
end

if nargin<2
  if size(data,1)==1 | size(data,2)==1
   fprintf('erpimage(): either specify times vector or size-(frames,trials) data.\n')
   return
  end
  times = 1:size(data,1);
  NoTimesPassed= 1;
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
    '\nerpimage(); length of sortvar doesnt divide no. of data elements.\n')
  return
end

if nargin < 6
  decfactor = 0;
end
if nargin < 5
  avewidth = 0;
end
if nargin<4
  titl = '';
end
if nargin<3
  times = NO;
end
if length(times) == 1 | times == NO,
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
   % could use default srate read from icadefs here...
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
if decfactor < 1
  help erpimage
  fprintf('\nerpimage(): Variable decfactor cannot be < 1.\n')
  return
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
	  elseif Alignflag == YES
		  if length(Arg) ~= 1 
			  help erpimage
			  fprintf('\nerpimage(): align arg must be a scalar msec.\n');
			  return
		  end
		  aligntime = Arg(1);
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
              error('erpimage(): Invalid negative argument for keyword ''phasesort''');
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
	  elseif Ampflag == YES
          n = length(Arg);
          if n > 4
			  error('erpimage(): Too many arguments for keyword ''ampsort''');
          end
          ampargs = Arg;
		  
          if ampargs(3) < 0
              error('erpimage(): Invalid negative argument for keyword ''ampsort''');
          end
          if n>=4
			  		if ampargs(4) < 0
				  		error('erpimage(): Invalid negative argument for keyword ''ampsort''');
			  		end
          end
          
          if min(ampargs(1)) < times(1) | max(ampargs(1)) > times(end)
			  		error('erpimage(): time for amplitude sorting filter out of bounds.');
          end

          if ampargs(2) >= 100 | ampargs(2) < -100
			  		error('%-argument for keyword ''ampsort'' must be (-100;100)');
          end
          
          if length(ampargs) == 4 & ampargs(3) > ampargs(4)
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
          if n >= 1 & valargs(1) > valargs(2)
			  		error('erpimage(): Value sorting time range must be increasing.');
          end
          if n==3 & (~isnum(valargs(3)) | valargs(3)==0)
			  		error('erpimage(): Value sorting direction must be +1 or -1.');
          end
          Valflag = NO;
	  elseif strcmp(Arg,'nosort')
		  Nosort = YES;
	  elseif strcmp(Arg,'renorm')
		  Renormflag = YES;
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
	  elseif strcmp(Arg,'valsort')
		  Valflag = YES;
	  elseif strcmp(Arg,'auxvar')
		  Auxvarflag = YES;
	  elseif strcmp(Arg,'yerplabel')
		  yerplabelflag = YES;
	  elseif strcmp(Arg,'srate')
		  Srateflag = YES;
	  elseif strcmp(Arg,'vert') |  strcmp(Arg,'verttimes')
		  Vertflag = YES;
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
	whos auxvar
	if size(auxvar,1)>size(auxvar,2)  % make (N,frames)
		auxvar = auxvar';               % (assuming N < frames)
	end
	if size(auxvar,2) ~= frames
		fprintf('erpimage(): auxvar size should be (N,nframes), e.g., (N,%d)\n',frames);
		return
	end
	if exist('auxcolors')==YES % if specified
		if isa(auxcolors,'cell')==NO
			fprintf('erpimage(): auxcolors argument to auxvar flag must be a cell array. See help.\n');
			return
		end
	end
end
if exist('phargs')
	if phargs(3) > srate/2
		fprintf('erpimage(): Phase-sorting frequency must be less than Nyquist rate.');
	end
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
	if ampargs(3) > srate/2
		fprintf('erpimage(): amplitude-sorting frequency must be less than Nyquist rate.');
	end
	if frames < DEFAULT_CYCLES*srate/ampargs(3)
		fprintf('\nerpimage(): amplitude-sorting freq. (%g) too low: epoch length < %d cycles.\n',...
				ampargs(3),DEFAULT_CYCLES);
		return
	end
	if length(ampargs)==4 & ampargs(4) > srate/2
		ampargs(4) = srate/2;
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
	if frames < DEFAULT_CYCLES*srate/coherfreq(1)
		fprintf('\nerpimage(): coher freq. (%g) too low:  epoch length < %d cycles.\n',...
				coherfreq(1),DEFAULT_CYCLES);
		return
	end
end
          
if isnan(timelimits)
   timelimits = [min(times) max(times)];
end
if ~isnan(aligntime)
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
fprintf('Plotting input data as %d epochs of %d frames sampled at %3.1f Hz.\n',...
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
	fprintf('Removing %d trials having NaN values for sorting variable', length(nanlocs));
	data(:,nanlocs) = [];
	sortvar(nanlocs) = [];
	if exist('data2') == 1
		data2(:,nanlocs) = [];
	end;
	if ~isempty(auxvar)
		auxvar(nanlocs,:) = [];
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
  locx = findstr('x', lower(renorm))
  if length(locx) ~= 1, error('erpimage: unrecognize renormalazing formula'); end;
  eval( [ 'sortvar =' renorm(1:locx-1) 'sortvar' renorm(locx+1:end) ';'] );
end;
%
%%%%%%%%%%%%%%%%%%% Align data to sortvar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~isnan(aligntime)
  if isinf(aligntime)
    aligntime= median(sortvar);
    fprintf('Aligning data to median sortvar.\n'); 
    % Alternative below: trimmed median - ignore top/bottom 5%
    %   ssv = sort(sortvar); % ssv = 'sorted sortvar'
    %   aligntime= median(ssv(ceil(ntrials/20)):floor(19*ntrials/20)); 
  end
  fprintf('Realigned sortvar plotted at %g ms.\n',aligntime);

  aligndata=0*ones(frames,ntrials); % begin with matrix of zeros()
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
end 
%
%%%%%%%%%%%%%%% Sort the data trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if exist('phargs') == 1 % if phase-sort
	if length(phargs) >= 4 % find max frequency in specified band
		[pxx,freqs] = psd(data(:),1024,srate,frames,0);
		
		%gf = gcf;
		% figure;plot(freqs,pxx);
		%xx=axis;
		%axis([phargs(3) phargs(4) xx(3) xx(4)]);
		%figure(gf);
		
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
	
	[dummy minx] = min(abs(times-phargs(1)));
	winlen = floor(3*srate/freq);
	%winloc = minx-[winlen:-1:0]; % ending time version
	winloc = minx-linspace(floor(winlen/2), floor(-winlen/2), winlen+1);
	winloc = winloc(find(winloc>0 & winloc<=frames));
	
	[phaseangles phsamp] = phasedet(data,frames,srate,winloc,freq);
	
	fprintf('Sorting data epochs by phase at %.2f Hz, window centered at %f3. ms.\n',...  
			freq,phargs(1));
	fprintf('Phase is computed using a filter of length %d frames.\n',...
			length(winloc));
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
	
	fprintf('Size of data = [%d,%d]\n',size(data,1),size(data,2));
	sortidx = n(sortidx); % return original trial indices in final sorted order
	if ~isempty(auxvar)
		auxvar = auxvar(:,sortidx);
	end

if exist('ampargs') == 1 % if amplitude-sort
	if length(ampargs) == 4 % find max frequency in specified band
		[pxx,freqs] = psd(data(:),1024,srate,frames,0);
		
		pxx = 10*log10(pxx);
		n = find(freqs >= ampargs(3) & freqs <= ampargs(4));
		if ~length(n)
			freq = ampargs(3);
		end
		[dummy maxx] = max(pxx(n));
		freq = freqs(n(maxx));
	else
		freq = ampargs(3); % else use specified frequency
	end
	
	[dummy minx] = min(abs(times-ampargs(1)));
	winlen = floor(3*srate/freq);
	%winloc = minx-[winlen:-1:0]; % ending time version
	winloc = minx-linspace(floor(winlen/2), floor(-winlen/2), winlen+1);
	winloc = winloc(find(winloc>0 & winloc<=frames));
	
	[phaseangles phsamp] = phasedet(data,frames,srate,winloc,freq);
	
	fprintf('Sorting data epochs by amplitude at %.2f Hz in window centered at %f3. ms.\n',...  
			freq,ampargs(1));
	fprintf('Amplitude is computed using a filter of length %d frames.\n',...
			length(winloc));
	%
	% Reject small (or large) phsamp trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	ampargs(2) = ampargs(2)/100; % convert rejection rate from % to fraction
	[tmp n] = sort(phsamp); % sort amplitudes
	if ampargs(2)>=0
		n = n(ceil(ampargs(2)*length(n))+1:end); % if rej 0, select all trials
		fprintf('Retaining %d epochs (%g percent) with largest power at the analysis frequency,\n',...
			length(n),100*(1-ampargs(2)));
		
	else % ampargs(2) < 0
		ampargs(2) = 1+ampargs(2); % subtract from end
		n = n(1:floor(ampargs(2)*length(n)));
		fprintf('Retaining %d epochs (%g percent) with smallest power at the analysis frequency,\n',...
                      length(n),ampargs(2)*100);
	end
	%
	% Remove low|high-amplitude trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	data = data(:,n); % amp-sort the data, removing rejected-amp trials
	phsamp = phsamp(n);           % amp-sort the amps
	phaseangles = phaseangles(n); % amp-sort the phaseangles
	sortvar = sortvar(n);         % amp-sort the trial indices
	ntrials = length(n);          % number of trials retained
	%
	% Sort remaining data by amplitude %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	[phsamp sortidx] = sort(phsamp); % sort trials on amplitude
	data    =  data(:,sortidx);                % sort data by amp
	phaseangles  =  phaseangles(sortidx);      % sort angles by amp
	sortvar = sortvar(sortidx);                % sort input sortvar by amp
	
	fprintf('Size of data = [%d,%d]\n',size(data,1),size(data,2));
	sortidx = n(sortidx); % return original trial indices in final sorted order
	if ~isempty(auxvar)
		auxvar = auxvar(:,sortidx);
	end

elseif Nosort == YES
  fprintf('Not sorting data on input sortvar.\n');
  sortidx = 1:ntrials;	
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
%
%%%%%%%%%%%%%%%%%%%%%% Sort trials on sortvar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
else 
  [sttime stframe] = min(times-valargs(1));
  if length(valargs)>1
     [endtime endframe] = min(times-valargs(2));
  else
     endframe = stframe;
  end
  if length(valargs)==1
     fprintf('Sorting data on value at time %f ms.\n',sttime);
  elseif length(valargs)>1
     fprintf('Sorting data on mean value between %f and %f ms.\n',...
            sttime,endtime);
  end
  sortval = mean(data(stframe:endframe,:));
  [sortval,sortidx] = sort(sortval);
  data = data(:,sortidx);
  if ~isempty(auxvar)
    auxvar = auxvar(:,sortidx);
  end
end

if max(sortvar)<0
   fprintf('Changing the sign of sortvar: making it positive.\n');
   sortvar = -sortvar;
end
%
%%%%%%%%%%%%%%%%%%% Adjust decfactor if phargs or ampargs %%%%%%%%%%%%%%%%%%%%%
%
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
  decfactor = length(n)/decfactor;
end
% 
%%%%%%%%%%%%%%%%%% Smooth data using moving average %%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% if ~isnan(coherfreq)
  urdata = data; % compute amp, coher on unsmoothed data
% end
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
  [data,outtrials]    = movav(data,1:ntrials,avewidth,decfactor); 
                        % Note: movav here sorts using square window
  [outsort,outtrials] = movav(sortvar,1:ntrials,avewidth,decfactor); 
  if ~isempty(auxvar)
         [auxvar,tmp] = movav(auxvar,1:ntrials,avewidth,decfactor); 
  end
  fprintf('Output data will be %d frames by %d smoothed trials.\n',...
                          frames,length(outtrials));
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
  fprintf('Using the specified caxis range of [%g,%g].\n',...
                                           mindat,maxdat);
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
	[pxx,tmpfreq] = psd(data(:),1024,srate,frames,0);

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

if ~Allampsflag & ~exist('data2') %%%%%%%%%%%%%%% Plot ERP image %%%%%%%%%%%%%%

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
elseif Allampsflag %%%%%%%%%%%%%%%% Plot allamps instead of data %%%%%%%%%%%%%%

 if freq > 0 
    coherfreq = freq; % use phase-sort frequency
 end

 if ~isnan(signifs) % plot received significance levels
   fprintf('Computing and plotting received amp and ITC significance levels...\n');
   [amps,cohers,cohsig,ampsig,allamps] = ...
     phasecoher(urdata,length(times),srate,coherfreq,DEFAULT_CYCLES,0);
     % need to receive cohsig and ampsig to get allamps
   ampsig = signifs([1 2]); % assume these already in dB
   cohsig = signifs(3);

 elseif alpha>0 % compute significance levels
   fprintf('Computing and plotting %g amp and ITC significance level...\n',alpha);
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
   else
    fprintf(...
       'Smoothing the sorted amplitude epochs with a %g-epoch moving window.',...
            avewidth);
   end
   fprintf('\n');
   fprintf('  and a decimation factor of %g\n',decfactor);
   %fprintf('4 Size of allamps = [%d %d]\n',size(allamps,1),size(allamps,2));
   [allamps,outtrials] = movav(allamps,1:ntrials,avewidth,decfactor); 
                                            % Note: using square window
   %fprintf('5 Size of allamps = [%d %d]\n',size(allamps,1),size(allamps,2));
   [outsort,outtrials] = movav(sortvar,1:ntrials,avewidth,decfactor); 
   fprintf('Output data will be %d frames by %d smoothed trials.\n',...
                          frames,length(outtrials));
 else
  outtrials = 1:ntrials;
  outsort = sortvar;
 end
 allamps = 20*log10(allamps);

 if isnan(baseamp) % if not specified in 'limits'
    [amps,baseamp] = rmbase(amps,length(times),base); % remove (log) baseline
 	allamps = allamps - baseamp; % divide by (non-log) baseline amplitude
 else
    amps = amps-baseamp; % use specified (log) baseamp
	allamps = allamps - baseamp; % divide by (non-log) baseline amplitude
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

elseif exist('data2') %%%%%% Plot allcohers instead of data %%%%%%%%%%%%%%%%%%%

 if freq > 0 
    coherfreq = freq; % use phase-sort frequency
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

end %%%%%%%%%%%%%%%%%%%%%%%%%%% End image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(verttimes)
 if size(verttimes,1) ~= 1 & size(verttimes,1) ~= ntrials
    fprintf('\nerpimage(): vert arg matrix must have 1 or %d rows\n',ntrials);
    return
 end;
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
           plot([vt vt],[0 max(outtrials)],DOTSTYLE,'Linewidth',VERTWIDTH);
         elseif length(vt)==ntrials
           [outvt,ix] = movav(vt,1:ntrials,avewidth,decfactor); 
           plot(outvt,outtrials,DOTSTYLE,'Linewidth',VERTWIDTH);
         end
       else
         if length(vt)==1
           plot([0 max(outtrials)],[vt vt],DOTSTYLE,'Linewidth',VERTWIDTH);
         elseif length(vt)==ntrials
           [outvt,ix] = movav(vt,1:ntrials,avewidth,decfactor); 
           plot(outtrials,outvt,DOTSTYLE,'Linewidth',VERTWIDTH);
         end
       end
     else                % re-aligned data
       if TIMEX          % overplot vt on image
         if length(vt)==ntrials
           [outvt,ix] = movav(vt,1:ntrials,avewidth,decfactor); 
           plot(aligntime+outvt-outsort,outtrials,DOTSTYLE,'LineWidth',VERTWIDTH); 
         elseif length(vt)==1
           plot(aligntime+vt-outsort,outtrials,DOTSTYLE,'LineWidth',VERTWIDTH); 
         end
       else
         if length(vt)==ntrials
           [outvt,ix] = movav(vt,1:ntrials,avewidth,decfactor); 
           plot(outtrials,aligntime+outvt-outsort,DOTSTYLE,'LineWidth',VERTWIDTH); 
         elseif length(vt)==1
           plot(outtrials,aligntime+vt-outsort,DOTSTYLE,'LineWidth',VERTWIDTH); 
         end
       end
     end
  end
 %end
  fprintf('\n');
end

set(gca,'FontSize',TICKFONT)
hold on;
%
%%%%%%%%%%% plot vertical line at 0 or align time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~isnan(aligntime) % if trials time-aligned 
 if times(1) <= aligntime & times(frames) >= aligntime
  plot([aligntime aligntime],[0 ntrials],'k','Linewidth',ZEROWIDTH); 
     % plot vertical line at aligntime
 end
else % trials not time-aligned 
 if times(1) <= 0 & times(frames) >= 0
  plot([0 0],[0 ntrials],'k','Linewidth',ZEROWIDTH); % plot vertical line at time 0
 end
end

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

if TIMEX
 if Nosort == YES
  l=ylabel('Trial Number');
 else
  if exist('phargs')
    l=ylabel('Phase-sorted Trials');
    l=ylabel('Trials');
  elseif exist('ampargs')
    l=ylabel('Amplitude-sorted Trials');
    l=ylabel('Trials');
  else
    l=ylabel('Sorted Trials');
  end
 end
else % if switch x<->y axes
 if Nosort == YES & NoTimeflag==NO
  l=xlabel('Trial Number');
 else
  if exist('phargs')
    l=ylabel('Phase-sorted Trials');
    l=ylabel('Trials');
  elseif NoTimeflag == NO
    l=xlabel('Sorted Trials');
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
     l=xlabel('Time (ms)');
  else
     l=xlabel('Frames');
  end
  set(l,'Fontsize',LABELFONT);
end
%
%%%%%%%%%%%%%%%%%%%% Overplot sortvar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if Noshow == YES
  fprintf('Not overplotting sorted sortvar on data.\n');

elseif isnan(aligntime) % plot sortvar on un-aligned data

  if Nosort == NO;
    fprintf('Overplotting sorted sortvar on data.\n');
  end
  hold on; 
  if TIMEX      % overplot sortvar
    plot(outsort,outtrials,'k','LineWidth',SORTWIDTH); 
  else
    plot(outtrials,outsort,'k','LineWidth',SORTWIDTH);
  end                                                 
  drawnow
else % plot re-aligned zeros on sortvar-aligned data
  if Nosort == NO;
    fprintf('Overplotting sorted sortvar on data.\n');
  end
  hold on; 
  if TIMEX      % overplot aligned sortvar on image
    plot([aligntime aligntime],[0 ntrials],'k','LineWidth',SORTWIDTH);
  else
    plot([[0 ntrials],aligntime aligntime],'k','LineWidth',SORTWIDTH);
  end
  fprintf('Overplotting realigned times-zero on data.\n');
  hold on; 

  if TIMEX      % overplot realigned 0-time on image
    plot(0+aligntime-outsort,outtrials,'k','LineWidth',ZEROWIDTH); 
  else
    plot(0+outtrials,aligntime-outsort,'k','LineWidth',ZEROWIDTH); 
  end                                                 
  drawnow
end
%
%%%%%%%%%%%%%%%%%%%% Overplot auxvar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~isempty(auxvar)
 % here first smooth auxvar and apply sorting!!!
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
  if isnan(aligntime) % plot sortvar on un-aligned data
    auxcolor = auxcolors{c};
   if TIMEX      % overplot auxvar
     plot(auxvar(c,:)',auxtrials',auxcolor,'LineWidth',SORTWIDTH); 
   else
     plot(auxtrials',auxvar(c,:)',auxcolor,'LineWidth',SORTWIDTH);
   end                                                 
   drawnow
  else % plot re-aligned zeros on sortvar-aligned data
   if TIMEX      % overplot realigned 0-time on image
     plot(0+aligntime-auxvar(c,:)',auxtrials',auxcolor,'LineWidth',ZEROWIDTH); 
   else
     plot(0+auxtrials',aligntime-auxvar(c,:)',auxcolor,'LineWidth',ZEROWIDTH); 
   end                                                 
   drawnow
  end % aligntime
 end % c
end % auxvar
%
%%%%%%%%%%%%%%%%%%%%%%%% Plot colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if Colorbar == YES
   pos=get(ax1,'Position');
   axcb=axes('Position',...
       [pos(1)+pos(3)+0.02 pos(2) ...
        0.03 pos(4)]);
   cbar(axcb,0,[mindat,maxdat]); % plot colorbar to right of image
   set(axcb,'fontsize',TICKFONT);
   % drawnow
   axes(ax1); % reset current axes to the erpimage
end
%
%%%%%%%%%%%%%%%%%%%%%%% Compute ERP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if Erpflag == YES
 axes(ax1); % reset current axes to the erpimage
 xtick = get(ax1,'Xtick');     % remember x-axis tick locations
 xticklabel = get(ax1,'Xticklabel');     % remember x-axis tick locations
 set(ax1, 'xticklabel', []);
 widthxticklabel = size(xticklabel,2);
 xticklabel = cellstr(xticklabel);
 for tmpindex = 1:length(xticklabel)
    if length(xticklabel{tmpindex}) < widthxticklabel
        spaces = char( ones(1,ceil((widthxticklabel - length(xticklabel{tmpindex}))/2) )*32);
        xticklabel{tmpindex} = [spaces xticklabel{tmpindex}];
    end;
 end;
 xticklabel = strvcat(xticklabel);
 erp=nan_mean(urdata');           % compute erp average, ignoring nan's
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
 end
 if isnan(minerp)
  fac = 1;
  minerp = 0;
  while minerp == 0
    minerp = round(fac*YEXPAND*min(erp))/fac; % minimal decimal places
    fac = 10*fac;
  end
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
 plot1erp(ax2,times,erp,limit); % plot ERP
 if ~isnan(aligntime)
   line([aligntime aligntime],[limit(3:4)*1.1],'Color','k'); % x=median sort value
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
     l=xlabel('Time (ms)');
  else
     l=xlabel('Frames');
  end
   set(l,'FontSize',LABELFONT);
  else
   if exist('NoTimesPassed')~=1
     l=ylabel('Time (ms)');
   else
     l=ylabel('Frames');
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
         plot([vt vt],[limit(3:4)],DOTSTYLE,'Linewidth',VERTWIDTH);
       else
         plot([0 max(outtrials)],[limit(3:4)],DOTSTYLE,'Linewidth',VERTWIDTH);
       end
     else
       if TIMEX      % overplot realigned vt on ERP
         plot(repmat(median(aligntime+vt-outsort),1,2),[limit(3),limit(4)],...
                                  DOTSTYLE,'LineWidth',VERTWIDTH); 
       else
         plot([limit(3),limit(4)],repmat(median(aligntime+vt-outsort),1,2),...
                                  DOTSTYLE,'LineWidth',VERTWIDTH); 
       end                                                 
     end
  end
end

 ydelta = 1/10*(limit(2)-limit(1)); 
 ytextoffset = limit(1)-1.1*ydelta;
 ynumoffset = limit(1)-0.3*ydelta; 

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
else
  erp = [];
end
%
%%%%%%%%%%%%%%%%%%%%% Plot amp, coher time series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~isnan(coherfreq) 
   if freq > 0 
      coherfreq = freq; % use phase-sort frequency
   end
   %
   %%%%%% Plot amp axis below ERP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   if ~Allampsflag %%%% don't repeat computation if already done for 'allamps'

    fprintf('Computing and plotting amplitude at %g Hz.\n',coherfreq);

    if ~isnan(signifs) | Cohsigflag==NO % don't compute or plot signif. levels
     [amps,cohers] = ...
       phasecoher(urdata,size(times,2),srate,coherfreq,DEFAULT_CYCLES);
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

   fprintf('Min, max plotting amplitudes: [%g, %g] dB\n',minamp,maxamp);
   fprintf('     relative to baseamp: %g dB\n',baseamp);
   plot1erp(ax3,times,amps,[timelimits minamp(1) maxamp(1)]); % plot AMP

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
         plot([vt vt],[minamp(1) maxamp(1)],DOTSTYLE,...
                  'Linewidth',VERTWIDTH);
       else
         plot([0 max(outtrials)],[minamp(1) maxamp(1)],DOTSTYLE,...
                  'Linewidth',VERTWIDTH);
       end
      else
       if TIMEX      % overplot realigned vt on amp
         plot(repmat(median(aligntime+vt-outsort),1,2),[minamp(1),maxamp(1)],DOTSTYLE,...
                  'LineWidth',VERTWIDTH); 
       else
         plot([minamp,maxamp],repmat(median(aligntime+vt-outsort),1,2),DOTSTYLE,...
                  'LineWidth',VERTWIDTH); 
       end                                                 
      end
     end
   end

   if Cohsigflag % plot amplitude significance levels
     hold on
      plot([timelimits(1) timelimits(2)],[ampsig(1) ampsig(1)] - mean(ampsig),'r',...
                  'linewidth',SIGNIFWIDTH);
      plot([timelimits(1) timelimits(2)],[ampsig(2) ampsig(2)] - mean(ampsig),'r',...
                  'linewidth',SIGNIFWIDTH);
   end

   t=text(ynumoffset,maxamp, num2str(maxamp,3));
   set(t,'HorizontalAlignment','right','FontSize',TICKFONT);

   t=text(ynumoffset,0, num2str(0));
   set(t,'HorizontalAlignment','right','FontSize',TICKFONT);

   t=text(ynumoffset,minamp, num2str(minamp,3));
   set(t,'HorizontalAlignment','right','FontSize',TICKFONT);

   t=text(ytextoffset,(maxamp+minamp)/2,'dB','Rotation',90);
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
   coh_handle = plot1erp(ax4,times,cohers,[timelimits mincoh maxcoh]); % plot COHER
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
         plot([vt vt],[mincoh maxcoh],DOTSTYLE,'Linewidth',VERTWIDTH);
       else
         plot([0 max(outtrials)],[mincoh maxcoh],DOTSTYLE,'Linewidth',VERTWIDTH);
       end
      else
       if TIMEX      % overplot realigned vt on coher
         plot(repmat(median(aligntime+vt-outsort),1,2),[mincoh,maxcoh],DOTSTYLE,'LineWidth',VERTWIDTH); 
       else
         plot([mincoh,maxcoh],repmat(median(aligntime+vt-outsort),1,2),DOTSTYLE,'LineWidth',VERTWIDTH); 
       end                                                 
      end
     end
   end

   t=text(ynumoffset,0, num2str(0));
   set(t,'HorizontalAlignment','right','FontSize',TICKFONT);

   t=text(ynumoffset,maxcoh, num2str(maxcoh));
   set(t,'HorizontalAlignment','right','FontSize',TICKFONT);

   t=text(ytextoffset,maxcoh/2,'ITCoh','Rotation',90);
   set(t,'HorizontalAlignment','center','FontSize',LABELFONT);
    drawnow

   if Cohsigflag % plot coherence significance level
     hold on
     plot([timelimits(1) timelimits(2)],[cohsig cohsig],'r',...
           'linewidth',SIGNIFWIDTH);
   end

   set(ax4,'Box','off','color',BACKCOLOR);
   set(ax4,'Fontsize',TICKFONT);
   if NoTimeflag==NO
     if exist('NoTimesPassed')~=1
      l=xlabel('Time (ms)');
     else
      l=xlabel('Frames');
     end
     set(l,'Fontsize',LABELFONT);
   end
   axtmp = axis;
   hztxt=text(10/13*(axtmp(2)-axtmp(1))+axtmp(1), ...
        8/13*(axtmp(4)-axtmp(3))+axtmp(3), ...
        [num2str(coherfreq,4) ' Hz']);
   set(hztxt,'fontsize',TICKFONT);
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
if (~isempty(topomap)) 
    h(12)=axes('Position',...
    [gcapos(1)+0.10*gcapos(3) gcapos(2)+0.86*gcapos(4),...
        0.20*gcapos(3) 0.14*gcapos(4)]);
    % h(12) = subplot('Position',[.10 .86 .20 .14]); 
    fprintf('Plotting a topo map in upper left.\n');
	if length(topomap) == 1
		topoplot(topomap,eloc_file,'electrodes','off', ...
				 'style', 'blank', 'emarkersize1chan', 10)
	else
		topoplot(topomap,eloc_file,'electrodes','off')
	end;
    axis('square')
    axhndls = [axhndls h(12)];
end 
%   
%%%%%%%%%%%%%%% plot a spectrum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
SPECFONT = 10;
if (~isempty(lospecHz)) 
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
    [Pxx,Pxxc,F] = psd(reshape(urdata,1,size(urdata,1)*size(urdata,2)),...
                              512,srate,winlength,0,0.05);
    plot(F,10*log10(Pxx));
    goodfs = find(F>= lospecHz & F <= hispecHz);
    maxgfs = max(10*log10(Pxx(goodfs)));
    mingfs = min(10*log10(Pxx(goodfs)));
    axis('square')
    axis([lospecHz hispecHz mingfs-1 maxgfs+1]);
    set(h(13),'Box','off','color',BACKCOLOR);
    set(h(13),'Fontsize',SPECFONT);
    l=ylabel('dB');
    set(l,'Fontsize',SPECFONT);
    if ~isnan(coherfreq)
      hold on; plot([coherfreq,coherfreq],[mingfs maxgfs],'r');
   end
   axhndls = [axhndls h(13)];
end 

%   
%%%%%%%%%%%%%%% save plotting limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
limits = [timelimits(1:2) minerp maxerp minamp maxamp mincoh maxcoh];
limits = [limits baseamp coherfreq];  % add coherfreq to output limits array
%   
%%%%%%%%%%%%%%% turn on axcopy() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
axcopy; % turn on popup zoom windows

%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End erpimage() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
axes('position',gcapos);
axis off
return
%
%%%%%%%%%%%%%%%%%%% function plot1erp() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [plot_handle] = plot1erp(ax,Time,erp,axlimits)

  ERPDATAWIDTH = 2;
  ERPZEROWIDTH = 2;
  [plot_handle] = plot(Time,erp,'LineWidth',ERPDATAWIDTH); hold on
  if sum(isnan(axlimits))==0
    if axlimits(2)>axlimits(1) & axlimits(4)>axlimits(3)
      axis([axlimits(1:2) 1.1*axlimits(3:4)])
    end
    l1=line([axlimits(1:2)],[0 0],    'Color','k',...
                 'linewidth',ERPZEROWIDTH); % y=zero-line
    l2=line([0 0],[axlimits(3:4)*1.1],'Color','k',...
                 'linewidth',ERPZEROWIDTH); % x=zero-line
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
win = exp(2i*pi*freq(:)*[1:length(nwin)]/srate);
win = win .* repmat(hanning(length(nwin))',length(freq),1);
resp = win * data(nwin,:);
ang = angle(resp);
amp = abs(resp);
if ~exist('allamps')
   allamps = [];
end
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
% nan_mean() - take the column means of a matrix, ignoring NaN values
% 
function out = nan_mean(in)

   nans = find(isnan(in));
   in(nans) = 0;
   sums = sum(in);
   nonnans = ones(size(in));
   nonnans(nans) = 0;
   nonnans = sum(nonnans);
   nononnans = find(nonnans==0);
   nonnans(nononnans) = 1;
   out = sum(in)./nonnans;
   out(nononnans) = NaN;
