% envtopo() - Plot the envelope of a data epoch plus envelopes and scalp maps of specified 
%             or largest contributing components. If a 3-D input matrix, operating on the
%             epoch mean. Click on individual topoplots to examine them in detail (by axcopy()).
% Usage:
%            >> envtopo(data,weights);
%            >> [compvarorder,compvars,compframes,comptimes,compsplotted,pvaf] ...
%                    = envtopo(data, weights, 'key1', val1, ...);
% Inputs:
%  data         = single data epoch (chans,frames) or a 3-D data matrix 
%                 (chans,frames,epochs). If data are 3-D, process the mean data epoch.
%  weights      = ICA weight matrix (e.g., icaweights*icasphere)
%
% Optional inputs:
%  'chanlocs'  = [string] channel location file or EEG.chanlocs structure. 
%                  See >> topoplot example 
%  'compnums'  = [integer array] vector of component numbers to plot {default|0 -> all}
%                  Else if int < 0, the number of largest contributing components to plot 
%                  {default|[] -> 7}
%  'timerange' = start and end input data latencies (in ms) {default: from 'limits' if any}
%  'limits'    = 0 or [minms maxms] or [minms maxms minuV maxuV]. Specify start/end plotting
%                  limits (in ms) and min/max data limits (in uV). If 0, or if both
%                  minmx & maxms == 0 -> use latencies from 'timerange' (else 0:frames-1).
%                  If both minuV and maxuV == 0 -> use data uV limits {default: 0}
%  'limcontrib' = [minms maxms]  time range (in ms) in which to rank component contribution
%                  (boundaries shown with thin dotted lines) {default|[]|[0 0] -> data limits}
%  'sortvar'    = ['mv'|'pv'|'rv'] if 'mv', sort components by back-projected variance; if 'pv', 
%                  sort by percent variance accounted for (pvaf). If 'rv', sort by relative 
%                  variance, where:                                      {default: 'mv'}
%                     pvaf(component) = 100-100*variance(data-component))/variance(data)
%                     rv(component)   = 100*variance(component)/variance(data) 
%  'title'      = [string] plot title {default|[] -> none}
%  'plotchans'  = [integer array] data channels to use in computing contributions and envelopes,
%                  and also for making scalp topo plots {default|[] -> all}
%  'voffsets'   = [float array] vertical line extentions above the data max to disentangle
%                  plot lines (left->right heads, values in y-axis units) {def|[] -> none}
%  'colors'     = [string] filename of file containing colors for envelopes, 3 chars
%                  per line, (. = blank). First color should be "w.." (white)
%                  Else, 'bold' -> plot default colors in thick lines.
%                  {default|[] -> standard Matlab color order}
%  'fillcomp'   = int_vector>0 -> fill the numbered component envelope(s) with 
%                  solid color. Ex: [1] or [1 5] {default|[]|0 -> no fill}
%  'vert'       = vector of times (in ms) at which to plot vertical dashed lines {default|[] -> none}
%  'icawinv'    = [float array] inverse weight matrix. By default computed by inverting
%                  the weight matrix (but if some components have been removed, then
%                  weight's pseudo-inverse matrix does not represent component's maps).
%  'icaact'     = [float array] ICA component activations. By default these are computed 
%                  from the input weight matrix.
%  'envmode'    = ['avg'|'rms'] compute the average envelope or the root mean square
%                  envelope {default: 'avg'}
%  'subcomps'   = [integer vector] indices of components to remove from data before 
%                  plotting. {default: none}
%  'sumenv'     = ['on'|'off'|'fill'] 'fill' -> show the filled envelope of the summed projections 
%                  of the selected components; 'on' -> show the envelope only {default: 'fill'}
%  'actscale'   = ['on'|'off'] scale component scalp maps by maximum component activity in the
%                  designated (limcontrib) interval. 'off' -> scale scalp maps individually using
%                  +/-max(abs(map value)) {default: 'off'}
%  'dispmaps'   = ['on'|'off'] display component number and scalp maps. {default: 'on'}
% Outputs:
%  compvarorder = component numbers in decreasing order of max variance in data
%  compvars     = component max variances
%  compframes   = frames of max variance
%  comptimes    = times of max variance
%  compsplotted = components plotted
%  mv|pvaf|rv   = max variance, percent variance accounted for, or relative variance (see 'sortvar')
%
% Notes:
%  To label maps with other than component numbers, put four-char strings into
%  a local pwd file named 'envtopo.labels' (. = space) in time-order of their projection maxima
%
% Authors: Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 3/1998 
%
% See also: timtopo()

% Copyright (C) 3-10-98 from timtopo.m Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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
% Revision 1.83  2004/11/29 23:13:45  scott
% compute variances accounted for ONLY from g.plotchans, if any
%
% Revision 1.82  2004/11/20 22:48:15  scott
% fixed error detection using old-style argument list (no longer documented)
%
% Revision 1.81  2004/11/18 19:11:11  scott
% adjust chanlocs = g.chanlocs
%
% Revision 1.80  2004/11/18 03:11:55  scott
% chanlocs -> subsample channels for topoplots
%
% Revision 1.79  2004/11/17 18:06:33  scott
% debug 'vert' (specify latencies in ms, documented), limcontrib and vert line styles, 'pvaf'
%
% Revision 1.78  2004/11/17 02:27:10  scott
% made xmin xmax in 'limits' arg control the times plotted
% changed key 'pvaf' to 'sortvar' ('pvaf' retained for backwards compatiability)
% documented the mean(data,3) option - was undocumented
%
% Revision 1.77  2004/11/13 19:29:10  scott
% print which components subtracted
%
% Revision 1.76  2004/11/13 19:27:12  scott
% added test for 'subcomps' argument
%
% Revision 1.75  2004/11/13 06:47:25  scott
% debugged yaxis limits selection
%
% Revision 1.74  2004/11/12 18:38:35  scott
% help msg details
%
% Revision 1.73  2004/11/12 18:35:34  scott
% reformed the limit setting code. Adjusted component label text, [ymin,ymax].
% Added xunitframes and ylimset variables.
%
% Revision 1.72  2004/11/12 07:21:02  scott
% debug imerange
%
% Revision 1.71  2004/11/12 07:14:53  scott
% debugging 'limcontrib'; changed 'times' keyword to 'timerange'
%
% Revision 1.70  2004/11/11 15:12:33  scott
% clean up
%
% Revision 1.69  2004/11/11 14:41:48  scott
% made new 'times' option [minms, maxms]
%
% Revision 1.68  2004/11/11 14:25:24  scott
% Added 'times' input (for pop_envtopo() use); made default component selection
% method 'mv' (max variance). Changed 'pvaf' options to 'mv' 'rv' 'pv' (but left
% the backwards compatible meanings of 'on' 'off' and 'pvaf'. Tuned the help message
% and printed messages.
%
% Revision 1.67  2004/11/11 07:20:57  scott
% rescale [ymin,ymax] to include largest components considered -sm
%
% Revision 1.66  2004/11/11 06:58:53  scott
% debugging limits -sm
%
% Revision 1.65  2004/11/05 03:39:53  scott
% nothing
%
% Revision 1.64  2004/08/03 22:01:47  arno
% time in seconds
%
% Revision 1.63  2004/07/30 01:26:14  arno
% changed the figure to be the current figure
%
% Revision 1.62  2004/07/30 01:03:52  arno
% changed text to double
%
% Revision 1.61  2004/05/06 23:41:25  scott
% typo
%
% Revision 1.60  2004/05/04 05:35:25  scott
% finished debugging limits
%
% Revision 1.59  2004/05/04 05:24:10  scott
% debug setting g.limits(2) to ms
%
% Revision 1.58  2004/04/28 06:00:53  scott
% debug 'compnums',[2 3 4 6] etc
%
% Revision 1.57  2004/04/25 16:44:41  scott
% converted limits and limcontrib inputs back to ms (for compatibility with pop_envtopo)
%
% Revision 1.56  2004/04/25 16:20:41  scott
% clarified that input times are in sec, not ms; cleaned up commandline info -sm
%
% Revision 1.55  2004/04/24 17:21:37  scott
% fixed pvaf computation and printout. Added 'sumenv' mode as default. Deprecated
% 'colorfile' for 'colors' (backward compatible). Cleaned up printout.
% Edited help message.
%
% Revision 1.54  2004/03/03 21:44:00  arno
% rv for residual variance
%
% Revision 1.53  2004/03/03 21:26:02  arno
% debugging sorting of components; relative variance; text output
%
% Revision 1.52  2004/03/03 18:54:29  arno
% retreive version 1.50
%
% Revision 1.50  2004/02/03 15:58:00  arno
% no interpreter for title
%
% Revision 1.49  2004/01/29 16:56:30  scott
% same
%
% Revision 1.45  2004/01/29 16:50:45  scott
% printout edit
%
% Revision 1.44  2004/01/29 16:44:56  scott
% rm pvaf printout for now
%
% Revision 1.43  2004/01/29 16:39:45  scott
% test for chanlocs file location and size
%
% Revision 1.42  2004/01/26 02:22:13  scott
% same
%
% Revision 1.28  2004/01/26 00:45:14  scott
% improved listing of pvaf in Matlab command window
%
% Revision 1.27  2003/12/03 18:31:35  scott
% percentage -> percent
%
% Revision 1.26  2003/09/17 02:00:06  arno
% debuging pvaf assignment when using compnums
%
% Revision 1.25  2003/07/30 01:56:06  arno
% adding 'actscale' option and cbar
%
% Revision 1.24  2003/04/15 16:55:10  arno
% allowing to plot up to 20 components
%
% Revision 1.23  2003/03/23 20:14:50  scott
% fill msg edit
%
% Revision 1.22  2003/03/23 20:07:38  scott
% making data env overplot bold
%
% Revision 1.21  2003/03/23 20:05:59  scott
% overplot data envelope on filled component projection
%
% Revision 1.20  2003/03/14 16:18:37  arno
% plot pvaf in topoplot
%
% Revision 1.19  2003/03/14 15:23:29  arno
% pvaf -> * 100
%
% Revision 1.18  2003/03/14 02:50:18  arno
% help for pvaf
%
% Revision 1.17  2003/03/14 02:45:48  arno
% now printing pvaf
%
% Revision 1.16  2003/03/11 01:45:03  arno
% first capital letter in labels
%
% Revision 1.15  2003/03/10 18:24:47  arno
% ploting only one contribution
%
% Revision 1.14  2003/03/08 21:02:57  arno
% debugging
%
% Revision 1.13  2003/03/08 00:11:06  arno
% allowing input of ICA component activity
%
% Revision 1.12  2003/03/05 03:23:44  scott
% minor
%
% Revision 1.11  2003/03/03 22:28:03  arno
% update text header
%
% Revision 1.10  2003/02/27 02:52:24  arno
% typo
%
% Revision 1.9  2002/10/25 18:47:33  luca
% sign of comparison typo - ad
%
% Revision 1.8  2002/10/09 21:57:25  arno
% documenting limcontrib
%
% Revision 1.7  2002/10/09 21:32:18  arno
% documenting option subcomp
%
% Revision 1.6  2002/10/05 01:52:20  arno
% debug envmode
%
% Revision 1.5  2002/10/05 01:50:34  arno
% new function with 'key', 'val' args, extra params: envmode, limcontrib, icawinv...
%
% Revision 1.4  2002/09/05 00:57:22  arno
% colorbar->cbar for removing menu bug
%
% Revision 1.3  2002/04/25 17:22:33  scott
% editted help msg -sm
%
% Revision 1.2  2002/04/09 02:13:22  arno
% make the color file internal
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% Edit History:
% 3-18-98 fixed bug in LineStyle for fifth component, topoplot maxproj with 
%         correct orientations, give specified component number labels -sm
% 4-28-98 plot largest components, ranked by max projected variance -sm
% 4-30-98 fixed bug found in icademo() -sm
% 5-08-98 fixed bug found by mw () -sm
% 5-23-98 made vert. line styles for comps 6 & 11 correct -sm
% 5-30-98 added 'envtopo.labels' option -sm
% 5-31-98 implemented plotchans arg -sm
% 7-13-98 gcf->gca to allow plotting in subplots -sm
% 9-18-98 worked more to get plotting in subplot to work -- no luck yet! -sm
% 2-22-99 draw oblique line to max env value if data clipped at max proj -sm
% 2-22-99 added colorfile -sm
% 4-17-99 added support for drawing in subplots -t-pj
% 10-29-99 new versions restores search through all components for max 7 and adds 
%          return variables (>7 if specified. Max of 7 comp envs still plotted. -sm
% 11-17-99 debugged new version -sm
% 12-27-99 improved help msg, moved new version to distribution -sm
% 01-21-00 added 'bold' option for colorfile arg -sm
% 02-28-00 added fill_comp_env arg -sm
% 03-16-00 added axcopy() -sm & tpj
% 05-02-00 added vert option -sm
% 05-30-00 added option to show "envelope" of only 1 channel -sm
% 09-07-00 added [-n] option for compnums, added BOLD_COLORS as default -sm
% 12-19-00 updated icaproj() args -sm
% 12-22-00 trying 'axis square' for topoplots -sm
% 02-02-01 fixed bug in printing component 6 env line styles -sm
% 04-11-01 added [] default option for args -sm
% 01-25-02 reformated help & license, added links -ad 
% 03-15-02 added readlocs and the use of eloc input structure -ad 
% 03-16-02 added all topoplot options -ad

function [compvarorder,compvars,compframes,comptimes,compsplotted,pvaf] = envtopo(data,weights,varargin)

% icadefs;    % read toolbox defaults

pvaf = []; % note previous arguments 'on' -> 'pv', 'off' -> 'rv'
           % for backwards compatibility 11/2004 -sm
all_bold = 0;
BOLD_COLORS = 1;    % 1 = use solid lines for first 5 components plotted
                    % 0 = use std lines according to component rank only
FILL_COMP_ENV = 0;  % default no fill
FILLCOLOR   = [.815 .94 1]; % use lighter blue for better env visibility
% FILLCOLOR = [.66 .76 1];
MAXTOPOS = 20;      % max topoplots to plot
VERTWEIGHT = 2.0;  % lineweight of specified vertical lines
LIMCONTRIBWEIGHT = 1.2; % lineweight of limonctrib vertical lines

myfig =gcf;         % remember the current figure (for Matlab 7.0.0 bug)
xmin = 0; xmax = 0;
    
if nargin < 2
   help envtopo
   return
end

if nargin <= 2 | isstr(varargin{1})
	% 'key' 'val' sequences
	fieldlist = { 'chanlocs'      ''         []                       [] ;
				  'title'         'string'   []                       '';
				  'limits'        'real'     []                       0;
				  'timerange'     'real'     []                       [];
				  'plotchans'     'integer'  [1:size(data,1)]         [] ;
				  'icawinv'       'real'     []                       pinv(weights) ;
				  'icaact'        'real'     []                       [] ;
				  'voffsets'      'real'     []                       [] ;
				  'vert'          'real'     []                       [] ;
				  'fillcomp'      'integer'  []                       0 ; 
				  'colorfile'     'string'   []                       '' ; 
				  'colors'        'string'   []                       '' ;
				  'compnums'      'integer'  []                       -7; 
				  'subcomps'      'integer'  []                       []; 
				  'envmode'       'string'   {'avg' 'rms'}            'avg'; 
				  'dispmaps'      'string'   {'on' 'off'}             'on'; 
				  'pvaf'          'string'   {'mv' 'rv' 'pv' 'pvaf' 'on' 'off' ''}  ''; 
				  'sortvar'       'string'   {'mv' 'rv' 'pv' 'pvaf' 'on' 'off'} 'mv';  
				  'actscale'      'string'   {'on' 'off'}             'off'; 
				  'limcontrib'    'real'     []                       0;  
                                  'sumenv'        'string'    {'on' 'off' 'fill'}     'fill'};
	
	[g varargin] = finputcheck( varargin, fieldlist, 'envtopo', 'ignore');
	if isstr(g), error(g); end;

else % dprecated - old style input args
	if nargin > 3,    g.chanlocs = varargin{1};
	else              g.chanlocs = [];
	end;
        if isstr(varargin{2}), help envtopo; return; end
	if nargin > 4,	  g.limits = varargin{2};
	else              g.limits = 0; % [];
	end;
        if isstr(varargin{3}), help envtopo; return; end
	if nargin > 5,    g.compnums = varargin{3};
	else              g.compnums = [];
	end;
        if ~isstr(varargin{4}), help envtopo; return; end
	if nargin > 6,    g.title = varargin{4};
	else              g.title = '';
	end;
        if isstr(varargin{5}), help envtopo; return; end
	if nargin > 7,    g.plotchans = varargin{5};
	else              g.plotchans = [];
	end;
        if isstr(varargin{6}), help envtopo; return; end
	if nargin > 8,    g.voffsets = varargin{6};
	else              g.voffsets = [];
	end;
        if isstr(varargin{7}), help envtopo; return; end
	if nargin > 9,    g.colorfile = varargin{7};
	else              g.colorfile = '';
	                  g.colors = '';
	end;
        if isstr(varargin{8}), help envtopo; return; end
	if nargin > 10,   g.fillcomp = varargin{8};
	else              g.fillcom = 0;
	end;
        if isstr(varargin{9}), help envtopo; return; end
	if nargin > 11,   g.vert = varargin{9};
	else              g.vert = [];
	end;
    if nargin > 12, varargin =varargin(10:end); end;
    g.sumenv = 'on';
    g.sortvar = 'mv';
    g.timerange = [];
    g.pvaf = ''; 
    g.pvaf = ''; 
    g.icaact = [];
    g.limcontrib = 0;
    g.icawinv = pinv(weights);
    g.subcomps = [];
    g.envmode = 'avg';
    g.dispmaps = 'on';
end;

if ~isempty(g.colors)
    g.colorfile = g.colors; % retain old usage 'colorfile' for 'colors' -sm 4/04
end

if ndims(data) == 3
    data = mean(data,3); % average the data if 3-D
end;
[chans,frames] = size(data);

if ~isempty(g.vert)
   g.vert = g.vert/1000; % convert from ms to s
end
%
%%%%%% Collect information about the gca, then delete it %%%%%%%%%%%%%
%
uraxes = gca; % the original figure or subplot axes
pos=get(uraxes,'Position');
axcolor = get(uraxes,'Color');
delete(gca)

%
%%% Convert g.timerange, g.limits and g.limcontrib to sec from ms %%%%
%
g.timerange = g.timerange/1000;   % the time range of the input data
g.limits(1) = g.limits(1)/1000;   % the time range to plot
if length(g.limits) == 1   % make g.limits at least of length 2
    g.limits(1) = 0; g.limits(2) = 0;
else
    g.limits(2) = g.limits(2)/1000;  % 
end;
g.limcontrib = g.limcontrib/1000; % the time range in which to select largest components

%
%%%%%%%%%%%% Collect time range information %%%%%%%%%%%%%%%%%%%%%%%%%%
%
if length(g.limits) == 3 | length(g.limits) > 4 % if g.limits wrong length
   fprintf('envtopo: limits should be 0, [minms maxms], or [minms maxms minuV maxuV].\n');
end


xunitframes = 0; % flag plotting if xmin & xmax are in frames instead of sec
if ~isempty(g.timerange)   % if 'timerange' given
    if g.limits(1)==0 & g.limits(2)==0
         g.limits(1) = min(g.timerange); % if no time 'limits
         g.limits(2) = max(g.timerange); % plot whole data epoch
    end
else % if no 'timerange' given
    if g.limits(1)==0 & g.limits(2)==0 % if no time limits as well, 
         fprintf('\nNOTE: No time limits given: using 0 to %d frames\n',frames-1);
         g.limits(1) = 0;
         g.limits(2) = frames-1;
         xunitframes     = 1; % mark frames as time unit instead of sec
    end
end
 
if isempty(g.timerange)
  xmin = g.limits(1); % (xmin, xmax) are data limits in sec
  xmax = g.limits(2);
else 
  xmin = g.timerange(1); % (xmin, xmax) are data limits in sec
  xmax = g.timerange(2);
end
pmin = g.limits(1); % plot min and max sec
pmax = g.limits(2);

dt = (xmax-xmin)/(frames-1);  % sampling interval in sec
times=xmin*ones(1,frames)+dt*(0:frames-1); % time points in sec

%
%%%%%%%%%%%% Collect y-axis range information %%%%%%%%%%%%%%%%%%%%%%%%
%
ylimset = 0; % flag whether hard limits have been set by the user
ymin = min(min(data)); % begin by setting limits from data
ymax = max(max(data));
if length(g.limits) == 4 
     if g.limits(3)~=0 | g.limits(4)~=0 % collect plotting limits from 'limits'
	 ymin = g.limits(3);
	 ymax = g.limits(4);
         ylimset = 1;
     end
end

%
%%%%%%%%%%%%%%% Find limits of the component selection window %%%%%%%%%
% 
if any(g.limcontrib ~= 0) 
    if xunitframes
       g.limcontrib = g.limcontrib*1000; % if no time limits, interpret
    end                                  % limcontrib as frames
    if g.limcontrib(1)<xmin
          g.limcontrib(1) = xmin;
    end
    if g.limcontrib(2)>xmax
          g.limcontrib(2) = xmax;
    end
    srate = (frames-1)/(xmax-xmin);
    limframe1  = round((g.limcontrib(1)-xmin)*srate)+1;
    limframe2  = round((g.limcontrib(2)-xmin)*srate)+1;
    g.vert(end+1) =  g.limcontrib(1);
    g.vert(end+1) =  g.limcontrib(2);
else
    limframe1 = 1;
    limframe2 = frames;
end;

%
%%%%%%%%%%%%%%%%%%%%% Read line color information %%%%%%%%%%%%%%%%%%%%%
%
ENVCOLORS = strvcat('w..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..');

if isempty(g.colorfile)
    g.colorfile = ENVCOLORS; % use default color order above
elseif ~isstr(g.colorfile)
      error('Color file name must be a string.');
end
if strcmpi(g.colorfile,'bold')
       all_bold = 1;
       g.colorfile = ENVCOLORS; % default colors 
end
if exist(g.colorfile) == 2  % if an existing file
        cid = fopen(g.colorfile,'r');
        if cid <3,
            error('cannot open color file');
        else
            colors = fscanf(cid,'%s',[3 MAXENVPLOTCHANS]);
            colors = colors';
        end;
else
        colors = g.colorfile;
end
[r c] = size(colors);
      for i=1:r
        for j=1:c
          if colors(i,j)=='.',
            if j==1
              error('Color file should have color letter in 1st column.');
            elseif j==2
              colors(i,j)='-';
            elseif j>2
              colors(i,j)=' ';
            end;
          end;
        end;
      end;
colors(1,1) = 'k'; % make sure 1st color (for data envelope) is black

%
%%%%%%%%%%%%%%%% Check other input variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[wtcomps,wchans] = size(weights);
if wchans ~= chans
     error('Sizes of weights and data do not agree');
end

if isempty(g.voffsets) | ( size(g.voffsets) == [1,1] & g.voffsets(1) == 0 )
    g.voffsets = zeros(1,MAXTOPOS); 
end
if isempty(g.plotchans) | g.plotchans(1)==0 
    g.plotchans = 1:chans;
end
if max(g.plotchans) > chans | min(g.plotchans) < 1
    error('invalid ''plotchan'' index');
end
if isempty(g.compnums) | g.compnums(1) == 0
    g.compnums = 1:wtcomps; % by default, all components
end
if min(g.compnums) < 0
    if length(g.compnums) > 1
	     error('Negative compnums must be a single integer.');
    end
    if -g.compnums > MAXTOPOS
	    fprintf('Can only plot a maximum of %d components.\n',MAXTOPOS);
	    return
    else
	    MAXTOPOS = -g.compnums;
	    g.compnums = 1:wtcomps;
    end
end

g.subcomps = abs(g.subcomps); % don't pass negative channel numbers
if min(g.subcomps)<1 | max(g.subcomps) > chans
    error('Keyword ''subcomps'' argument out of bounds');
end

%
%%%%%%%%%%%%%%% Subtract components from data if requested %%%%%%%%%%%%%
%
if ~isempty(g.subcomps)
	    fprintf('Subtracting requested components from plotting data: ');
            for c=g.subcomps
              fprintf('%d ',c);
            end
            fprintf('\n');
	    proj = icaproj(data,weights,g.subcomps); % updated arg list 12/00 -sm
	    data = data - proj;
end;
	    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Process components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
ncomps = length(g.compnums);
for i=1:ncomps-1
    if g.compnums(i) == 0
	       fprintf('Removing component number 0 in compnums.\n');
	       g.compnums(i)=[];
    elseif g.compnums(i)>wtcomps
	       fprintf('compnums(%d) > number of comps (%d)?\n',i,wtcomps);
	       return
    end
    for j=i+1:ncomps
	    if g.compnums(i)==g.compnums(j)
	       fprintf('Removing repeated component number (%d) in compnums.\n',g.compnums(i));
	       g.compnums(j)=[];
            end
	  end
    end

    limitset = 0;
    if isempty(g.limits)
      g.limits = 0;
    end
    if length(g.limits)>1
      limitset = 1;
    end

    %
    %%%%%%%%%%%%%%% Compute plotframes and envdata %%%%%%%%%%%%%%%%%%%%%
    %
    ntopos = length(g.compnums);
    if ntopos > MAXTOPOS
      ntopos = MAXTOPOS; % limit the number of topoplots to display
    end

    if max(g.compnums) > wtcomps | min(g.compnums)< 1
      fprintf(...
       'envtopo(): one or more compnums out of range (1,%d).\n',wtcomps);
      return
    end

    plotframes = ones(ncomps);
    maxproj = zeros(chans,ncomps);
    %
    % first, plot the data envelope
    %
    envdata = zeros(2,frames*(ncomps+1));
    envdata(:,1:frames) = envelope(data(g.plotchans,:), g.envmode); 

    fprintf('Data epoch is from %.0f ms to %.0f ms.\n',1000*xmin,1000*xmax);
    fprintf('Plotting data from %.0f ms to %.0f ms.\n',1000*xmin,1000*xmax);
    fprintf('Comparing projection sizes for components:  ');
        if ncomps>32
           fprintf('\n');
        end
    compvars = zeros(1,ncomps);

    %
    %%%%%%%%%%%%%% find max variances and their frame indices %%%%%%%%%%%
    %
    for c = 1:ncomps 
      if ~rem(c,5) 
          fprintf('%d ... ',g.compnums(c)); % c is index into compnums
      end
      if ~rem(c,100)
        fprintf('\n');
      end

      if isempty(g.icaact)
          proj = g.icawinv(g.plotchans,g.compnums(c))*weights(g.compnums(c),:)*data; % updated -sm 11/04
      else 
          proj = g.icawinv(g.plotchans,g.compnums(c))*g.icaact(g.compnums(c),:);     % updated -sm 11/04
      end;                                                % now proj has only g.plotchans and g.compnums     

      envdata(:,c*frames+1:(c+1)*frames) = envelope(proj(:,:), g.envmode);

      [maxval,maxi] = max(sum(proj(:,limframe1:limframe2).*proj(:,limframe1:limframe2))); 
                                  % find max variance
      compvars(c)   = maxval;
      %
      %%%%%% find variance in interval after removing component %%%%%%%%%%%
      %
      if isempty(g.pvaf)
                 g.pvaf = g.sortvar;
      end
      if strcmpi(g.pvaf, 'pvaf') | strcmpi('pvaf','on') | strcmpi(g.pvaf,'mv')
              pvaf(c) = mean(mean((data(g.plotchans,limframe1:limframe2)-proj(:,limframe1:limframe2)).^2)); 
      else % if 'pvaf' is 'rv' (or, formerly, 'off')
              pvaf(c) = mean(mean(proj(:,limframe1:limframe2).^2));      
      end;

      maxi = maxi+limframe1-1;
      plotframes(c) = maxi;
      maxproj(:,c)  = proj(:,maxi); % Note: maxproj contains only g.plotchans -sm 11/04
      
      %    ix = find(envdata(1,c*frames+1:(c+1)*frames) > ymax);
      %    [val,ix] = max(envdata(1,c*frames+ix));
      %    plotframes(c) = ix; % draw line from max non-clipped env maximum
      %    maxproj(:,c)  = proj(:,ix);
      % else draw line from max envelope value at max projection time point

end % component c
fprintf('\n');

%
%%%%%%%%%%%%%%% Compute component selection criterion %%%%%%%%%%%%%%%%%%%%%%%%%%
%

% compute pvaf
if ~xunitframes
  fprintf('  from the interval %3.0f ms to %3.0f ms.\n',1000*times(limframe1),1000*times(limframe2));
end
vardat = mean(mean((data(g.plotchans,limframe1:limframe2).^2))); % find data variance in interval

if strcmpi(g.pvaf, 'pvaf') | strcmpi(g.pvaf,'on') | strcmpi(g.pvaf,'mv')
    pvaf = 100-100*pvaf/vardat;
    ot   = g.pvaf;
    if strcmpi(ot,'on'), ot = 'pvaf'; end
elseif strcmpi(g.pvaf, 'rv')| strcmpi(g.pvaf,'off')
    pvaf = 100*pvaf/vardat;
    ot   = 'rv';
end;
[sortpvaf spx] = sort(pvaf);
sortpvaf = sortpvaf(end:-1:1);
spx      = spx(end:-1:1);
npercol = ceil(ncomps/3);

% for index =1:npercol
%     try, fprintf('   IC%d\t%s: %6.2f%%\t', spx(index)          , ot, sortpvaf(index)); catch, end;
%     try, fprintf('   IC%d\t%s: %6.2f%%\t', spx(index+npercol)  , ot, sortpvaf(index+npercol)); catch, end;
%     try, fprintf('   IC%d\t%s: %6.2f%%\t', spx(index+2*npercol), ot, sortpvaf(index+2*npercol)); catch, end;
%     fprintf('\n');
% end;

%
%%%%%%%%%%%%%%%%%%%%%%%%% Sort by max variance in data %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
sampint  = (xmax-xmin)/(frames-1);     % sampling interval in sec
times    = xmin:sampint:xmax;          % make vector of data time values

[v minf] = min(abs(times-pmin));
[v maxf] = max(abs(times-pmax));
pframes  = minf:maxf;         % frames to plot
ptimes   = times(pframes);    % times to plot

if strcmpi(g.pvaf,'mv') | strcmpi(g.pvaf,'on')
    [compvars,compx] = sort(compvars');  % sort compnums on max variance
else % if 'pvaf' is 'pvaf' or 'rv'
    [tmp,compx]  = sort(pvaf');   % sort compnums on pvaf or rv
end

compx        = compx(ncomps:-1:1);    % reverse order of sort
compvars     = compvars(ncomps:-1:1)';% reverse order of sort (output var)
compvarorder = g.compnums(compx);     % actual component numbers (output var)
plotframes   = plotframes(compx);     % plotted comps have these max frames 
compframes   = plotframes';           % frame of max variance in each comp (output var)
comptimes    = times(plotframes(compx));  % time of max variance in each comp (output var)
compsplotted = compvarorder(1:ntopos);% (output var)
%
%%%%%%%%%%%%%%%%%%%%%%%% Reduce to ntopos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[plotframes,ifx] = sort(plotframes(1:ntopos));% sort plotframes on their temporal order
plottimes  = times(plotframes);       % convert to times in ms
compx      = compx(ifx);              % indices into compnums, in plotting order
maporder   = g.compnums(compx);       % comp. numbers, in plotting order (l->r)
maxproj    = maxproj(:,compx);        % maps in plotting order 

vlen = length(g.voffsets); % extend voffsets if necessary
while vlen< ntopos
  g.voffsets = [g.voffsets g.voffsets(vlen)]; % repeat last offset given
  vlen=vlen+1;
end

head_sep = 1.2;
topowidth = pos(3)/(ntopos+(ntopos-1)/5); % width of each topoplot
if topowidth > 0.20    % adjust for maximum height
    topowidth = 0.2;
end

if rem(ntopos,2) == 1  % odd number of topos
   topoleft = pos(3)/2 - (floor(ntopos/2)*head_sep + 0.5)*topowidth;
else % even number of topos
   topoleft = pos(3)/2 - (floor(ntopos/2)*head_sep)*topowidth;
end

%
%%%%%%%%%%%%%%%%%%%% Print times and frames of comp maxes %%%%%%%%%%%%%%
%

% fprintf('\n');
fprintf('Plotting envelopes of %d component projections.\n',ntopos);
if length(g.plotchans) ~= chans
  fprintf('Envelopes computed from %d specified data channels.\n',...
      length(g.plotchans));
end
fprintf('Topo maps will show components: ');
for t=1:ntopos
  fprintf('%4d  ',maporder(t));
end
fprintf('\n');
if ~xunitframes
  fprintf('    with max var at times (ms): ');
  for t=1:ntopos
    fprintf('%4.0f  ',1000*plottimes(t));
  end
  fprintf('\n');
end

fprintf('                  epoch frames: ');
for t=1:ntopos
   fprintf('%4d  ',limframe1-1+plotframes(t));
end
fprintf('\n');

if strcmpi(g.pvaf,'on') | strcmpi(g.pvaf,'pvaf')
   fprintf('    Component pvaf in interval:  ');
   for t=1:ntopos
     fprintf('%4.2f ',pvaf(t));
   end
   fprintf('\n');
end

sumproj = zeros(size(data(g.plotchans,:)));
for n = 1:ntopos
  if isempty(g.icaact)
      sumproj = sumproj + g.icawinv(g.plotchans,maporder(n))*weights(maporder(n),:)*data; % updated -sm 11/04
  else 
      sumproj = sumproj + g.icawinv(g.plotchans,maporder(n))*g.icaact(maporder(n),:);     % updated -sm 11/04
  end;                                                              % Note: sumproj here has only g.plotchans
end
varproj = mean(mean((data(g.plotchans,limframe1:limframe2).^2))); % find data variance in interval
if strcmpi(g.pvaf, 'on')
      sumpvaf = mean(mean((data(g.plotchans,limframe1:limframe2)-sumproj(:,limframe1:limframe2)).^2)); 
  else
      sumpvaf = mean(mean(sumproj(:,limframe1:limframe2).^2));      
  end;
if strcmpi(g.pvaf, 'on') | strcmpi(g.pvaf,'pv') | strcmpi(g.pvaf,'pvaf') | strcmpi(g.pvaf,'mv')
    sumpvaf = 100-100*sumpvaf/varproj;
    ot   = 'pvaf';
else % if strcmpi(g.pvaf, 'off') | strcmpi(g.pvaf,'rv')
    sumpvaf = 100*sumpvaf/varproj;
    ot   = 'rv';
end;
if ~xunitframes
   fprintf('    Summed component %s in interval [%4g %4g] ms: %4.2f%%\n',ot, 1000*times(limframe1),1000*times(limframe2), sumpvaf);
end
%
%%%%%%%%%%%%%%%%%%%%% Plot the data envelopes %%%%%%%%%%%%%%%%%%%%%%%%%
%
BACKCOLOR = [0.7 0.7 0.7];
newaxes=axes('position',pos);
axis off
set(newaxes,'FontSize',16,'FontWeight','Bold','Visible','off');
set(newaxes,'Color',BACKCOLOR); % set the background color
delete(newaxes) %XXX

% site the plot at bottom of the current axes
axe = axes('Position',[pos(1) pos(2) pos(3) 0.6*pos(4)],...
           'FontSize',16,'FontWeight','Bold');

g.limits = get(axe,'Ylim');
set(axe,'GridLineStyle',':')
set(axe,'Xgrid','off')
set(axe,'Ygrid','on')
axes(axe)
set(axe,'Color',axcolor);

% fprintf('    Plot limits (sec, sec, uV, uV) [%g,%g,%g,%g]\n\n',xmin,xmax,ymin,ymax);

%
%%%%%%%%%%%%%%%%% Plot the envelope of the summed selected components %%%%%%%%%%%%%%%%%
%
if BOLD_COLORS==1
  mapcolors = 1:ntopos+1;
else
  mapcolors = [1 maporder+1];
end

if strcmpi(g.sumenv,'on')  | strcmpi(g.sumenv,'fill')
 sumenv = envelope(sumproj(:,:), g.envmode);
 if ~ylimset & max(sumenv) > ymax, ymax = max(curenv); end
 if ~ylimset & min(sumenv) < ymin, ymin = min(curenv); end
 if strcmpi(g.sumenv,'fill')  
    %
    % Plot the summed projection filled 
    %
    mins = matsel(sumenv,frames,0,2,0);
    p=fill([times times(frames:-1:1)],...
        [matsel(sumenv,frames,0,1,0) mins(frames:-1:1)],FILLCOLOR);
    set(p,'EdgeColor',FILLCOLOR);
    hold on
    %
    % Overplot the data envelope again so it is not covered by the fill()'d component
    %
    p=plot(times,matsel(envdata,frames,0,1,1),colors(mapcolors(1),1));% plot the max
    set(p,'LineWidth',2);                % component order (if BOLD_COLORS==0)
    p=plot(times,matsel(envdata,frames,0,2,1),colors(mapcolors(1),1));% plot the min
    set(p,'LineWidth',2);                % component order (if BOLD_COLORS==0)

 else % if no 'fill'
    tmp = matsel(sumenv,frames,0,2,0);
    p=plot(times,tmp);% plot the min
    hold on
    set(p,'color',FILLCOLOR);
    set(p,'linewidth',2);
    p=plot(times,matsel(sumenv,frames,0,1,0));% plot the max
    set(p,'linewidth',2);
    set(p,'color',FILLCOLOR);
 end
end
if strcmpi(g.pvaf,'pvaf')
    t = text(double(xmin+0.1*(xmax-xmin)), ...
             double(ymin+0.1*(ymax-ymin)), ...
             ['pvaf ' num2str(sumpvaf,'%4.2f') '%']);
   set(t,'fontsize',12,'fontweight','bold')
elseif strcmpi(g.pvaf,'rv') | strcmpi(g.pvaf,'off')
    t = text(double(xmin+0.1*(xmaxxpmin)), ...
             double(ymin+0.1*(ymax-ymin)), ...
             ['rv ' num2str(sumpvaf,'%4.2f') '%']);
   set(t,'fontsize',12,'fontweight','bold')
end

%
% %%%%%%%%%%%%%%%%%%%%%%%% Plot the computed component envelopes %%%%%%%%%%%%%%%%%%
%
    envx = [1;compx+1];
    for c = 1:ntopos+1   
        curenv = matsel(envdata,frames,0,1,envx(c));
        if ~ylimset & max(curenv) > ymax, ymax = max(curenv); end
        p=plot(times,curenv,colors(mapcolors(c),1));% plot the max
        set(gca,'FontSize',12,'FontWeight','Bold')
        if c==1                                % Note: use colors in original
            set(p,'LineWidth',2);              %       component order (if BOLD_COLORS==0)
        else
            set(p,'LineWidth',1);
        end
        if all_bold > 0
                set(p,'LineStyle','-','LineWidth',3);
        elseif mapcolors(c)>15                            % thin/dot 16th-> comp. envs.
                set(p,'LineStyle',':','LineWidth',1);
        elseif mapcolors(c)>10                            % 
                set(p,'LineStyle',':','LineWidth',2);
        elseif mapcolors(c)>6                             % dot 6th-> comp. envs.
                set(p,'LineStyle',':','LineWidth',3);
        elseif mapcolors(c)>1
            set(p,'LineStyle',colors(mapcolors(c),2),'LineWidth',1);
            if colors(mapcolors(c),2) == ':'
                set(l1,'LineWidth',2);  % embolden dotted env lines
            end
        end
        hold on
        curenv = matsel(envdata,frames,0,2,envx(c));
        if ~ylimset & min(curenv) < ymin, ymin = min(curenv); end
        p=plot(times,curenv,colors(mapcolors(c),1));% plot the min

        if c==1
            set(p,'LineWidth',2);
        else
            set(p,'LineWidth',1);
        end
        if all_bold > 0
                set(p,'LineStyle','-','LineWidth',3);
        elseif mapcolors(c)>15                            % thin/dot 11th-> comp. envs.
                set(p,'LineStyle',':','LineWidth',1);
        elseif mapcolors(c)>10                            
                set(p,'LineStyle',':','LineWidth',2);
        elseif mapcolors(c)>6                             % dot 6th-> comp. envs.
                set(p,'LineStyle',':','LineWidth',3);
        elseif mapcolors(c)>1
            set(p,'LineStyle',colors(mapcolors(c),2),'LineWidth',1);
            if colors(mapcolors(c),2) == ':'
                set(l1,'LineWidth',2);  % embolden dotted env lines
            end
        end
        if c==1 & ~isempty(g.vert)
            for v=1:length(g.vert)
                vl=plot([g.vert(v) g.vert(v)], [-1e10 1e10],'k--'); % plot specified vertical lines
                if any(g.limcontrib ~= 0) & v>= length(g.vert)-1;
                    set(vl,'linewidth',LIMCONTRIBWEIGHT);
                    set(vl,'linestyle',':');
                else
                    set(vl,'linewidth',VERTWEIGHT);
                    set(vl,'linestyle','--');
                end
            end
        end
        if g.limits(1) <= 0 & g.limits(2) >= 0    % plot vertical line at time zero
                vl=plot([0 0], [-1e10 1e10],'k');
                    set(vl,'linewidth',2);
        end
 
        %
        % plot the n-th component filled 
        %
        if g.fillcomp(1)>0 & find(g.fillcomp==c-1) 
            fprintf('filling the envelope of component %d\n',c-1);
            mins = matsel(envdata,frames,0,2,envx(c));
            p=fill([times times(frames:-1:1)],...
                   [matsel(envdata,frames,0,1,envx(c)) mins(frames:-1:1)],...
                   colors(mapcolors(c),1));
            %
            % Overplot the data envlope again so it is not covered by the fill()'d component
            %
            p=plot(times,matsel(envdata,frames,0,1,envx(1)),colors(mapcolors(1),1));% plot the max
            set(p,'LineWidth',2);                % component order (if BOLD_COLORS==0)
            p=plot(times,matsel(envdata,frames,0,2,envx(1)),colors(mapcolors(1),1));% plot the min
            set(p,'LineWidth',2);                % component order (if BOLD_COLORS==0)
        end
    end  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%%%%%%%%%%%%%%%%%%%%%%% Extend y limits by 5% %%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~ylimset
  datarange = ymax-ymin;
  ymin = ymin-0.05*datarange;
  ymax = ymax+0.05*datarange;
end
axis([pmin pmax ymin ymax]);

%
%%%%%%%%%%%%%%%%%%%%%% Label axes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
set(axe,'Color',axcolor);
if ~xunitframes
   l= xlabel('Time (s)');
else % xunitframes == 1
   l= xlabel('Data (time points)');
end
set(l,'FontSize',14,'FontWeight','Bold');
if strcmpi(g.envmode, 'avg')
    l=ylabel('Potential (uV)');
else 
    l=ylabel('RMS of uV');
end;    
set(l,'FontSize',14,'FontWeight','Bold');
%
%%%%%%%%%%%%%% Draw maps and oblique/vertical lines %%%%%%%%%%%%%%%%%%%%%
%
% axall = axes('Units','Normalized','Position',pos,...
axall = axes('Position',pos,...
    'Visible','Off','Fontsize',16); % whole-figure invisible axes
axes(axall)
set(axall,'Color',axcolor);
axis([0 1 0 1])

width  = xmax-xmin;
pwidth  = pmax-pmin;
height = ymax-ymin;

if strcmpi(g.dispmaps, 'on')
    for t=1:ntopos % draw oblique lines from max env vals (or plot top)
                 % to map bases, in left to right order
        %
        %%%%%%%%%%%%%%%%%%% draw oblique lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if BOLD_COLORS==1
            linestyles = 1:ntopos;
        else
            linestyles = maporder;
        end
        axes(axall) 
        axis([0 1 0 1]);
        set(axall,'Visible','off');
        maxenv = matsel(envdata,frames,plotframes(t),1,compx(t)+1); 
        % max env val
        data_y = 0.6*(g.voffsets(t)+maxenv-ymin)/height;
        if (data_y > pos(2)+0.6*pos(4)) 
            data_y = pos(2)+0.6*pos(4);
        end
        l1 = plot([(plottimes(t)-pmin)/pwidth  ...
                   topoleft + 1/pos(3)*(t-1)*1.2*topowidth + (topowidth*0.6)],...
                  [data_y 0.68], ...
                  colors(linestyles(t)+1)); % 0.68 is bottom of topo maps
        if all_bold > 0
                set(l1,'LineStyle','-','LineWidth',3);
        elseif linestyles(t)>15                        % thin/dot 11th-> comp. envs.
                set(l1,'LineStyle',':','LineWidth',1);
        elseif linestyles(t)>10 
                set(l1,'LineStyle',':','LineWidth',2);
        elseif linestyles(t)>5                     % dot 6th-> comp. envs.
                set(l1,'LineStyle',':','LineWidth',3);
        elseif linestyles(t)>1
            set(l1,'LineStyle',colors(linestyles(t)+1,2),'LineWidth',1);
            if colors(linestyles(t)+1,2) == ':'
                set(l1,'LineStyle',colors(linestyles(t)+1,2),'LineWidth',2);
            end
        end
        hold on
        %
        %%%%%%%%%%%%%%%%%%%% add specified vertical lines %%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if g.voffsets(t) > 0                    
            l2 = plot([(plottimes(t)-xmin)/width  ...
                       (plottimes(t)-xmin)/width],...
                      [0.6*(maxenv-ymin)/height ...
                       0.6*(g.voffsets(t)+maxenv-ymin)/height],...
                      colors(linestyles(t)+1));
            if all_bold > 0
                    set(l2,'LineStyle','-','LineWidth',3);
            elseif linestyles(t)>15                      % thin/dot 11th-> comp. envs.
                    set(l2,'LineStyle',':','LineWidth',1);
            elseif linestyles(t)>10                   
                    set(l2,'LineStyle',':','LineWidth',2);
            elseif linestyles(t)>5                   % dot 6th-> comp. envs.
                    set(l2,'LineStyle',':','LineWidth',3);
            else
                set(l1,'LineStyle',colors(linestyles(t)+1,2),'LineWidth',1);
                if colors(linestyles(t)+1,2) == ':'
                    set(l1,'LineWidth',2);
                end
            end
        end
        set(gca,'Visible','off');
        axis([0 1 0 1]);
    end % t
end; % if g.dispmaps == on

%
%%%%%%%%%%%%%%%%%%%%%%%%% Plot the topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%
%

if strcmpi(g.dispmaps, 'on')
    
    % common scale for colors
    % -----------------------
    if strcmpi(g.actscale, 'on')
        maxvolt = 0;
        for n=1:ntopos
            maxvolt = max(max(abs(maxproj(:,n))), maxvolt);
        end;
    end;
    
    [tmp tmpsort] = sort(maporder);
    [tmp tmpsort] = sort(tmpsort);
    if isstr(g.chanlocs)
        if exist(g.chanlocs) ~= 2  % if no such file
            fprintf('envtopo(): named channel location file ''%s'' not found.\n',g.chanlocs);
            return
        end
        eloc = readlocs(g.chanlocs);
        if length(eloc) ~= chans
            fprintf(...
              'envtopo(): locations for the %d data channels not in the channel location file.\n',chans);
            return
        end
    end
    chanlocs = g.chanlocs(g.plotchans); % topoplot based on g.plotchans only! -sm 11/04

    for t=1:ntopos % left to right order 
                   % axt = axes('Units','Normalized','Position',...
        axt = axes('Units','Normalized','Position',...
                   [pos(3)*topoleft+pos(1)+(t-1)*head_sep*topowidth pos(2)+0.66*pos(4) ...
                    topowidth topowidth*head_sep]);
        axes(axt)                             % topoplot axes
        cla
        
        if ~isempty(chanlocs)
            if ~isempty(varargin) 
                figure(myfig);topoplot(maxproj(:,t),chanlocs, varargin{:}); 
            else 
                figure(myfig);topoplot(maxproj(:,t),chanlocs,'style','both','emarkersize',3);
            end
            axis square
            if strcmpi(g.pvaf, 'on')
                set(gca, 'userdata', ['text(-0.6, -0.6, ''pvaf: ' sprintf('%6.2f', pvaf(tmpsort(t))) ''');'] );
            else
                set(gca, 'userdata', ['text(-0.6, -0.6, ''rv: ' sprintf('%6.2f', pvaf(tmpsort(t))) ''');'] );
            end;
        else axis off;
        end;

        %
        %%%%%%%%%%%%% Scale colors %%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if strcmpi(g.actscale, 'on')
            caxis([-maxvolt maxvolt]);
        end;
        %
        %%%%%%%%%%%%%%%%%%%%%%%% label components %%%%%%%%%%%%%%%%%%%%%%%
        %
        if t==1
            chid = fopen('envtopo.labels','r');
            if chid <3,
                numlabels = 1;
            else
                fprintf('Will label scalp maps with labels from file %s\n','envtopo.labels');
                compnames = fscanf(chid,'%s',[4 MAXPLOTDATACHANS]);
                compnames = compnames';
                [r c] = size(compnames);
                for i=1:r
                    for j=1:c
                        if compnames(i,j)=='.',
                            compnames(i,j)=' ';
                        end;
                    end;
                end;
                numlabels=0;
            end
        end
        if numlabels == 1
            complabel = int2str(maporder(t));        % label comp. numbers
        else
            complabel = compnames(t,:);              % use labels in file
        end
        text(0.00,0.80,complabel,'FontSize',14,...
             'FontWeight','Bold','HorizontalAlignment','Center');
        % axt = axes('Units','Normalized','Position',[0 0 1 1],...
        axt = axes('Position',[0 0 1 1],...
                   'Visible','Off','Fontsize',16);
        set(axt,'Color',axcolor);           % topoplot axes
        drawnow
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Plot a colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % axt = axes('Units','Normalized','Position',[.88 .58 .03 .10]);
    axt = axes('Position',[pos(1)+pos(3)*1.02 pos(2)+0.6*pos(4) pos(3)*.02 pos(4)*0.09]);
    if strcmpi(g.actscale, 'on')
        h=cbar(axt, [1:64],[-maxvolt maxvolt],3);
    else
        h=cbar(axt);                        % colorbar axes
        set(h,'Ytick',[]);
        
        axes(axall)
        set(axall,'Color',axcolor);
        tmp = text(0.50,1.05,g.title,'FontSize',16,'HorizontalAlignment','Center','FontWeight','Bold');
        set(tmp, 'interpreter', 'none');
        text(0.98,0.68,'+','FontSize',16,'HorizontalAlignment','Center');
        text(0.98,0.62,'-','FontSize',16,'HorizontalAlignment','Center');
    end;
    axes(axall)
    set(axall,'layer','top'); % bring component lines to top
    
end;
%
%%%%%%%%%%%%%%%%%%%%%%%%% turn on axcopy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end;');

return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function envdata = envelope(data, envmode)  % also in release as env()
  if nargin < 2
      envmode = 'avg';
  end;
  if strcmpi(envmode, 'rms');
      warning off;
      negflag = (data < 0);
      dataneg = negflag.* data;
      dataneg = -sqrt(sum(dataneg.*dataneg,1) ./ sum(negflag,1));
      posflag = (data > 0);
      datapos = posflag.* data;
      datapos = sqrt(sum(datapos.*datapos,1) ./ sum(posflag,1)); 
      envdata = [datapos;dataneg];
      warning on;
  else    
      if size(data,1)>1
          maxdata = max(data); % max at each time point
          mindata = min(data); % min at each time point
          envdata = [maxdata;mindata];
      else
          maxdata = max([data;data]); % max at each time point
          mindata = min([data;data]); % min at each time point
          envdata = [maxdata;mindata];
      end
  end;
  
return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
