%
% envtopo() - Plot the envelope of a multichannel data epoch, plus envelopes and 
%             scalp maps of specified or largest-contributing components. If a 3-D 
%             input matrix, operates on the mean of the data epochs. Click on 
%             individual axes to examine them in detail. The black lines represent 
%             the max and min values across all channels at each time point. 
%             The blue shading represents the max and min contributions of 
%             the selected components tothose channels. The paired colored lines 
%             represent the max and min contributions of the indicated component 
%             across all channels.
% Usage:
%             >> envtopo(data,weights,'chanlocs',file_or_struct);
%             >> [cmpvarorder,cmpvars,cmpframes,cmptimes,cmpsplotted,sortvar] ...
%                                     = envtopo(data, weights, 'key1', val1, ...);
% Inputs:
%  data        = single data epoch (chans,frames), else it a 3-D epoched data matrix 
%                (chans,frames,epochs) ->  processes the mean data epoch. 
%  weights     = linear decomposition (unmixing) weight matrix. The whole matrix 
%                should be passed to the function here.  
%                    Ex: (EEG.icaweights*EEG.icasphere)
% Required keyword:
%  'chanlocs'  = [string] channel location filename or EEG.chanlocs structure. 
%                For more information, see >> topoplot example 
%
% Optional inputs:
%  'compnums'  = [integer array] vector of component indices to use in the 
%                  calculations and to select plotted components from 
%                  {default|[]: all}
%  'compsplot' = [integer] the number of largest contributing components to plot.
%                  compnums in the latter case is restricted in size by the 
%                  internal variable MAXTOPOS (20) {default|[] -> 7}
%  'subcomps'  = [integer vector] indices of comps. to remove from the whole data 
%                  before plotting. 0 -> Remove none. [] -> If 'compnums' 
%                  also listed, remove *all* except 'compnums' {default: 0}
%  'limits'    = 0 or [minms maxms] or [minms maxms minuV maxuV]. Specify 
%                  start/end plot (x) limits (in ms) and min/max y-axis limits 
%                  (in uV). If 0, or if both minmx & maxms == 0 -> use latencies 
%                  from 'timerange' (else 0:frames-1). If both minuV and 
%                  maxuV == 0 -> use data uV limits {default: 0}
%  'timerange' = start and end input data latencies (in ms) {default: from 'limits' 
%                  if any}. Note: Does NOT select a portion of the input data, 
%                  just makes time labels.
%  'limcontrib' = [minms maxms]  time range (in ms) in which to rank component 
%                  contribution (boundaries shown with thin dotted lines) 
%                  {default|[]|[0 0] -> plotting limits}
%  'sortvar'   = ['mp'|'pv'|'pp'|'rp'] {default:'mp'} 
%                  'mp', sort components by maximum mean back-projected power 
%                  in the 'limcontrib' time range: mp(comp) = max(Mean(back_proj.^2));
%                  where back_proj = comp_map * comp_activation(t) for t in 'limcontrib'
%                  'pv', sort components by percent variance accounted for (eeg_pvaf())
%                    pvaf(comp) = 100-100*mean(var(data - back_proj))/mean(var(data));
%                  'pp', sort components by percent power accounted for (ppaf) 
%                    ppaf(comp) = 100-100*Mean((data - back_proj).^2)/Mean(data.^2);
%                  'rp', sort components by relative power 
%                    rp(comp) = 100*Mean(back_proj.^2)/Mean(data.^2);
%  'title'     = [string] plot title {default|[] -> none}
%  'plotchans' = [integer array] data channels to use in computing contributions and 
%                  envelopes, and also for making scalp topo plots
%                  {default|[] -> all}, by calling topoplot().
%  'voffsets'  = [float array] vertical line extentions above the data max to 
%                  disentangle plot lines (left->right heads, values in y-axis units) 
%                  {default|[] -> none}
%  'colors'    = [string] filename of file containing colors for envelopes, 3 chars
%                  per line, (. = blank). First color should be "w.." (white)
%                  Else, 'bold' -> plot default colors in thick lines.
%                  {default|[] -> standard Matlab color order}
%  'fillcomp'  = int_vector>0 -> fill the numbered component envelope(s) with 
%                  solid color. Ex: [1] or [1 5] {default|[]|0 -> no fill}
%  'vert'      = vector of times (in ms) at which to plot vertical dashed lines 
%                  {default|[] -> none}
%  'icawinv'   = [float array] inverse weight matrix. Normally computed by inverting
%                  the weights*sphere matrix (Note: If some components have been 
%                  removed, the pseudo-inverse does not represent component maps 
%                  accurately).
%  'icaact'    = [float array] component activations. {default: computed from the 
%                  input weight matrix}
%  'envmode'   = ['avg'|'rms'] compute the average envelope or the root mean square
%                  envelope {default: 'avg'}
%  'sumenv'    = ['on'|'off'|'fill'] 'fill' -> show the filled envelope of the summed 
%                  projections of the selected components; 'on' -> show the envelope 
%                  only {default: 'fill'}
%  'actscale'  = ['on'|'off'] scale component scalp maps by maximum component activity 
%                  in the designated (limcontrib) interval. 'off' -> scale scalp maps 
%                  individually using +/- max(abs(map value)) {default: 'off'}
%  'dispmaps'  = ['on'|'off'] display component numbers and scalp maps {default: 'on'}
%  'topoplotkey','val' = optional additional topoplot() arguments {default: none}
%
% Outputs:
% 
%  cmpvarorder = component numbers in decreasing order of max variance in data
%  cmpvars     = ('sortvar') max power for all components tested, sorted 
%                   from highest to lowest. See 'sortvar' 'mp'
%  cmpframes   = frames of comvars for each component plotted
%  cmptimes    = times of compvars for each component plotted
%  cmpsplotted = indices of components plotted. 
%                       unsorted_compvars = compvars(compsplotted);
%  sortvar      = The computed data used to sort components with. 
%                 See 'sortvar' option above.
%
% Notes:
%  To label maps with other than component numbers, put four-char strings into 
%  a local (pwd) file named 'envtopo.labels' (using . = space) in time-order 
%  of their projection maxima
%
% Authors: Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 3/1998 
%
% See also: timtopo() axcopy()

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

function [compvarorder,compvars,compframes,comptimes,compsplotted,sortvar] = envtopo(data,weights,varargin)

% icadefs;    % read toolbox defaults

sortvar = []; 
all_bold = 0;
BOLD_COLORS = 1;    % 1 = use solid lines for first 5 components plotted
                    % 0 = use std lines according to component rank only
FILL_COMP_ENV = 0;  % default no fill
FILLCOLOR   = [.815 .94 1]; % use lighter blue for better env visibility
% FILLCOLOR = [.66 .76 1];
MAXTOPOS = 20;      % max topoplots to plot
VERTWEIGHT = 2.0;  % lineweight of specified vertical lines
LIMCONTRIBWEIGHT = 1.2; % lineweight of limonctrib vertical lines
MAX_FRAMES = 10000; % maximum frames to plot
MAXENVPLOTCHANS = 10;  
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
				  'figure'        'integer'  []                       [] ;                   
				  'colorfile'     'string'   []                       '' ; 
				  'colors'        'string'   []                       '' ;
				  'compnums'      'integer'  []                       []; 
				  'compsplot'     'integer'  []                       7; 
				  'subcomps'      'integer'  []                       0; 
				  'envmode'       'string'   {'avg','rms'}            'avg'; 
				  'dispmaps'      'string'   {'on','off'}             'on'; 
				  'pvaf'          'string'   {'mp','mv','on','rp','rv','pv','pvaf','pp','off',''} ''; 
				  'sortvar'       'string'   {'mp','mv','rp','rv','pv','pvaf','pp'} 'mp';  
				  'actscale'      'string'   {'on','off'}             'off'; 
				  'limcontrib'    'real'     []                       0;  
				  'topoarg'    'real'     []                       0;  
				  'sumenv'        'string'    {'on','off','fill'}     'fill'};

	% Note: Above, previous 'pvaf' arguments 'on' -> 'pv', 'off' -> 'rv'
	%       for backwards compatibility 11/2004 -sm
	
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
    g.sortvar = 'mp';
    g.pvaf = []; 
    g.timerange = [];
    g.icaact = [];
    g.limcontrib = 0;
    g.icawinv = pinv(weights);
    g.subcomps = 0;
    g.envmode = 'avg';
    g.dispmaps = 'on';
end;
if ~isempty(g.figure), figure(g.figure); end;        % remember the current figure (for Matlab 7.0.0 bug)

if ~isempty(g.pvaf) 
	g.sortvar = g.pvaf; % leave deprecated g.pvaf behind. 
end

if strcmpi(g.sortvar,'on') | strcmpi(g.sortvar,'pvaf') | strcmpi(g.sortvar,'mv')
   g.sortvar = 'mp'; % update deprecated flags
end
if strcmpi(g.sortvar,'off') | strcmp(g.sortvar,'rv') 
   g.sortvar = 'rp';
end


%
% Check input flags and arguments
% 
if ndims(data) == 3
    data = mean(data,3); % average the data if 3-D
end;
[chans,frames] = size(data); 

if frames > MAX_FRAMES
   error('number of data frames to plot too large!');
end

if isstr(g.chanlocs)
    g.chanlocs = readlocs(g.chanlocs);  % read channel location information
    if length(g.chanlocs) ~= chans
     fprintf(...
      'envtopo(): locations for the %d data channels not in the channel location file.\n', ...
        chans);
        return
    end
end
% Checking dimension of chanlocs
if length(g.chanlocs) > size(data,1)
    fprintf(2,['\n envtopo warning: Dimensions of data to plot and channel location are not consistent\n' ...
        '                 As a result, the plots generated can be WRONG\n'...
        '                 If you are providing the path to your channel location file. Make sure\n'...
        '                 to load the file first, and to provide as an input only the channels that\n'...
        '                 correspond to the data to be plotted, otherwise just provide a valid chanlocs\n ']);
elseif length(g.chanlocs) > size(data,1)
        fprintf(2,['\n envtopo error: Dimensions of data to plot and channel location are not consistent\n' ...
        '                 Check the dimension of the channel location provided\n'...
        '                 Aborting plot ...\n']);
    exit;
end

if ~isempty(g.colors)
    g.colorfile = g.colors; % retain old usage 'colorfile' for 'colors' -sm 4/04
end

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
if pmin < xmin
   pmin = xmin;     % don't allow plotting beyond the data limits
end
pmax = g.limits(2);
if pmax > xmax
   pmax = xmax;
end

dt = (xmax-xmin)/(frames-1);  % sampling interval in sec
times=xmin*ones(1,frames)+dt*(0:frames-1); % time points in sec

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
if wtcomps ~= chans
   fprintf('Number of components not the same as number of channels.\n'); 
   fprintf('  - component scalp maps and time courses may not be correct.\n');
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

if g.compsplot < 0
   g.compsplot = abs(g.compsplot);
end

if g.compnums < 0 % legacy syntax
   g.compsplot = abs(g.compnums);
   g.compnums = [];
end
if isempty(g.compnums) | g.compnums(1) == 0
    g.compnums = 1:wtcomps; % by default, select from all components
end

if g.compsplot > MAXTOPOS
    fprintf('Can only plot a maximum of %d components.\n',MAXTOPOS);
    return
else
    MAXTOPOS = g.compsplot;
end

if max(g.compnums) > wtcomps | min(g.compnums)< 1
    error('Keyword ''compnums'' out of range (1 to %d)', wtcomps);
end

g.subcomps = abs(g.subcomps); % don't pass negative channel numbers
if max(g.subcomps) > wtcomps
    error('Keyword ''subcomps'' argument out of bounds');
end

%
%%%%%%%%%%%%%%% Subtract components from data if requested %%%%%%%%%%%%%
%

ncomps = length(g.compnums);

if isempty(g.subcomps) % remove all but compnums
     g.subcomps = 1:wtcomps;
     g.subcomps(g.compnums) = [];
else
  % g.subcomps 0    -> subtract no comps
  %            list -> subtract subcomps list
  if min(g.subcomps) < 1
      if length(g.subcomps) > 1
         error('Keyword ''subcomps'' argument incorrect.');
      end
      g.subcomps = [];   % if subcomps contains a 0 (or <1), don't remove components
  elseif max(g.subcomps) > wtcomps
      error('Keyword ''subcomps'' argument out of bounds.');
  end
end

g.icaact = weights*data;
if ~isempty(g.subcomps)
  fprintf('Subtracting requested components from plotting data: ');
  for k = 1:length(g.subcomps)
      fprintf('%d ',g.subcomps(k));
      if ~rem(k,32)
         fprintf('\n');
      end
  end
  fprintf('\n');
  g.icaact(g.subcomps,:) = zeros(length(g.subcomps),size(data,2));
end;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Process components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

for i=1:ncomps-1
    if g.compnums(i) == 0
	       fprintf('Removing component number 0 from compnums.\n');
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

% toby 2.16.2006: maxproj will now contain all channel info, in case
%                 plotgrid is called in topoplot.
%                 NOT maxproj = zeros(length(g.plotchans),ncomps);
maxproj = zeros(chans,ncomps);

%
% first, plot the data envelope
%
envdata = zeros(2,frames*(ncomps+1));
envdata(:,1:frames) = envelope(g.icawinv(g.plotchans,:)*g.icaact, g.envmode); 

fprintf('Data epoch is from %.0f ms to %.0f ms.\n',1000*xmin,1000*xmax);
fprintf('Plotting data from %.0f ms to %.0f ms.\n',1000*xmin,1000*xmax);
fprintf('Comparing maximum projections for components:  \n');
        if ncomps>32
           fprintf('\n');
        end

compvars = zeros(1,ncomps);
mapsigns = zeros(1,ncomps);

%
% Compute frames to plot
%
sampint  = (xmax-xmin)/(frames-1);     % sampling interval in sec
times    = xmin:sampint:xmax;          % make vector of data time values

[v minf] = min(abs(times-pmin));
[v maxf] = min(abs(times-pmax));
pframes  = minf:maxf;         % frames to plot
ptimes   = times(pframes);    % times to plot
if limframe1 < minf
   limframe1 = minf;
end
if limframe2 > maxf
   limframe2 = maxf;
end

%
%%%%%%%%%%%%%% find max variances and their frame indices %%%%%%%%%%%
%

if strcmp(g.sortvar,'pv') %Changed -Jean
	% Variance of the data in the interval, for calculating sortvar. 
	vardat = mean(var(data(g.plotchans,limframe1:limframe2),1));
else 
	% Compute data rms for sortvar
	powdat = mean(mean(data(g.plotchans,limframe1:limframe2).^2));
end

for c = 1:ncomps
    
      if isempty(g.icaact) % make the back-projection of component c

          % Changed to include all channels in computation for use in 
          % topoplot, particularly with plotgrid option. toby 2.16.2006
          proj = g.icawinv(:,g.compnums(c))*weights(g.compnums(c),:)*data;
      else 
          proj = g.icawinv(:,g.compnums(c))*g.icaact(g.compnums(c),:);
      end;    

      % save the comp envelope for plotting component waveforms
      envdata(:,c*frames+1:(c+1)*frames) = envelope(proj(g.plotchans,:), g.envmode); 

      % Find the frame (timepoint) of largest rms component value 
      % and the relative value to those channels defined by plotchans. 
      if length(g.plotchans) > 1
        [maxval,maxi] = max(mean((proj(g.plotchans,limframe1:limframe2)).^2));
      else % g.plotchans == 1 --> find absmax value
        [maxval,maxi] = max((abs(proj(g.plotchans,limframe1:limframe2)))); 
      end
      maxi = maxi+limframe1-1;

      % plotframes and compvars are needed for plotting the lines indicating 
      % the timepoint a topoplot refers to.
      plotframes(c) = maxi;
      compvars(c)   = maxval;       % find value of max variance for comp c
      maxproj(:,c)  = proj(:,maxi); % maxproj now contains all channels, to handle 
                                    % the 'plotchans'/topoplot'gridplot' conflict. 
									% Toby 2.17.2006
      %
      %%%%%% Calculate sortvar, used to sort the components %%%%%%%%%%%
      %
      if strcmpi(g.sortvar,'mp')  % Maximum Power of backproj
          sortvar(c) = maxval;
          fprintf('IC%1.0f maximum mean power of back-projection: %g\n',c,sortvar(c));
          
      elseif strcmpi(g.sortvar, 'pv')   % Percent Variance
          % toby 2.19.2006: Changed to be consistent with eeg_pvaf().
          sortvar(c) = 100-100*mean(var(data(g.plotchans,limframe1:limframe2)...
                        - proj(g.plotchans,limframe1:limframe2),1))/vardat;
          fprintf('IC%1.0f percent variance accounted for(pvaf): %2.2f%%\n',c,sortvar(c));

      elseif strcmpi(g.sortvar,'pp')    % Percent Power
          sortvar(c) = 100-100*mean(mean((data(g.plotchans,limframe1:limframe2)...
                        - proj(g.plotchans,limframe1:limframe2)).^2))/powdat;
          fprintf('IC%1.0f percent power accounted for(ppaf): %2.2f%%\n',c,sortvar(c));
          
      elseif strcmpi(g.sortvar,'rp')    % Relative Power
          sortvar(c) = 100*mean(mean((proj(g.plotchans,limframe1:limframe2)).^2))/powdat;
          fprintf('IC%1.0f relative power of back-projection: %g\n',c,sortvar(c));
      else
          error('''sortvar'' argument unknown');
      end;

end % component c
fprintf('\n');

%
%%%%%%%%%%%%%%% Compute component selection criterion %%%%%%%%%%%%%%%%%%%%%%%%%%
%

% compute sortvar
if ~xunitframes
  fprintf('  in the interval %3.0f ms to %3.0f ms.\n',...
					1000*times(limframe1),1000*times(limframe2));
end

[sortsortvar spx] = sort(sortvar);
sortsortvar = sortsortvar(end:-1:1);
spx      = spx(end:-1:1);
npercol = ceil(ncomps/3);

%
%%%%%%%%%%%%%%%%%%%%%%%%% Sort the components %%%%%%%%%%%%%%%%%%%%%%%%%%%
%

[tmp,compx]  = sort(sortvar');           % sort compnums on sortvar (as defined by input 
                                      % 'sortvar', default is 'mp').
compx        = compx(ncomps:-1:1);    % reverse order of sort
compvars     = compvars(ncomps:-1:1)';% reverse order of sort (output var)
compvarorder = g.compnums(compx);     % actual component numbers (output var)
plotframes   = plotframes(compx);     % plotted comps have these max frames 
compframes   = plotframes';           % frame of max variance in each comp (output var)
comptimes    = times(plotframes(compx)); % time of max variance in each comp (output var)
compsplotted = compvarorder(1:ntopos); % (output var)

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

fprintf('    Component sortvar in interval:  ');
for t=1:ntopos
     fprintf('%4.2f ',sortvar(t));
end
fprintf('\n');

sumproj = zeros(size(data(g.plotchans,:))); % toby 2.21.2006 REDUNDANT CALCULATIONS!
for n = 1:ntopos
  if isempty(g.icaact)
      sumproj = sumproj + ...
		g.icawinv(g.plotchans,maporder(n))*weights(maporder(n),:)*data; 
  else 
      sumproj = sumproj + g.icawinv(g.plotchans,maporder(n))*g.icaact(maporder(n),:);     
       % updated -sm 11/04
  end; % Note: sumproj here has only g.plotchans
end
rmsproj = mean(mean((data(g.plotchans,limframe1:limframe2).^2))); % find data rms in interval

if strcmpi(g.sortvar,'rp') 
 	sumppaf = mean(mean(sumproj(:,limframe1:limframe2).^2));      
    sumppaf = 100*sumppaf/rmsproj;
    ot   = 'rp';
elseif strcmpi(g.sortvar,'pv') 
    sumppaf = mean(var((data(g.plotchans,limframe1:limframe2) ...
                                  - sumproj(:,limframe1:limframe2)),1));
    sumppaf = 100-100*sumppaf/vardat;
    ot   = 'pvaf';
else
   sumppaf = mean(mean((data(g.plotchans,limframe1:limframe2) ...
                                  - sumproj(:,limframe1:limframe2)).^2)); 
    sumppaf = 100-100*sumppaf/rmsproj;
    ot   = 'ppaf';
end;

if ~xunitframes
   fprintf('    Summed component ''%s'' in interval [%4g %4g] ms: %4.2f%%\n',...
					ot, 1000*times(limframe1),1000*times(limframe2), sumppaf);
end

%
% Collect user-supplied Y-axis information
% Edited and moved here from 'Collect Y-axis information' section below -Jean
%
if length(g.limits) == 4 
     if g.limits(3)~=0 | g.limits(4)~=0 % collect plotting limits from 'limits'
	 ymin = g.limits(3);
	 ymax = g.limits(4);
         ylimset = 1;
     end
else
  ylimset = 0; % flag whether hard limits have been set by the user
  ymin = min(min(g.icawinv(g.plotchans,:)*g.icaact(:,pframes))); % begin by setting limits from plotted data
  ymax = max(max(g.icawinv(g.plotchans,:)*g.icaact(:,pframes)));
end

fprintf('    Plot limits (sec, sec, uV, uV) [%g,%g,%g,%g]\n\n',xmin,xmax, ymin,ymax);

%
%%%%%%%%%%%%%%%%%%%%% Plot the data envelopes %%%%%%%%%%%%%%%%%%%%%%%%%
%
BACKCOLOR = [0.7 0.7 0.7];
FONTSIZE=12;
FONTSIZEMED=10;
FONTSIZESMALL=8;
newaxes=axes('position',pos);
axis off
set(newaxes,'FontSize',FONTSIZE,'FontWeight','Bold','Visible','off');
set(newaxes,'Color',BACKCOLOR); % set the background color
delete(newaxes) %XXX

% site the plot at bottom of the current axes
axe = axes('Position',[pos(1) pos(2) pos(3) 0.6*pos(4)],...
           'FontSize',FONTSIZE,'FontWeight','Bold');

g.limits = get(axe,'Ylim');
set(axe,'GridLineStyle',':')
set(axe,'Xgrid','off')
set(axe,'Ygrid','on')
axes(axe)
set(axe,'Color',axcolor);

%
%%%%%%%%%%%%%%%%% Plot the envelope of the summed selected components %%%%%%%%%%%%%%%%%
%
if BOLD_COLORS==1
  mapcolors = 1:ntopos+1;
else
  mapcolors = [1 maporder+1];
end

if strcmpi(g.sumenv,'on')  | strcmpi(g.sumenv,'fill') %%%%%%%% if 'sunvenv' %%%%%%%%%
 sumenv = envelope(sumproj(:,:), g.envmode);
 if ~ylimset & max(sumenv) > ymax, ymax = max(sumenv(1,:)); end
 if ~ylimset & min(sumenv) < ymin, ymin = min(sumenv(2,:)); end
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
    % Overplot the data envelope so it is not covered by the fill()'d component
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

if strcmpi(g.sortvar,'rp')
	t = text(double(xmin+0.1*(xmax-xmin)), ...
         double(ymin+0.1*(ymax-ymin)), ...
         ['rp ' num2str(sumppaf,'%4.2f') '%']);
elseif strcmpi(g.sortvar,'pv')
    	t = text(double(xmin+0.1*(xmax-xmin)), ...
         double(ymin+0.1*(ymax-ymin)), ...
         ['pvaf ' num2str(sumppaf,'%4.2f') '%']);
else
	t = text(double(xmin+0.1*(xmax-xmin)), ...
         double(ymin+0.1*(ymax-ymin)), ...
         ['ppaf ' num2str(sumppaf,'%4.2f') '%']);
end
set(t,'fontsize',FONTSIZESMALL,'fontweight','bold')

%
% %%%%%%%%%%%%%%%%%%%%%%%% Plot the computed component envelopes %%%%%%%%%%%%%%%%%%
%
    envx = [1;compx+1];
    for c = 1:ntopos+1   
        curenv = matsel(envdata,frames,0,1,envx(c));
        if ~ylimset & max(curenv) > ymax, ymax = max(curenv); end
        p=plot(times,curenv,colors(mapcolors(c),1));% plot the max
        set(gca,'FontSize',FONTSIZESMALL,'FontWeight','Bold')
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
set(l,'FontSize',FONTSIZEMED,'FontWeight','Bold');
if strcmpi(g.envmode, 'avg')
    l=ylabel('Potential (uV)');
else 
    l=ylabel('RMS of uV');
end;    
set(l,'FontSize',FONTSIZEMED,'FontWeight','Bold');
%
%%%%%%%%%%%%%% Draw maps and oblique/vertical lines %%%%%%%%%%%%%%%%%%%%%
%
% axall = axes('Units','Normalized','Position',pos,...
axall = axes('Position',pos,...
    'Visible','Off','Fontsize',FONTSIZE); % whole-figure invisible axes
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

    for t=1:ntopos % left to right order  (maporder)
                   % axt = axes('Units','Normalized','Position',...
        axt = axes('Units','Normalized','Position',...
                   [pos(3)*topoleft+pos(1)+(t-1)*head_sep*topowidth pos(2)+0.66*pos(4) ...
                    topowidth topowidth*head_sep]);
        axes(axt)                             % topoplot axes
        cla
        
        if ~isempty(g.chanlocs)  % plot the component scalp maps
            if ~isempty(varargin)
                topoplot(maxproj(g.plotchans,t),g.chanlocs(g.plotchans), varargin{:});
            else  % if no varargin specified
                topoplot(maxproj(g.plotchans,t),g.chanlocs(g.plotchans),...
							'style','both','emarkersize',3);
            end
            axis square
            set(gca, 'userdata', ...
   ['text(-0.6, -0.6, ''' g.sortvar ': ' sprintf('%6.2f', sortvar(tmpsort(t))) ''');']);
        else 
			axis off;
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
                fprintf('Will label scalp maps with labels from file %s\n',...
					'envtopo.labels');
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
        text(0.00,0.80,complabel,'FontSize',FONTSIZEMED,...
             'FontWeight','Bold','HorizontalAlignment','Center');
        % axt = axes('Units','Normalized','Position',[0 0 1 1],...
        axt = axes('Position',[0 0 1 1],...
                   'Visible','Off','Fontsize',FONTSIZE);
        set(axt,'Color',axcolor);           % topoplot axes
        drawnow
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Plot a colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % axt = axes('Units','Normalized','Position',[.88 .58 .03 .10]);
    axt = axes('Position',[pos(1)+pos(3)*1.015 pos(2)+0.6055*pos(4) ...
				pos(3)*.02 pos(4)*0.09]);
    if strcmpi(g.actscale, 'on')
        h=cbar(axt, [1:64],[-maxvolt maxvolt],3);
    else
        h=cbar(axt);                        % colorbar axes
        set(h,'Ytick',[]);
        
        axes(axall)
        set(axall,'Color',axcolor);
        tmp = text(0.50,1.05,g.title,'FontSize',FONTSIZE,...
						'HorizontalAlignment','Center',...
						'FontWeight','Bold');
        set(tmp, 'interpreter', 'none');
        text(1,0.68,'+','FontSize',FONTSIZE,'HorizontalAlignment','Center');
        % text(1,0.637,'0','FontSize',12,'HorizontalAlignment','Center',...
		%	'verticalalignment','middle');
        text(1,0.61,'-','FontSize',FONTSIZE,'HorizontalAlignment','Center');
    end;
    axes(axall)
    set(axall,'layer','top'); % bring component lines to top
    
end;
%
%%%%%%%%%%%%%%%%%%%%%%%%% turn on axcopy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
axcopy(gcf, ...
  'if ~isempty(get(gca,''''userdata'''')), eval(get(gca,''''userdata''''));end;');

return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function envdata = envelope(data, envmode)  % also in release as env()
  if nargin < 2
      envmode = 'avg';
  end;
  if strcmpi(envmode, 'rms'); % return rms of pos and neg vals at each time point 
      warning off;
      datnaeg = (data < 0).*data; % zero out positive values
      dataneg = -sqrt(sum(dataneg.^2,1) ./ sum(negflag,1));

      datapos = (data > 0).*data; % zero out negative values
      datapos = sqrt(sum(datapos.^2,1) ./ sum(posflag,1)); 

      envdata = [datapos;dataneg];
      warning on;
  else  % find max and min value at each time point
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
