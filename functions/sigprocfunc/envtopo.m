% envtopo() - Plot a data epoch envelope with envelopes and scalp maps of 
%             specified components. Click on individual plots to examine
%             separately (with zoom feature).
% Usage:
%     >> [compvarorder,compvars,compframes,comptimes,compsplotted] ...
%           = envtopo(data,weights,'chan_locs',[limits],[compnums],...
%             'title',[plotchans],[voffsets],'colorfile',fill_comp,[vert]);
% Inputs:
%  data       = single data epoch (chans,frames)
%  weights    = final weight matrix from runica() (=weights*sphere)
%
% Optional inputs:
%  'chan_locs' = channel location file. See >> topoplot example {def|[]->see icadefs)
%  limits      = [xmin xmax ymin ymax]  x values in ms 
%                {def|[] or both y's 0 -> y data limits}
%  compnums    = vector of component numbers to plot {default|0 -> all}
%                ELSE n<0, the number largest-comp. maps to plot {default|[] -> 7}
%  'title'     = plot title {default|[] -> none}
%  plotchans   = data channels to use in computing envelopes {default|[] -> all}
%                (Can specify >= 1 data channels).
%  voffsets    = vert. line extentions above the data max to disentangle plot
%                lines (left->right heads, values in y-axis units) {def|[] -> none}
% 'colorfile'  = filename of file containing colors for envelopes, 3 chars
%                per line, (. = blank). First color should be "w.." (white)
%                Colorfile argument 'bold' uses default colors, all thick lines.
%                {default|[] -> standard color order}
%  fill_comp   = int_vector>0 -> fill the numbered component envelope(s) with 
%                solid color. Ex: [1] or [1 5] {default|[]|0 -> no fill}
%  vert        = vector of times to plot vertical lines {default|[] -> none}
%
% Outputs:
%  compvarorder  = component numbers in decreasing order of max variance in data
%  compvars      = component max variances
%  compframes    = frames of max variance
%  comptimes     = times of max variance
%  compsplotted  = components plotted
%
% Notes:
%  To label maps with other than component numbers, put 4-char strings into
%  file 'envtopo.labels' (. = space) in time-order of their projection maxima
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 3/1998 
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

function [compvarorder,compvars,compframes,comptimes,compsplotted] = envtopo(data,weights,chan_locs,limits,compnums,titl,plotchans,voffsets,colorfile,fill_comp_env,vert, varargin)

if nargin < 2
   help envtopo
   return
end

uraxes = gca; % the original figure or subplot axes
pos=get(gca,'Position');
axcolor = get(gca,'Color');
delete(gca)

all_bold = 0;
BOLD_COLORS = 1;  % 1 = use solid lines for first 5 components plotted
                  % 0 = use std lines according to component rank only
FILL_COMP_ENV = 0;  % default no fill
MAXTOPOS = 7;  % max topoplots to plot

[chans,frames] = size(data);
[wtcomps,wchans] = size(weights);
if wchans ~= chans
   fprintf('envtopo(): sizes of weights and data do not agree.\n');
   return
end

icadefs;    % read toolbox defaults
ENVCOLORS = strvcat('w..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..','m..','c..','r..','g..','b..');

if nargin < 11
    vert = [];
end
if nargin <10
    fill_comp_env = [];
end
if isempty(fill_comp_env) | fill_comp_env(1) == 0
    fill_comp_env = FILL_COMP_ENV; % default
end
if nargin < 9
    colorfile = [];
end
if isempty(colorfile) | colorfile(1) == 0
    colorfile = ENVCOLORS; % filename read from icadefs
end

if nargin < 8
  voffsets = [];
end
if isempty(voffsets) | ( size(voffsets) == [1,1] & voffsets(1) == 0 )
  voffsets = zeros(1,MAXTOPOS); 
end
if nargin < 7
   plotchans=[];
end
if isempty(plotchans) | plotchans(1)==0 
   plotchans = 1:chans;
end
if max(plotchans) > chans | min(plotchans) < 1
   help envtopo
   return
end

if nargin < 6,
   titl = [];
end
if isempty(titl) | titl(1) == 0
   titl = '';     % DEFAULT TITLE
end

if nargin < 5
    compnums = [];
end
if isempty(compnums) | compnums(1) == 0
    compnums = 1:wtcomps; % by default, all components
end
if min(compnums) < 0
  if length(compnums) > 1
     fprintf('envtopo(): negative compnums must be a single integer.\n');
     return
  end
  if -compnums > MAXTOPOS
    fprintf('Can only plot a maximum of %d components.\n',MAXTOPOS);
    return
  else
    MAXTOPOS = -compnums;
    compnums = 1:wtcomps;
  end
end
ncomps = length(compnums);
for i=1:ncomps-1
  for j=i+1:ncomps
    if compnums(i)==compnums(j)
       fprintf('Cannot repeat component number (%d) in compnums.\n',compnums(i));
       return
    end
  end
end

limitset = 0;
if nargin < 4,
  limits = 0;
elseif length(limits)>1
  limitset = 1;
end
if isempty(limits)
  limits = 0;
end

if nargin< 3
  chan_locs = [];
end
if isempty(chan_locs) | (isnumeric(chan_locs) & (chan_locs(1) == 0))
  chan_locs = DEFAULT_ELOC;
% chan_spline = [];  % for headplot() alternative (commented out below)
end

%
%%%%%%%%%%%%%%%%%%%% Read and adjust limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if limits==0,      % == 0 or [0 0 0 0]
  xmin=0;
  xmax=frames-1;
  ymin=min(min(data));
  ymax=max(max(data));
  datarange = ymax-ymin;
  ymin = ymin-0.05*datarange;
  ymax = ymax+0.05*datarange;
else
  if length(limits)~=4,
    fprintf('envtopo: limits should be 0 or an array [xmin xmax ymin ymax].\n');
    return
  end;
  xmin = limits(1);
  xmax = limits(2);
  ymin = limits(3);
  ymax = limits(4);
end;

if xmax == 0 & xmin == 0,
  x = (0:1:frames-1);
  xmin = 0;
  xmax = frames-1;
else
  dx = (xmax-xmin)/(frames-1);
  x=xmin*ones(1,frames)+dx*(0:frames-1); % compute x-values
end;
if xmax<=xmin,
  fprintf('envtopo() - xmax must be > xmin.\n')
  return
end

if ymax == 0 & ymin == 0,
  ymax=max(max(data));
  ymin=min(min(data));
  datarange = ymax-ymin;
  ymin = ymin-0.05*datarange;
  ymax = ymax+0.05*datarange;
end
if ymax<=ymin,
  fprintf('envtopo() - ymax must be > ymin.\n')
  return
end
%
%%%%%%%%%%%%%%%%%%%% Read the color names %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~isstr(colorfile)
  fprintf('envproj(): color file name must be a string.\n');
  return
end
if length(colorfile)== 4 & colorfile == 'bold'
   all_bold = 1;
   colorfile = ENVCOLORS; % filename read from icadefs
end
if length(colorfile(:)) < 50
	cid = fopen(colorfile,'r');
	if cid <3,
		fprintf('envproj(): cannot open file %s.\n',colorfile);
		return
	else
		colors = fscanf(cid,'%s',[3 MAXENVPLOTCHANS]);
		colors = colors';
	end;
else
	colors = colorfile;
end;
[r c] = size(colors);
for i=1:r
	for j=1:c
		if colors(i,j)=='.',
            if j==1
				fprintf(...
					'envtopo(): colors file should have color letter in 1st column.\n');
				return
            elseif j==2
				colors(i,j)='-';
            elseif j>2
				colors(i,j)=' ';
            end
		end;
	end;
end;
colors(1,1) = 'k'; % make sure 1st color (for data envelope) is black
% [rr cc] = size(colors);
% colors

%
%%%%%%%%%%%%%%% Compute plotframes and envdata %%%%%%%%%%%%%%%%%%%%%
%

ntopos = length(compnums);
if ntopos > MAXTOPOS
  ntopos = MAXTOPOS; % limit the number of topoplots to display
end

if max(compnums) > wtcomps | min(compnums)< 1
  fprintf(...
'envtopo(): one or more compnums out of range (1,%d).\n',wtcomps);
  return
end

if size(weights,1) == size(weights,2)
  winv = inv(weights);
else
  fprintf('Using pseudo inverse pinv().\n');
  winv = pinv(weights);
end
plotframes = ones(ncomps);
maxproj = zeros(chans,ncomps);
envdata = zeros(2,frames*(ncomps+1));
envdata(:,1:frames) = envelope(data(plotchans,:)); % first, plot the data envelope
fprintf('Comparing projection sizes for components: ');
compvars = zeros(1,ncomps);

for c = 1:ncomps %%% find max variances and their frame indices %%%%%

  fprintf('%d ',compnums(c)); % c is index into compnums
  if rem(c,31)==15
    fprintf('\n');
  end
  proj = icaproj(data,weights,compnums(c)); % updated arg list 12/00 -sm
  envdata(:,c*frames+1:(c+1)*frames) = envelope(proj(plotchans,:));

  [val,i] = max(sum(proj.*proj)); % find max variance
  compvars(c)   = val;

  if envdata(1,c*frames+i) > ymax % if envelop max at max variance clipped
    ix = find(envdata(1,c*frames+1:(c+1)*frames) <= ymax);
    [val,ix] = max(envdata(1,c*frames+ix));
    plotframes(c) = ix; % draw line from max non-clipped env maximum
    maxproj(:,c)  = proj(:,ix);
  else  % draw line from max envelope value at max projection time point
    plotframes(c) = i;
    maxproj(:,c)  = proj(:,i);
  end
end 
fprintf('\n');
%
%%%%%%%%%%%%%%%%%%%%%%%%% Sort by max variance in data %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
sampint = (xmax-xmin)/(frames-1);     % sampling interval = 1000/srate;
x = xmin:sampint:xmax;                % make vector of x-values

[compvars,compx] = sort(compvars');   % sort compnums on max variance
compx        = compx(ncomps:-1:1);    % reverse order of sort
compvarorder = compnums(compx);       % actual component numbers (output var)
compvars     = compvars(ncomps:-1:1)';% reverse order of sort (output var)
plotframes   = plotframes(compx);     % plotted comps have these max frames 
compframes   = plotframes';           % frame of max variance in each comp (output var)
comptimes    = x(plotframes(compx));  % time of max variance in each comp (output var)
compsplotted = compvarorder(1:ntopos);% (output var)
%
%%%%%%%%%%%%%%%%%%%%%%%% Reduce to ntopos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[plotframes,ifx] = sort(plotframes(1:ntopos));% sort plotframes on their temporal order
plottimes  = x(plotframes);           % convert to times in ms
compx      = compx(ifx);              % indices into compnums, in plotting order
maporder   = compnums(compx);         % comp. numbers, in plotting order (l->r)
maxproj    = maxproj(:,compx);        % maps in plotting order 

vlen = length(voffsets); % extend voffsets if necessary
while vlen< ntopos
  voffsets = [voffsets voffsets(vlen)]; % repeat last offset given
  vlen=vlen+1;
end

head_sep = 1.2;
topowidth = pos(3)/(ntopos+(ntopos-1)/5); % width of each topoplot
if topowidth*head_sep + pos(3) > 0.90    % adjust for maximum height
  topowidth = (0.90-0.68)/head_sep;
end
topowidth = 0.10;
if rem(ntopos,2) == 1  % odd number of topos
   topoleft = pos(3)/2 - (floor(ntopos/2)*head_sep + 0.5)*topowidth;
else % even number of topos
   topoleft = pos(3)/2 - (floor(ntopos/2)*head_sep)*topowidth;
end
if 0
 if ntopos == 3
  topoleft = 0.22;
 elseif ntopos == 2
  topoleft = 0.36;
 elseif ntopos == 1
  topoleft = 0.5; % center single topomap
 else
  topoleft = 0;
 end
end

%
%%%%%%%%%%%%%%%%%%%% Print times and frames of comp maxes %%%%%%%%%%%%%%
%

fprintf('\n');
fprintf('Plotting envelopes of %d component projections.\n',ntopos);
if length(plotchans) ~= chans
  fprintf('Envelopes computed from %d specified data channels.\n',...
      length(plotchans));
end
fprintf('Topo maps will show components: ');
for t=1:ntopos
  fprintf('%4d ',maporder(t));
end
fprintf('\n');
fprintf('    with max variance at times: ');
for t=1:ntopos
  fprintf('%4.0f ',plottimes(t));
end
fprintf('\n');
fprintf('                      = frames: ');
for t=1:ntopos
  fprintf('%4d ',plotframes(t));
end
fprintf('\n');

%
%%%%%%%%%%%%%%%%%%%%% Plot the data envelopes %%%%%%%%%%%%%%%%%%%%%%%%%
%
BACKCOLOR = [0.7 0.7 0.7];
newaxes=axes('position',pos);
axis off
%set(newaxes,'Units','Normalized','Position',...
%           [0 0 1 1],'FontSize',16,'FontWeight','Bold','Visible','off');
set(newaxes,'FontSize',16,'FontWeight','Bold','Visible','off');
set(newaxes,'Color',BACKCOLOR); % set the background color
delete(newaxes) %XXX

% site the plot at bottom of the current axes
%axe = axes('Units','Normalized','Position',...
axe = axes('Position',...
               [pos(1) pos(2) pos(3) 0.6*pos(4)],...
               'FontSize',16,'FontWeight','Bold');
limits = get(axe,'Ylim');
set(axe,'GridLineStyle',':')
set(axe,'Xgrid','off')
set(axe,'Ygrid','on')
axes(axe)
set(axe,'Color',axcolor);

fprintf('Using limits [%g,%g,%g,%g]\n',xmin,xmax,ymin,ymax);

if BOLD_COLORS==1
  mapcolors = 1:ntopos+1;
else
  mapcolors = [1 maporder+1];
end
envx = [1;compx+1];

for c = 1:ntopos+1   % plot the computed component envelopes %%%%%%%%%%%%%%%%%%

  p=plot(x,matsel(envdata,frames,0,1,envx(c)),colors(mapcolors(c),1));% plot the max
  set(gca,'FontSize',12,'FontWeight','Bold')
  if c==1                                % Note: use colors in original
    set(p,'LineWidth',2);                % component order (if BOLD_COLORS==0)
  else
    set(p,'LineWidth',1);
  end
  if mapcolors(c)>15                                % thin/dot 16th-> comp. envs.
    set(p,'LineStyle',':','LineWidth',1);
    if all_bold
      set(p,'LineStyle','-','LineWidth',3);
    end
  elseif mapcolors(c)>10                            % 
    set(p,'LineStyle',':','LineWidth',2);
    if all_bold
      set(p,'LineStyle','-','LineWidth',3);
    end
  elseif mapcolors(c)>6                             % dot 6th-> comp. envs.
    set(p,'LineStyle',':','LineWidth',3);
    if all_bold
      set(p,'LineStyle','-','LineWidth',3);
    end
  elseif mapcolors(c)>1
    set(p,'LineStyle',colors(mapcolors(c),2),'LineWidth',1);
      if colors(mapcolors(c),2) == ':'
        set(l1,'LineWidth',2);  % embolden dotted env lines
      end
  end
  hold on
  p=plot(x,matsel(envdata,frames,0,2,envx(c)),colors(mapcolors(c),1));% plot the min
  if c==1
    set(p,'LineWidth',2);
  else
    set(p,'LineWidth',1);
  end
  if mapcolors(c)>15                                % thin/dot 11th-> comp. envs.
    set(p,'LineStyle',':','LineWidth',1);
    if all_bold
      set(p,'LineStyle','-','LineWidth',3);
    end
  elseif mapcolors(c)>10                            
    set(p,'LineStyle',':','LineWidth',2);
    if all_bold
      set(p,'LineStyle','-','LineWidth',3);
    end
  elseif mapcolors(c)>6                             % dot 6th-> comp. envs.
    set(p,'LineStyle',':','LineWidth',3);
    if all_bold
      set(p,'LineStyle','-','LineWidth',3);
    end
  elseif mapcolors(c)>1
    set(p,'LineStyle',colors(mapcolors(c),2),'LineWidth',1);
      if colors(mapcolors(c),2) == ':'
        set(l1,'LineWidth',2);  % embolden dotted env lines
      end
  end
  if c==1 & ~isempty(vert)
   for v=vert
      vl=plot([v v], [ymin ymax],'k--'); % plot specified vertical lines
      set(vl,'linewidth',2.5);           % if any
   end
  end
  %
  % plot the n-th component filled 
  %
  if fill_comp_env(1)>0 & find(fill_comp_env==c-1) 
    fprintf('filling component %d\n',c);
     mins = matsel(envdata,frames,0,2,envx(c));
     p=fill([x x(frames:-1:1)],...
       [matsel(envdata,frames,0,1,envx(c)) mins(frames:-1:1)],...
       colors(mapcolors(c),1));
  end
  axis([xmin xmax ymin ymax]);
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(axe,'Color',axcolor);
l= xlabel('time (ms)');
% l= xlabel('Time (ms)');
set(l,'FontSize',14,'FontWeight','Bold');
l=ylabel('potential (uV)');
% l=ylabel('Potential (uV)');
set(l,'FontSize',14,'FontWeight','Bold');
%
%%%%%%%%%%%%%%%%%%%%%% Draw oblique/vertical lines %%%%%%%%%%%%%%%%%%%%
%
% axall = axes('Units','Normalized','Position',pos,...
axall = axes('Position',pos,...
    'Visible','Off','Fontsize',16); % whole-figure invisible axes
axes(axall)
set(axall,'Color',axcolor);
axis([0 1 0 1])

width  = xmax-xmin;
height = ymax-ymin;

for t=1:ntopos % draw oblique lines from max env vals (or plot top)
               % to map bases, in left to right order
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
  data_y = 0.6*(voffsets(t)+maxenv-ymin)/height;
  if (data_y > pos(2)+0.6*pos(4)) 
      data_y = pos(2)+0.6*pos(4);
  end
  l1 = plot([(plottimes(t)-xmin)/width  ...
               topoleft+1/pos(3)*(t-1)*6*topowidth/5+topowidth*0.6],...
                 [data_y 0.68], ...
            colors(linestyles(t)+1)); % 0.68 is bottom of topo maps
  if linestyles(t)>15                        % thin/dot 11th-> comp. envs.
    set(l1,'LineStyle',':','LineWidth',1);
    if all_bold
     set(l1,'LineStyle','-','LineWidth',3);
    end
  elseif linestyles(t)>10 
    set(l1,'LineStyle',':','LineWidth',2);
    if all_bold
     set(l1,'LineStyle','-','LineWidth',3);
    end
  elseif linestyles(t)>5                     % dot 6th-> comp. envs.
    set(l1,'LineStyle',':','LineWidth',3);
    if all_bold
     set(l1,'LineStyle','-','LineWidth',3);
    end
  elseif linestyles(t)>1
    set(l1,'LineStyle',colors(linestyles(t)+1,2),'LineWidth',1);
      if colors(linestyles(t)+1,2) == ':'
        set(l1,'LineStyle',colors(linestyles(t)+1,2),'LineWidth',2);
      end
  end
  hold on

  if voffsets(t) > 0                    % if needed add vertical lines
    l2 = plot([(plottimes(t)-xmin)/width  ...
               (plottimes(t)-xmin)/width],...
              [0.6*(maxenv-ymin)/height ...
               0.6*(voffsets(t)+maxenv-ymin)/height],...
               colors(linestyles(t)+1));
    if linestyles(t)>15                      % thin/dot 11th-> comp. envs.
      set(l2,'LineStyle',':','LineWidth',1);
      if all_bold
        set(l2,'LineStyle','-','LineWidth',3);
      end
    elseif linestyles(t)>10                   
      set(l2,'LineStyle',':','LineWidth',2);
      if all_bold
        set(l2,'LineStyle','-','LineWidth',3);
      end
    elseif linestyles(t)>5                   % dot 6th-> comp. envs.
      set(l2,'LineStyle',':','LineWidth',3);
      if all_bold
        set(l2,'LineStyle','-','LineWidth',3);
      end
    else
      set(l1,'LineStyle',colors(linestyles(t)+1,2),'LineWidth',1);
      if colors(linestyles(t)+1,2) == ':'
        set(l1,'LineWidth',2);
      end
    end
  end
  set(gca,'Visible','off');
  axis([0 1 0 1]);
end

%
%%%%%%%%%%%%%%%%%%%%%%%%% Plot the topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%
%

for t=1:ntopos % left to right order 
  % axt = axes('Units','Normalized','Position',...
  axt = axes('Units','Normalized','Position',...
       [pos(3)*topoleft+pos(1)+(t-1)*head_sep*topowidth pos(2)+0.66*pos(4) ...
        topowidth topowidth*head_sep]);
  axes(axt)                             % topoplot axes
  cla

  if ~isempty(varargin)
     topoplot(maxproj(:,t),chan_locs, varargin{:}); 
  else
     topoplot(maxproj(:,t),chan_locs,'style','both','emarkersize',3);
  end
  %ELSE headplot(winv(:,t),chan_spline);% make a 3-d headplot
                                          % if available
  axis square
  if t==1
    chid = fopen('envtopo.labels','r');
    if chid <3,
     numlabels = 1;
    else
     fprintf('Will label scalp maps with labels from pwd file %s\n','envtopo.labels');
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
  text(0.00,0.70,complabel,'FontSize',14,...
             'FontWeight','Bold','HorizontalAlignment','Center');
  % axt = axes('Units','Normalized','Position',[0 0 1 1],...
  axt = axes('Position',[0 0 1 1],...
               'Visible','Off','Fontsize',16);
  set(axt,'Color',axcolor);           % topoplot axes
  drawnow
end
axcopy(gcf);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot a colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%
%
% axt = axes('Units','Normalized','Position',[.88 .58 .03 .10]);
axt = axes('Position',[pos(1)+pos(3)*0.99 pos(2)+0.6*pos(4) pos(3)*.02 pos(4)*0.09]);
h=colorbar(axt);                        % colorbar axes
set(h,'Ytick',[]);

axes(axall)
set(axall,'Color',axcolor);
text(0.50,1.01,titl,'FontSize',16,'HorizontalAlignment','Center','FontWeight','Bold');
text(0.98,0.68,'+','FontSize',16,'HorizontalAlignment','Center');
text(0.98,0.62,'-','FontSize',16,'HorizontalAlignment','Center');

axes(axall)
set(axall,'layer','top'); % bring component lines to top

return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function envdata = envelope(data)  % also in release as env()
  if size(data,1)>1
    maxdata = max(data); % max at each time point
    mindata = min(data); % min at each time point
    envdata = [maxdata;mindata];
  else
    maxdata = max([data;data]); % max at each time point
    mindata = min([data;data]); % min at each time point
    envdata = [maxdata;mindata];
  end

return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
