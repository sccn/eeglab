% timtopo()   - plot all channels of a data epoch on the same axis 
%               and map its scalp map(s) at selected latencies.
% Usage:
%  >> timtopo(data, chan_locs);
%  >> timtopo(data, chan_locs, 'key', 'val', ...);
% Inputs:
%  data       = (channels,frames) single-epoch data matrix
%  chan_locs  = channel location file or EEG.chanlocs structure. 
%               See >> topoplot example for file format.
%
% Optional ordered inputs:
%  'limits'    = [minms maxms minval maxval] data limits for latency (in ms) and y-values
%                 (assumes uV) {default|0 -> use [0 npts-1 data_min data_max]; 
%                 else [minms maxms] or [minms maxms 0 0] -> use
%                [minms maxms data_min data_max]
%  'plottimes' = [vector] latencies (in ms) at which to plot scalp maps 
%                {default|NaN -> latency of maximum variance}
% 'title'      = [string] plot title {default|0 -> none}
% 'plotchans'  = vector of data channel(s) to plot. Note that this does not
%                affect scalp topographies {default|0 -> all}
% 'voffsets'   = vector of (plotting-unit) distances vertical lines should extend 
%                above the data (in special cases) {default -> all = standard}
%
% Optional keyword, arg pair inputs (must come after the above):
% 'topokey','val' = optional topoplot() scalp map plotting arguments. See >> help topoplot 
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 1-10-98 
%
% See also: envtopo(), topoplot()

% Copyright (C) 1-10-98 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 5-31-00 added o-time line and possibility of plotting 1 channel -sm & mw
% 11-02-99 added maplimits arg -sm
% 01-22-01 added to help message -sm
% 01-25-02 reformated help & license, added link -ad 
% 03-15-02 add all topoplot options -ad

function M = timtopo(data, chan_locs, varargin)

MAX_TOPOS = 24;

if nargin < 1 %should this be 2?
   help timtopo;
   return
end

[chans,frames] = size(data);
icadefs;   

% 
% if nargin > 2
%     if ischar(limits)
%         varargin = { limits,plottimes,titl,plotchans,voffsets varargin{:} };
%         
%     else
%         varargin = { 'limits' limits 'plottimes' plottimes 'title' titl 'plotchans' plotchans 'voffsets' voffsets varargin{:} };        
%     end;
% end;



if nargin > 2 && ~ischar(varargin{1})
   options = {};
   if length(varargin) > 0, options = { options{:} 'limits' varargin{1} }; end;
   if length(varargin) > 1, options = { options{:} 'plottimes' varargin{2} }; end;
   if length(varargin) > 2, options = { options{:} 'title'      varargin{3} }; end;
   if length(varargin) > 3, options = { options{:} 'plotchans' varargin{4} }; end;
   if length(varargin) > 4, options = { options{:} 'voffsets'     varargin{5} }; end;
   if length(varargin) > 5, options = { options{:} varargin{6:end} }; end;
else
   options = varargin;
end;

fieldlist = { 'limits'        'real'     []                       0;
              'plottimes'     'real'     []                       [];
              'title'         'string'   []                       '';
              'plotchans'     'integer'  [1:size(data,1)]         0;
              'voffsets'      'real'     []                       [] ;};
[g topoargs] = finputcheck(options, fieldlist, 'timtopo', 'ignore');

if ischar(g), error(g); end;
%Set Defaults
if isempty(g.title), g.title = ''; end;
if isempty(g.voffsets) || g.voffsets == 0, g.voffsets = zeros(1,MAX_TOPOS); end;
if isempty(g.plotchans) || isequal(g.plotchans,0), g.plotchans = 1:chans; end;
plottimes_set=1;   % flag variable
if isempty(g.plottimes) || any(isnan(g.plottimes)), plottimes_set = 0;end;
limitset = 0; %flag variable
if isempty(g.limits), g.limits = 0; end;
if length(g.limits)>1, limitset = 1; end;


% if nargin < 7 | voffsets == 0
%   voffsets = zeros(1,MAX_TOPOS);
% end
% 
% if nargin < 6
%    plotchans = 0;
% end
% 
% if plotchans==0
%    plotchans = 1:chans;
% end
% 
% if nargin < 5,
%    titl = '';     % DEFAULT NO TITLE
% end
% 
% plottimes_set=1;   % flag variable
% if nargin< 4 | isempty(plottimes) | any(isnan(plottimes))
%    plottimes_set = 0;
% end
% 
% limitset = 0;
% if nargin < 3,
%     limits = 0;
% elseif length(limits)>1
%     limitset = 1;
% end

if nargin < 2 %if first if-statement is changed to 2 should this be 3?
    chan_locs = 'chan.locs';  % DEFAULT CHAN_FILE
end
if isnumeric(chan_locs) && chan_locs == 0,
    chan_locs = 'chan.locs';  % DEFAULT CHAN_FILE
end


  %
  %%%%%%%%%%%%%%%%%%%%%%% Read and adjust limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % defaults: limits == 0 or [0 0 0 0]
  if ( length(g.limits) == 1 & g.limits==0) | (length(g.limits)==4 & ~any(g.limits))  
    xmin=0;
    xmax=frames-1;
    ymin=min(min(data));
    ymax=max(max(data));
  elseif length(g.limits) == 2  % [minms maxms] only
    ymin=min(min(data));
    ymax=max(max(data));
    xmin = g.limits(1);
    xmax = g.limits(2);
 elseif length(g.limits) == 4
    xmin = g.limits(1);
    xmax = g.limits(2);
    if any(g.limits([3 4]))
      ymin = g.limits(3);
      ymax = g.limits(4);
    else % both 0
      ymin=min(min(data));
      ymax=max(max(data));
    end
  else
    fprintf('timtopo(): limits format not correct. See >> help timtopo.\n');
    return
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
      fprintf('timtopo() - in limits, maxms must be > minms.\n')
      return
  end

  if ymax == 0 & ymin == 0,
      ymax=max(max(data));
      ymin=min(min(data));
  end
  if ymax<=ymin,
      fprintf('timtopo() - in limits, maxval must be > minmval.\n')
      return
  end

sampint = (xmax-xmin)/(frames-1); % sampling interval = 1000/srate;
x = xmin:sampint:xmax;   % make vector of x-values

%
%%%%%%%%%%%%%%% Compute plot times/frames %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if plottimes_set == 0
  [mx plotframes] = max(sum(data.*data)); 
                  % default plotting frame has max variance
  if nargin< 4 | isempty(g.plottimes)
	  g.plottimes = x(plotframes);
  else
	  g.plottimes(find(isnan(g.plottimes))) = x(plotframes);
  end;
  plottimes_set = 1;
end;

if plottimes_set == 1
  ntopos = length(g.plottimes);
  if ntopos > MAX_TOPOS
    fprintf('timtopo(): too many plottimes - only first %d will be shown!\n',MAX_TOPOS);
    g.plottimes = g.plottimes(1:MAX_TOPOS);
    ntopos = MAX_TOPOS;
  end

  if max(g.plottimes) > xmax | min(g.plottimes)< xmin
    fprintf(...
'timtopo(): at least one plottimes value outside of epoch latency range - cannot plot.\n');
    return
  end

  g.plottimes = sort(g.plottimes); % put map latencies in ascending order, 
                               % else map lines would cross.
  xshift = [x(2:frames) xmax+1]; % 5/22/2014 Ramon: '+1' was added to avoid errors when time== max(x) in line 231
  plotframes = ones(size(g.plottimes));
  for t = 1:ntopos
    time = g.plottimes(t);
    plotframes(t) = find(time>=x & time < xshift);
  end
end

vlen = length(g.voffsets); % extend voffsets if necessary
i=1;
while vlen< ntopos
        g.voffsets = [g.voffsets g.voffsets(i)];
        i=i+1;
        vlen=vlen+1;
end

%
%%%%%%%%%%%%%%%%  Compute title and axes font sizes %%%%%%%%%%%%%%%
%
pos = get(gca,'Position');
axis('off')
cla % clear the current axes
if pos(4)>0.70
   titlefont= 16;
   axfont = 16;
elseif pos(4)>0.40
   titlefont= 14;
   axfont = 14;
elseif pos(4)>0.30
   titlefont= 12;
   axfont = 12;
elseif pos(4)>0.22
   titlefont= 10;
   axfont = 10;
else
   titlefont= 8;
   axfont = 8;
end

%
%%%%%%%%%%%%%%%% Compute topoplot head width and separation %%%%%%%%%%%%%%%
%
head_sep = 0.2;
topowidth = pos(3)/((6*ntopos-1)/5); % width of each topoplot
if topowidth> 0.25*pos(4) % dont make too large (more than 1/4 of axes width)!
  topowidth = 0.25*pos(4);
end

halfn = floor(ntopos/2);
if rem(ntopos,2) == 1  % odd number of topos
   topoleft = pos(3)/2 - (ntopos/2+halfn*head_sep)*topowidth;
else % even number of topos
   topoleft = pos(3)/2 - ((halfn)+(halfn-1)*head_sep)*topowidth;
end
topoleft = topoleft - 0.01; % adjust left a bit for colorbar

if max(plotframes) > frames |  min(plotframes) < 1
    fprintf('Requested map frame %d is outside data range (1-%d)\n',max(plotframes),frames);
    return
end

%
%%%%%%%%%%%%%%%%%%%% Print times and frames %%%%%%%%%%%%%%%%%%%%%%%%%%
%

fprintf('Scalp maps will show latencies: ');
for t=1:ntopos
  fprintf('%4.0f ',g.plottimes(t));
end
fprintf('\n');
fprintf('                     at frames: ');
for t=1:ntopos
  fprintf('%4d ',plotframes(t));
end
fprintf('\n');

%
%%%%%%%%%%%%%%%%%%%%%%% Plot the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%% site the plot at bottom of the figure %%%%%%%%%%%%%%%%%%
%
axdata = axes('Units','Normalized','Position',[pos(1) pos(2) pos(3) 0.6*pos(4)],'FontSize',axfont);
set(axdata,'Color',BACKCOLOR);

g.limits = get(axdata,'Ylim');
set(axdata,'GridLineStyle',':')
set(axdata,'Xgrid','off')
set(axdata,'Ygrid','on')
axes(axdata)
axcolor = get(gcf,'Color');
set(axdata,'Color',BACKCOLOR);
pl=plot(x,data(g.plotchans,:)');    % plot the data
if length(g.plotchans)==1
  set(pl,'color','k');
  set(pl,'linewidth',2);
end
xl= xlabel('Latency (ms)');
set(xl,'FontSize',axfont);
yl=ylabel('Potential (\muV)');
set(yl,'FontSize',axfont,'FontAngle','normal');
axis([xmin xmax ymin ymax]);
hold on

%
%%%%%%%%%%%%%%%%%%%%%%%%% Plot zero time line %%%%%%%%%%%%%%%%%%%%%%%%%%
%

if xmin<0 & xmax>0
   plot([0 0],[ymin ymax],'k:','linewidth',1.5);
else
  fprintf('xmin %g and xmax %g do not cross time 0.\n',xmin,xmax)
end
%
%%%%%%%%%%%%%%%%%%%%%%%%% Draw vertical lines %%%%%%%%%%%%%%%%%%%%%%%%%%
%
width  = xmax-xmin;
height = ymax-ymin;
lwidth = 1.5;  % increment line thickness

for t=1:ntopos % dfraw vertical lines through the data at topoplot frames
 if length(g.plotchans)>1 | g.voffsets(t)
  l1 = plot([g.plottimes(t) g.plottimes(t)],...
       [min(data(g.plotchans,plotframes(t))) ...
       g.voffsets(t) + max(data(g.plotchans,plotframes(t)))],'w'); % white underline behind
  set(l1,'linewidth',2);
  l1 = plot([g.plottimes(t) g.plottimes(t)],...
       [min(data(g.plotchans,plotframes(t))) ...
       g.voffsets(t) + max(data(g.plotchans,plotframes(t)))],'b'); % blue line
  set(l1,'linewidth',lwidth);
 end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%% Draw oblique lines %%%%%%%%%%%%%%%%%%%%%%%%%%
%
axall = axes('Position',pos,...
             'Visible','Off','FontSize',axfont);   % whole-gca invisible axes
axes(axall)
set(axall,'Color',BACKCOLOR);
axis([0 1 0 1])
  axes(axall)
  axis([0 1 0 1]);
  set(gca,'Visible','off'); % make whole-figure axes invisible

for t=1:ntopos % draw oblique lines through to the topoplots 
  maxdata = max(data(:,plotframes(t))); % max data value at plotframe
  axtp = axes('Units','Normalized','Position',...
       [pos(1)+topoleft+(t-1)*(1+head_sep)*topowidth ...
              pos(2)+0.66*pos(4) ...
                  topowidth ...
                       topowidth*(1+head_sep)]); % this will be the topoplot axes
                       % topowidth]); % this will be the topoplot axes
  axis([-1 1 -1 1]);

  from = changeunits([g.plottimes(t),maxdata],axdata,axall); % data axes
  to   = changeunits([0,-0.74],axtp,axall);                % topoplot axes
  delete(axtp);
  axes(axall);                                             % whole figure axes
  l1 = plot([from(1) to(1)],[from(2) to(2)]);
  set(l1,'linewidth',lwidth);

  hold on
  set(axall,'Visible','off');
  axis([0 1 0 1]);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%% Plot the topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%
%
topoaxes = zeros(1,ntopos);
for t=1:ntopos
       % [pos(3)*topoleft+pos(1)+(t-1)*(1+head_sep)*topowidth ...
  axtp = axes('Units','Normalized','Position',...
       [topoleft+pos(1)+(t-1)*(1+head_sep)*topowidth ...
              pos(2)+0.66*pos(4) ...
                  topowidth topowidth*(1+head_sep)]);
  axes(axtp)                             % topoplot axes
  topoaxes(t) = axtp; % save axes handles
  cla

  if topowidth<0.12
      topoargs = { topoargs{:} 'electrodes' 'off' };
  end
  topoplot( data(:,plotframes(t)),chan_locs, topoargs{:});

  % Else make a 3-D headplot
  %
  % headplot(data(:,plotframes(t)),'chan.spline'); 
  
  timetext = [num2str(g.plottimes(t),'%4.0f')];
  % timetext = [num2str(plottimes(t),'%4.0f') ' ms']; % add ' ms'
  text(0.00,0.80,timetext,'FontSize',axfont-3,'HorizontalAlignment','Center'); % ,'fontweight','bold');
end

%
%%%%%%%%%%%%%%%%%%% Plot a topoplot colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
axcb = axes('Position',[pos(1)+pos(3)*0.995 pos(2)+0.62*pos(4) pos(3)*0.02 pos(4)*0.09]);
h=cbar(axcb);                        % colorbar axes
pos_cb = get(axcb,'Position');
set(h,'Ytick',[]);

axes(axall)
set(axall,'Color',axcolor);
%
%%%%%%%%%%%%%%%%%%%%% Plot the color bar '+' and '-' %%%%%%%%%%%%%%%%%%%%%%%%%%
%
text(0.986,0.695,'+','FontSize',axfont,'HorizontalAlignment','Center');
text(0.986,0.625,'-','FontSize',axfont,'HorizontalAlignment','Center');

%
%%%%%%%%%%%%%%%%%%%%%%%%% Plot the plot title if any %%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot title between data panel and topoplots (to avoid crowding at top of
% figure), on the left
ttl = text(0.03,0.635,g.title,'FontSize',titlefont,'HorizontalAlignment','left'); % 'FontWeight','Bold');

% textent = get(ttl,'extent');
% titlwidth = textent(3);
% ttlpos = get(ttl,'position');
% set(ttl,'position',[     ttlpos(2), ttlpos(3)]);

axes(axall)
set(axall,'layer','top'); % bring component lines to top
for t = 1:ntopos
  set(topoaxes(t),'layer','top'); % bring topoplots to very top
end

  if ~isempty(varargin)
    try,
		if ~isempty( strmatch( 'absmax', varargin))
			text(0.86,0.624,'0','FontSize',axfont,'HorizontalAlignment','Center');
		end;
	catch, end;
  end

%
% Turn on axcopy()
%

% clicking on ERP pop_up topoplot
% -------------------------------
disp('Click on ERP waveform to show scalp map at specific latency');

dat.times = x;
dat.erp   = data;
dat.chanlocs = chan_locs;
dat.options  = topoargs;
dat.srate    = (size(data,2)-1)/(x(end)-x(1))*1000;
dat.axes     = axtp;
dat.line     = l1;

cb_code = [ 'tmppos = get(gca, ''currentpoint'');' ...
            'dattmp = get(gcf, ''userdata'');' ...
            'set(dattmp.line, ''visible'', ''off'');' ...
            'axes(dattmp.axes); cla;' ...
            'latpoint = round((tmppos(1)-dattmp.times(1))/1000*dattmp.srate);' ...
            'topoplot(dattmp.erp(:,latpoint), dattmp.chanlocs, dattmp.options{:});' ...
            'title(sprintf(''%.0f ms'', tmppos(1)));' ...
            'clear latpoint dattmp;' ...
          ];
axcopy;

set(gcf, 'userdata', dat);
set(axdata, 'ButtonDownFcn', cb_code); %windowbuttondownfcn', cb_code);
set(pl, 'ButtonDownFcn', cb_code);

%axcopy(gcf, cb_code);
