% timtopo() - plot a data epoch and map its scalp map(s) at given time(s)
%
% Usage:
%  >> timtopo(data,'chan_locs')
%  >> timtopo(data,'chan_locs',[limits],[plottimes]','title',[plotchans], ...
%                 [voffsets], 'key', 'val', ...);
%
% Inputs:
%  data       = EEG/ERP data epoch (chans,frames)
%  chan_locs  = channel location file, See  >> topoplot example 
%
% Optional inputs:
%  [limits]   = [xmin xmax ymin ymax]  (x's in ms) {def|0 or both y's 0 -> data limits}
%  plottimes  = vector of times to topoplot {default|nan -> frame of max(var())}
% 'title'     = plot title {default|0 -> none}
%  plotchans  = data channel(s) to plot {default|0 -> all}
%  voffsets   = vertical lines extend above the data this much (plot units){default -> 0}
% 'key','val' = optional topoplot() arguments. See >> help topoplot 

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

% $Log: not supported by cvs2svn $
% Revision 1.26  2003/03/04 18:34:07  scott
% cbar text
%
% Revision 1.25  2003/03/04 18:32:59  scott
% cbar text -sm
%
% Revision 1.24  2003/03/04 18:28:38  scott
% cbar text -sm
%
% Revision 1.23  2003/03/04 18:24:36  scott
% final debugs? -sm
%
% Revision 1.22  2003/03/04 18:22:30  scott
% debug last -sm
%
% Revision 1.21  2003/03/04 18:20:01  scott
% debug last -sm
%
% Revision 1.20  2003/03/04 18:15:06  scott
% debug last -sm
%
% Revision 1.19  2003/03/04 18:13:58  scott
% debug last -sm
%
% Revision 1.18  2003/03/04 18:11:36  scott
% debug last -sm
%
% Revision 1.17  2003/03/04 18:09:45  scott
% debug last -sm
%
% Revision 1.16  2003/03/04 18:07:29  scott
% debug last -sm
%
% Revision 1.15  2003/03/04 18:06:21  scott
% using changeunits -sm
%
% Revision 1.14  2003/03/04 17:49:36  scott
% debug oblique lines -sm
%
% Revision 1.13  2003/03/04 17:41:55  scott
% debug last -sm
%
% Revision 1.12  2003/03/04 17:40:51  scott
% debug head size for sbplots -sm
%
% Revision 1.11  2003/03/04 17:37:47  scott
% debug last -sm
%
% Revision 1.10  2003/03/04 17:36:26  scott
% debug last -sm
%
% Revision 1.9  2003/03/04 17:21:43  scott
% debug last -sm
%
% Revision 1.8  2003/03/04 17:20:23  scott
% debug last -sm
%
% Revision 1.7  2003/03/04 17:16:01  scott
% debug last -sm
%
% Revision 1.6  2003/03/04 17:11:51  scott
% edit to work in subplot axes -sm
%
% Revision 1.5  2002/11/15 03:07:45  arno
% header for web
%
% Revision 1.4  2002/08/28 00:52:30  arno
% allow to plot NaN with other latencies
%
% Revision 1.3  2002/08/27 00:20:54  arno
% debugging colorbar->cbar (for menus)
%
% Revision 1.2  2002/08/12 23:48:05  arno
% debug absmax
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 5-31-00 added o-time line and possibility of plotting 1 channel -sm & mw
% 11-02-99 added maplimits arg -sm
% 01-22-01 added to help message -sm
% 01-25-02 reformated help & license, added link -ad 
% 03-15-02 add all topoplot options -ad

function M = timtopo(data,chan_locs,limits,plottimes,titl,plotchans,voffsets, varargin)

if nargin < 1
   help timtopo;
   return
end

[chans,frames] = size(data);
icadefs;   

if nargin < 7 | voffsets == 0
  voffsets = zeros(1,32);
end

if nargin < 6
   plotchans = 0;
end

if plotchans==0
   plotchans = 1:chans;
end

if nargin < 5,
   titl = '';     % DEFAULT NO TITLE
end

plottimes_set=1;   % flag variable
if nargin< 4 | isempty(plottimes) | any(isnan(plottimes))
   plottimes_set = 0;
end

limitset = 0;
if nargin < 3,
    limits = 0;
elseif length(limits)>1
    limitset = 1;
end

if nargin < 2
    chan_locs = 'chan.locs';  % DEFAULT CHAN_FILE
end
if isnumeric(chan_locs) & chan_locs == 0,
    chan_locs = 'chan.locs';  % DEFAULT CHAN_FILE
end

%
%%%%%%%%%%%%%%%%%%%%%%% Read and adjust limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  if limits==0,      % == 0 or [0 0 0 0]
    xmin=0;
    xmax=frames-1;
    ymin=min(min(data));
    ymax=max(max(data));
  else
    if length(limits)~=4,
      fprintf( ...
       'timtopo: limits should be 0 or an array [xmin xmax ymin ymax].\n');
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
      fprintf('timtopo() - xmax must be > xmin.\n')
      return
  end

  if ymax == 0 & ymin == 0,
      ymax=max(max(data));
      ymin=min(min(data));
  end
  if ymax<=ymin,
      fprintf('timtopo() - ymax must be > ymin.\n')
      return
  end

sampint = (xmax-xmin)/(frames-1); % sampling interval = 1000/srate;
x = xmin:sampint:xmax;   % make vector of x-values

%
%%%%%%%%%%%%%%% Compute plotframes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if plottimes_set == 0
  [mx plotframes] = max(sum(data.*data)); 
                  % default plotting frame has max variance
  if nargin< 4 | isempty(plottimes)
	  plottimes = x(plotframes);
  else
	  plottimes(find(isnan(plottimes))) = x(plotframes);
  end;
  plottimes_set = 1;
end;

if plottimes_set == 1
	plottimes = sort(plottimes);
  ntopos = length(plottimes);
  if ntopos > 32
    fprintf('timtopo(): too many plottimes!\n');
    return
  end
  if max(plottimes) > xmax | min(plottimes)< xmin
    fprintf('timtopo(): plottimes out of range - cannot plot.\n');
    return
  end
  if sort(plottimes) ~= plottimes
    fprintf('timtopo(): plottimes out of order - lines would cross.\n');
    return
  end
  xshift = [x(2:frames) xmax];
  plotframes = ones(size(plottimes));
  for t = 1:ntopos
    time = plottimes(t);
    plotframes(t) = find(time>=x & time < xshift);
  end
end

vlen = length(voffsets); % extend voffsets if necessary
i=1;
while vlen< ntopos
        voffsets = [voffsets voffsets(i)];
        i=i+1;
        vlen=vlen+1;
end

pos = get(gca,'Position');
axis('off')
cla % clear the current axes

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

% topowidth = 0.76/(ntopos+(ntopos-1)/5); % width of each topoplot
% if topowidth*1.2 + 0.68 > 0.90    % adjust for maximum height
    % topowidth = (0.90-0.68)/1.2;
% end
% if ntopos == 2
   % topoleft = 0.32;
% elseif ntopos == 1
   % topoleft = (1-(0.94-0.68)/1.2)/2; % center single topomap
% else
   % topoleft = 0.12;
% end

if max(plotframes) > frames
    fprintf('Plot frame %d is > frames in data (%d)\n',max(plotframes),frames);
    return
end
if min(plotframes) < 1
    fprintf('Plot frame %d is < 1\n',min(plotframes));
    return
end

%
%%%%%%%%%%%%%%%%%%%% Print times and frames %%%%%%%%%%%%%%%%%%%%%%%%%%
%

fprintf('Topo maps will show times: ');
for t=1:ntopos
  fprintf('%4.0f ',plottimes(t));
end
fprintf('\n');
fprintf('                   frames: ');
for t=1:ntopos
  fprintf('%4d ',plotframes(t));
end
fprintf('\n');

%
%%%%%%%%%%%%%%%%%%%%%%% Plot the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% site the plot at bottom of the figure
axdata = axes('Units','Normalized','Position',[pos(1) pos(2) pos(3) 0.6*pos(4)],'FontSize',16);
set(axdata,'Color',BACKCOLOR);

limits = get(axdata,'Ylim');
set(axdata,'GridLineStyle',':')
set(axdata,'Xgrid','off')
set(axdata,'Ygrid','on')
axes(axdata)
axcolor = get(gcf,'Color');
set(axdata,'Color',BACKCOLOR);
pl=plot(x,data(plotchans,:));    % plot the data
if length(plotchans)==1
  set(pl,'color','k');
  set(pl,'linewidth',2);
end
l= xlabel('Time (ms)');
set(l,'FontSize',14);
l=ylabel('Potential (uV)');
set(l,'FontSize',14);
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

for t=1:ntopos % dfraw vertical lines through the data at topoplot frames
 if length(plotchans)>1 | voffsets(t)
  l1 = plot([plottimes(t) plottimes(t)],...
       [min(data(plotchans,plotframes(t))) ...
       voffsets(t) + max(data(plotchans,plotframes(t)))],'b');
 end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%% Draw oblique lines %%%%%%%%%%%%%%%%%%%%%%%%%%
%
axall = axes('Position',pos,...
             'Visible','Off','Fontsize',16);    % whole-gca invisible axes
axes(axall)
set(axall,'Color',BACKCOLOR);
axis([0 1 0 1])

  axes(axall)
  axis([0 1 0 1]);
  set(gca,'Visible','off');

for t=1:ntopos % draw oblique lines through to the topoplots 
  maxdata = max(data(:,plotframes(t))); % max data value at plotframe

  head_sep = 1.2;
  axtp = axes('Units','Normalized','Position',...
       [pos(3)*topoleft+pos(1)+(t-1)*head_sep*topowidth ...
              pos(2)+0.66*pos(4) ...
                  topowidth ...
                       topowidth*head_sep]); % this will be the topoplot axes
  postp = axis(axtp);
  axis([-1 1 -1 1]);

%   >> large = sbplot(1,1,1);   % normal whole figure axes
%   >> % Now draw line from point [0.3 0.4] in small1 axes
%   >> %                 to point [0.6 0.7] in small2 axes
%   >> from = changeunits([0.3 0.4],small1,large); % point in small1 axes
%   >> to   = changeunits([0.6 0.7],small2,large); % point in small2 axes
%   >> plot([from(1) to(1)],[from(2) to(2)])
%   >> axis([0 1 0 1]); % set large axes limits
%   >> axis off % finally, hide large axes

  from = changeunits([plottimes(t),maxdata],axdata,axall);
  to   = changeunits([0,-0.8],axtp,axall);
  delete(axtp);
  axes(axall);
  l1 = plot([from(1) to(1)],[from(2) to(2)]);

  % data_y = 0.6*(voffsets(t)+maxdata-ymin)/height;
  % if (data_y > pos(2)+0.6*pos(4))
      % data_y = pos(2)+0.6*pos(4);
  % end
  % l1 = plot([(plottimes(t)-xmin)/width  ...
                 % pos(3)*topoleft+pos(1)+(t-1)*head_sep*topowidth], ...
                 % [data_y 0.70]); % 0.70 is bottom of topo maps
% 
% % x_to was topoleft+1/pos(3)*(t-1)*6*topowidth/5+topowidth*0.6],...


%  l1 = plot(...
%     [0.12+0.76*(plottimes(t)-xmin)/width  ...
%      topoleft+(t-1)*6*topowidth/5+topowidth/2],...
%      [0.44*(voffsets(t)/height) ...
%            + 0.12+(max(data(plotchans,plotframes(t)))-ymin)/height*0.44 ...
%      0.70],'b');
  hold on
  set(axall,'Visible','off');
  axis([0 1 0 1]);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%% Plot the topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%
%

for t=1:ntopos
%  axtp = axes('Units','Normalized','Position',...
%                 [topoleft+(t-1)*6*topowidth/5 0.68 topowidth topowidth*1.2]);

  axtp = axes('Units','Normalized','Position',...
       [pos(3)*topoleft+pos(1)+(t-1)*head_sep*topowidth ...
              pos(2)+0.66*pos(4) ...
                  topowidth ...
                       topowidth*head_sep]);
  axes(axtp)                             % topoplot axes
  cla

  % topoplot(data(:,plotframes(t)),chan_locs,'style','both'); % make a topoplot
  if ~isempty(varargin)
    topoplot(data(:,plotframes(t)),chan_locs, varargin{:}); 
  else
    topoplot(data(:,plotframes(t)),chan_locs); 
  end
  %
  % ELSE make a 3-D headplot
  %
  % headplot(data(:,plotframes(t)),'chan.spline'); 
  
  timetext = num2str(plottimes(t),'%4.0f');
  text(0.00,0.70,timetext,'FontSize',14,'HorizontalAlignment','Center');
end

%
%%%%%%%%%%%%%%%%%%%%%%%%% Make the colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%
%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot a colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%
%
axcb = axes('Position',[pos(1)+pos(3)*0.99 pos(2)+0.62*pos(4) pos(3)*.02 pos(4)*0.09]);
h=cbar(axcb);                        % colorbar axes
pos_cb = get(axcb,'Position');
set(h,'Ytick',[]);

axes(axall)
set(axall,'Color',axcolor);
text(0.50,1.01,titl,'FontSize',16,'HorizontalAlignment','Center','FontWeight','Bold');

text(0.966,0.695,'+','FontSize',16,'HorizontalAlignment','Center');
text(0.966,0.625,'-','FontSize',16,'HorizontalAlignment','Center');

axes(axall)
set(axall,'layer','top'); % bring component lines to top

%  axcb = axes('Units','Normalized','Position',[.88 .58 .03 .10]); % colorbar axes
%  h=cbar(axcb);  
%  set(h,'Ytick',[]);
%
%  axes(axall)
%  set(axall,'Color',BACKCOLOR);
%  text(0.50,0.96,titl,'FontSize',16,'HorizontalAlignment','Center');
%  text(0.86,0.67,'+','FontSize',16,'HorizontalAlignment','Center');

  if ~isempty(varargin)
    try,
		if ~isempty( strmatch( 'absmax', varargin))
			text(0.86,0.624,'0','FontSize',16,'HorizontalAlignment','Center');
		end;
	catch, end;
  end
%
% Turn on axcopy()
%
axcopy(gcf);
