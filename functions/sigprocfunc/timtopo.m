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
% Revision 1.52  2003/03/05 02:44:04  scott
% topowidth
% .,
%
% Revision 1.51  2003/03/05 02:41:49  scott
% topowidth
%
% Revision 1.50  2003/03/05 02:40:58  scott
% same
%
% Revision 1.49  2003/03/05 02:40:00  scott
% topowidth
%
% Revision 1.48  2003/03/05 02:36:31  scott
% same
%
% Revision 1.47  2003/03/05 02:36:01  scott
% eval topoplot
%
% Revision 1.46  2003/03/05 02:34:32  scott
% topoplot
%
% Revision 1.45  2003/03/05 02:33:01  scott
% topoplot
%
% Revision 1.44  2003/03/05 02:30:37  scott
% emarkersize
%
% Revision 1.43  2003/03/05 02:28:35  scott
% same
%
% Revision 1.42  2003/03/05 02:27:54  scott
% emarkersize
%
% Revision 1.41  2003/03/05 02:04:58  scott
% emarkersize
%
% Revision 1.40  2003/03/05 02:04:07  scott
% emarkersize -sm
%
% Revision 1.39  2003/03/05 01:56:34  scott
% topowidth -sm
%
% Revision 1.38  2003/03/05 01:55:20  scott
% topowidth -sm
%
% Revision 1.37  2003/03/05 01:51:57  scott
% topowidth -sm
%
% Revision 1.36  2003/03/04 21:17:52  scott
% axfont -sm
%
% Revision 1.35  2003/03/04 21:10:50  scott
% titlefont -sm
%
% Revision 1.34  2003/03/04 18:52:41  scott
% cleaning up -sm
%
% Revision 1.33  2003/03/04 18:48:25  scott
% test size -sm
%
% Revision 1.32  2003/03/04 18:44:58  scott
% title text debug -sm
%
% Revision 1.31  2003/03/04 18:43:45  scott
% title text -sm
%
% Revision 1.30  2003/03/04 18:41:09  scott
% title text -sm
%
% Revision 1.29  2003/03/04 18:40:34  scott
% title text -sm
%
% Revision 1.28  2003/03/04 18:38:57  scott
% title text
%
% Revision 1.27  2003/03/04 18:35:02  scott
% cbar text -sm
%
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

head_sep = 1.2;
pos
topowidth = pos(3)/(ntopos+(ntopos-1)/5); % width of each topoplot
topowidth
if topowidth> 0.25*pos(4) % dont make too high
  topowidth = 0.25*pos(4);
end
topowidth
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
axdata = axes('Units','Normalized','Position',[pos(1) pos(2) pos(3) 0.6*pos(4)],'FontSize',axfont);
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
set(l,'FontSize',axfont);
l=ylabel('Potential (uV)');
set(l,'FontSize',axfont);
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
             'Visible','Off','FontSize',axfont);    % whole-gca invisible axes
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

  from = changeunits([plottimes(t),maxdata],axdata,axall);
  to   = changeunits([0,-0.8],axtp,axall);
  delete(axtp);
  axes(axall);
  l1 = plot([from(1) to(1)],[from(2) to(2)]);

  hold on
  set(axall,'Visible','off');
  axis([0 1 0 1]);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%% Plot the topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%
%

for t=1:ntopos

  axtp = axes('Units','Normalized','Position',...
       [pos(3)*topoleft+pos(1)+(t-1)*head_sep*topowidth ...
              pos(2)+0.66*pos(4) ...
                  topowidth ...
                       topowidth*head_sep]);
  axes(axtp)                             % topoplot axes
  cla

  if ~isempty(varargin)
    topoargs = varargin{:};
  else
    topoargs = [];
  end
  if topowidth<0.12
    if ~isempty(topoargs)
        % topoargs = ['''emarkersize'',4,' topoargs];
        topoargs = ['''electrodes'', ''off'', ' topoargs];
    else
        % topoargs = ['''emarkersize'',4'];
        topoargs = ['''electrodes'', ''off'''];
topoargs
    end
  end
  eval(['topoplot(data(:,plotframes(t)),chan_locs,' topoargs ');']); % plot the scalp map 
  %
  % ELSE make a 3-D headplot
  %
  % headplot(data(:,plotframes(t)),'chan.spline'); 
  
  timetext = num2str(plottimes(t),'%4.0f');
  text(0.00,0.70,timetext,'FontSize',axfont-2,'HorizontalAlignment','Center');
end

%
%%%%%%%%%%%%%%%%%%%%%%%%% Make the colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%
%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot a colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%
%
axcb = axes('Position',[pos(1)+pos(3)*0.985 pos(2)+0.62*pos(4) pos(3)*0.02 pos(4)*0.09]);
h=cbar(axcb);                        % colorbar axes
pos_cb = get(axcb,'Position');
set(h,'Ytick',[]);

axes(axall)
set(axall,'Color',axcolor);
text(0.16,0.625,titl,'FontSize',titlefont,'HorizontalAlignment','Center','FontWeight','Bold');

text(0.966,0.695,'+','FontSize',axfont,'HorizontalAlignment','Center');
text(0.966,0.625,'-','FontSize',axfont,'HorizontalAlignment','Center');

axes(axall)
set(axall,'layer','top'); % bring component lines to top

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
axcopy(gcf);
