%  plottopo() - plot concatenated multichannel data epochs in a topographic or
%               rectangular array. Uses a channel location file with the same 
%               format as topoplot(), or else plots data on a rectangular grid. 
%               If data are all positive, they are assumed to be spectra.
% Usage:
%    >> plottopo(data,'chan_locs')
%    >> plottopo(data,[rows cols])
%    >> plottopo(data,'chan_locs',frames,limits,title,channels,...
%                                            axsize,colors,ydir,vert) 
%
% Inputs:
%   data       = data consisting of consecutive epochs of (chans,frames) 
%  'chan_locs' = file of channel locations as in >> topoplot example   {grid}
%                ELSE: [rows cols] grid size for rectangular matrix. Example: [6 4]
%   frames     = time frames (points) per epoch {def|0 -> data length}
%  [limits]    = [xmin xmax ymin ymax]  (x's in ms or Hz) {def|0 
%                 (or both y's 0) -> use data limits)
%  'title'     = plot title {def|0 -> none}
%   channels   = vector of channel numbers to plot & label {def|0 -> all}
%                   else, filename of ascii channel-name file
%   axsize     = [x y] axis size {default [.07 .07]}
%  'colors'    = file of color codes, 3 chars per line  
%                ( '.' = space) {0 -> default color order}
%   ydir       = y-axis polarity (pos-up = 1; neg-up = -1) {def -> pos-up}
%   vert       = [vector] of times (in ms or Hz) to plot vertical lines {def none}
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 3-2-98 
%
% See also: plotdata(), topoplot()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 3-2-98 from plotdata() Scott Makeig, SCCN/INC/UCSD,
% scott@sccn.ucsd.edu
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
% Revision 1.4  2002/04/24 18:21:35  scott
% added vertcolor -sm
%
% Revision 1.3  2002/04/24 18:19:52  scott
% [same] -sm
%
% Revision 1.2  2002/04/24 17:53:46  scott
% added 'vert' vertical line plotting -sm
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 5-11-98 added channels arg -sm
% 7-15-98 added ydir arg, made pos-up the default -sm
% 7-23-98 debugged ydir arg and pos-up default -sm
% 12-22-99 added grid size option, changed to sbplot() order -sm
% 03-16-00 added axcopy() feature -sm & tpj
% 08-21-00 debugged axheight/axwidth setting -sm
% 01-25-02 reformated help & license, added links -ad 
% 03-11-02 change the channel names ploting position and cutomize pop-up -ad 
% 03-15-02 added readlocs and the use of eloc input structure -ad 
% 03-15-02 debuging chanlocs structure -ad & sm 

function plottopo(data,chan_locs,frames,limits,plottitle,channels,axsize,colors,ydr,vert)

if nargin < 1 | nargin > 10
    help plottopo
    return
end

%
%%%%%%%%%%%%% Extend the size of the plotting area in the window %%%%%%%%%%%%
%
  curfig = gcf;
  h=figure(curfig);
  set(h,'PaperUnits','normalized'); % use percentages to avoid US/A4 difference
  set(h,'PaperPosition',[0.0235308 0.0272775 0.894169 0.909249]); % equivalent
  orient portrait
  axis('normal');
%
%%%%%%%%%%%%%%%%%%%%% Graphics Settings - can be customized %%%%%%%%%%%%%%%%%%
%
LINEWIDTH     = 2.0;     % data line widths (can be non-integer)
FONTSIZE      = 14;      % font size to use for labels
CHANFONTSIZE  = 12;      % font size to use for channel names
TICKFONTSIZE  = 10;      % font size to use for axis labels
TITLEFONTSIZE = 16;      % font size to use for the plot title
PLOT_WIDTH    = 0.75;    % width and height of plot array on figure
PLOT_HEIGHT   = 0.81;
gcapos = get(gca,'Position');
PLOT_WIDTH    = gcapos(3)*PLOT_WIDTH; % width and height of gca plot array on gca
PLOT_HEIGHT   = gcapos(4)*PLOT_HEIGHT;
MAXCHANS      = 256;     % can be increased
%
%%%%%%%%%%%%%%%%%%%% Default settings - use commandline to override %%%%%%%%%%%
%
DEFAULT_AXWIDTH  = 0.07;
DEFAULT_AXHEIGHT = 0.07;
DEFAULT_SIGN = 1;                         % Default - plot positive-up
ISRECT = 0;                               % default
ISSPEC = 0;                               % Default - not spectral data 

icadefs; % read BACKCOLOR, MAXPLOTDATACHANS constant from icadefs.m
set(gca,'Color',BACKCOLOR);               % set the background color

[chans,framestotal]=size(data);           % data size

axcolor= get(0,'DefaultAxesXcolor'); % find what the default x-axis color is
vertcolor = 'b';
plotfile = 'plottopo.ps';
ls_plotfile = 'ls -l plottopo.ps';
%
%%%%%%%%%%%%%%% Substitute defaults for missing parameters %%%%%%%%%%%%%%%%
%

if nargin < 10
   vert = [];
end
if vert==0
   vert = [];
end
SIGN = DEFAULT_SIGN;
if nargin < 9
   ydr = 0;
end
if ydr == -1
   SIGN = -1;
end
  
if nargin < 8
    colors = 0;
end

if nargin < 7,
  axwidth  = nan; % DEFAULT_AXWIDTH;
  axheight = nan; % DEFAULT_AXHEIGHT;
elseif size(axsize) == [1 1] & axsize(1) == 0
  axwidth  = nan; % DEFAULT_AXWIDTH;
  axheight = nan; % DEFAULT_AXHEIGHT;
elseif size(axsize) == [1 2]
  axwidth  = axsize(1);
  axheight = axsize(2);
  if axwidth > 1 | axwidth < 0 | axheight > 1 | axwidth < 0
    help plottopo
    return
  end
else
  help plottopo
  return
end
if nargin < 6
   channels = 0;
end
if channels == 0
   channelnos = 1:size(data,1);
elseif ~isstr(channels)
   channelnos = channels;
else
   channelnos = 1:size(data,1);
end
if nargin < 5
    plottitle = 0; %CJH
end
limitset = 0;
if nargin < 4,
    limits = 0;
elseif length(limits)>1
    limitset = 1;
end
if nargin < 3,
    frames = 0;
end
if nargin < 2
  chan_locs = '';
end
if (isnumeric(chan_locs) & (chan_locs == 0)) | isempty(chan_locs) 
  chan_locs = '';
end
if isempty(chan_locs) 
  n = ceil(sqrt(chans));
  chan_locs = [n n];
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Test parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  if frames <=0,
    frames = framestotal;    % default
    datasets=1;
  elseif frames==1,
    fprintf('plottopo: cannot plot less than 2 frames per trace.\n');
    return
    datasets=1;
  else
    datasets = fix(framestotal/frames);        % number of traces to overplot
  end;

  if max(channelnos) > chans
  channelnos
    fprintf('plottopo(): max channel index > %d channels in data.\n',...
                       chans);
    return
  end
  if min(channelnos) < 1
    fprintf('plottopo(): min channel index (%g) < 1.\n',...
                       min(channels));
    return
  end;
  if length(channelnos)>MAXPLOTDATACHANS,
    fprintf('plottopo(): not set up to plot more than %d channels.\n',...
                       MAXPLOTDATACHANS);
    return
  end;

  if datasets>MAXPLOTDATAEPOCHS 
      fprintf('plottopo: not set up to plot more than %d epochs.\n',...
                       MAXPLOTDATAEPOCHS);
    return
  end;
  if datasets<1
      fprintf('plottopo: cannot plot less than 1 epoch!\n');
      return
  end;

  if size(chan_locs,2)==2 % if grid plot
    if isnan(axheight) % if not specified
      axheight = gcapos(4)/(chan_locs(1)+1);
      axwidth  = gcapos(3)/(chan_locs(2)+1);
    end
    % if chan_locs(2) > 5
     %     axwidth = 0.66/(chan_locs(2)+1);
    % end
  else
     axheight = DEFAULT_AXHEIGHT;
     axwidth =  DEFAULT_AXWIDTH;
  end
    fprintf('\nPlotting data using axis size [%g,%g]\n',axwidth,axheight);
%
%%%%%%%%%%%%%%%%%%%% Read the channel names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  if ~isstr(channels) 
    % channames = zeros(MAXPLOTDATACHANS,4);
    % for c=1:length(channels),
    %     channames(c,:)= sprintf('%4d',channels(c));
    % end;
    channames = num2str(channels(:));                   %%CJH
  else % isstr(channels)
    if ~isstr(channels)
       fprintf('plottopo(): channel file name must be a string.\n');
       return
    end
    chid = fopen(channels,'r');
    if chid <3,
        fprintf('plottopo(): cannot open file %s.\n',channels);
        return
    else
        fprintf('plottopo(): opened file %s.\n',channels);
    end;

    %%%%%%%
    % fid=fopen('fgetl.m');
    % while 1
    %   line = fgetl(fid);
    %   if ~isstr(line), break, end
    %     disp(line)
    %   end
    % end
    % fclose(fid);
    %%%%%%%                       

    channames = fscanf(chid,'%s',[4 MAXPLOTDATACHANS]);
    channames = channames';
       [r c] = size(channames);
    for i=1:r
        for j=1:c
            if channames(i,j)=='.',
                channames(i,j)=' ';
            end;
        end;
    end;
  end; % setting channames
%
%%%%%%%%%%%%%%%%%%%%%%%%% Plot and label specified channels %%%%%%%%%%%%%%%%%%
%
data = data(channelnos,:);
chans = length(channelnos);
%
%%%%%%%%%%%%%%%%%%%%%%%%% Read the color names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  if colors ~=0,
    if ~isstr(colors)
       fprintf('plottopo(): color file name must be a string.\n');
       return
    end
    cid = fopen(colors,'r');
    % fprintf('cid = %d\n',cid);
    if cid <3,
        fprintf('plottopo: cannot open file %s.\n',colors);
        return
    end;
    colors = fscanf(cid,'%s',[3 MAXPLOTDATAEPOCHS]);
    colors = colors';
       [r c] = size(colors);
    for i=1:r
        for j=1:c
            if colors(i,j)=='.',
                colors(i,j)=' ';
            end;
        end;
    end;
  else % use default color order (no yellow!)
     colors =['r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  '];
     colors = [colors; colors];  % make > 64 available
  end;
  for c=1:length(colors)   % make white traces black unless axis color is white
    if colors(c,1)=='w' & axcolor~=[1 1 1]
         colors(c,1)='k';
    end
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
       'plottopo: limits should be 0 or an array [xmin xmax ymin ymax].\n');
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
      fprintf('plottopo() - xmax must be > xmin.\n')
      return
  end

  if ymax == 0 & ymin == 0,
      ymax=max(max(data));
      ymin=min(min(data));
  end
  if ymax<=ymin,
      fprintf('plottopo() - ymax must be > ymin.\n')
      return
  end

  xlabel = 'Time (ms)';
  if ISSPEC
    ISSPEC = 1;
    SIGN = 1;
    fprintf('\nPlotting positive up. Assuming data are spectra.\n');
    xlabel = 'Freq (Hz)';
    ymin = 0;                        % plot positive-up
  end;
%
%%%%%%%%%%%%%%%%%%%%%% Set up plotting environment %%%%%%%%%%%%%%%%%%%%%%%%%
%
  % h = gcf;
  % set(h,'YLim',[ymin ymax]);       % set default plotting parameters
  % set(h,'XLim',[xmin xmax]);
  % set(h,'FontSize',18);
  % set(h,'DefaultLineLineWidth',1); % for thinner postscript lines
%
%%%%%%%%%%%%%%%%%%%%%%%%%% Print plot info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  % clf;   % clear the current figure

  % print plottitle over (left) subplot 1
  if plottitle==0,
    plottitle = '';
  end
  h=gca;title(plottitle,'FontSize',TITLEFONTSIZE); % title plot and
  hold on
  msg = ['\nPlotting %d traces of %d frames with colors: '];

  for c=1:datasets
    msg = [msg  colors(c,:)];
  end
  msg = [msg ' -> \n'];    % print starting info on screen . . .
  fprintf(...
    '\nlimits: [xmin,xmax,ymin,ymax] = [%4.1f %4.1f %4.2f %4.2f]\n',...
                xmin,xmax,ymin,ymax);
  fprintf(msg,datasets,frames);

  set(h,'YLim',[ymin ymax]);            % set default plotting parameters
  set(h,'XLim',[xmin xmax]);
  set(h,'FontSize',FONTSIZE);           % choose font size

  set(h,'FontSize',FONTSIZE);           % choose font size
  set(h,'YLim',[ymin ymax]);            % set default plotting parameters
  set(h,'XLim',[xmin xmax]);

  axis('off')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Read chan_locs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if size(chan_locs,2)==2 % plot in a rectangular grid
   ISRECT = 1;
   ht = chan_locs(1);
   wd = chan_locs(2);
   if chans > ht*wd
      fprintf(...
    '\ntopoplot(): (%d) channels to be plotted > grid size [%d %d]\n',...
                           chans,ht,wd);
      return
   end
   halfht = (ht-1)/2;
   halfwd = (wd-1)/2;
   xvals = zeros(ht*wd,1);
   yvals = zeros(ht*wd,1);
   dist  = zeros(ht*wd,1);
   for j=1:ht  
     for i=1:wd
       xvals(i+(j-1)*wd) = -halfwd+(i-1);
       yvals(i+(j-1)*wd) =  halfht-(j-1);
     % dist(i+(j-1)*wd) =  sqrt(xvals(j+(i-1)*ht).^2+yvals(j+(i-1)*ht).^2);
     end
   end
   % maxdist = max(dist);
   maxxvals = max(xvals);
   maxyvals = max(yvals);
   for j=1:ht
     for i=1:wd
       % xvals(i+(j-1)*wd) = 0.499*xvals(i+(j-1)*wd)/maxdist; 
       % yvals(i+(j-1)*wd) = 0.499*yvals(i+(j-1)*wd)/maxdist; 
         xvals(i+(j-1)*wd) = 0.499*xvals(i+(j-1)*wd)/maxxvals; 
         yvals(i+(j-1)*wd) = 0.499*yvals(i+(j-1)*wd)/maxyvals; 
     end
   end
   if ~exist('channames')
     channames = repmat(' ',ht*wd,4);
     for i=1:ht*wd
       channum = num2str(i);
       channames(i,1:length(channum)) = channum;
     end
   end
  
else % read chan_locs file
	% read the channel location file
	% ------------------------------
	[tmp channames Th Rd] = readlocs(chan_locs);
	channames = strvcat(channames{channelnos});
	Th = pi/180*Th(channelnos);                 % convert degrees to radians
	Rd = Rd(channelnos); 
	
%  fid = fopen(chan_locs);
%  if fid<1,
%    fprintf('plottopo(): cannot open chan_locs file "%s"\n',chan_locs)
%    return
%  end
%  A = fscanf(fid,'%d %f %f %s',[7 MAXCHANS]);
%  fclose(fid);
%  A = A';
  %channames = setstr(A(channelnos,4:7));
  %idx = find(channames == '.');                     % some labels have dots
  %channames(idx) = setstr(abs(' ')*ones(size(idx)));% replace them with spaces
  %Th = pi/180*A(channelnos,2);                      % convert degrees to rads
  %Rd = A(channelnos,3);
  % ii = find(Rd <= 0.5); % interpolate on-head channels only
  % Th = Th(ii);
  % Rd = Rd(ii);

  if length(channelnos) > length(Th),
    error('plottopo(): data channels must be <= chan_locs channels')
  end

  [yvals,xvals] = pol2cart(Th,Rd); % translate from polar to cart. coordinates
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xvals = 0.5+PLOT_WIDTH*xvals;   % controls width of  plot array on page!
% yvals = 0.5+PLOT_HEIGHT*yvals;  % controls height of plot array on page!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xvals = gcapos(1)+gcapos(3)/2+PLOT_WIDTH*xvals;   % controls width of plot 
                                                  % array on current axes
yvals = gcapos(2)+gcapos(4)/2+PLOT_HEIGHT*yvals;  % controls height of plot 
                                                  % array on current axes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

  xdiff=xmax-xmin;
  ydiff=ymax-ymin;

  Axes = [];
  for P=0:datasets-1, %  for each data epoch
      fprintf('\ntrace %d: ',P+1);

    for c=1:chans, %%%%%%%% for each data channel %%%%%%%%%%%%%%%%%%%%%%%%%%

      if P>0 % subsequent pages (Axes specified)
        axes(Axes(c))
        hold on;                      % plot down left side of page first
        axis('off')
      else   % first page, specify axes
        xcenter = xvals(c);
        ycenter = yvals(c);
        Axes = [Axes axes('Units','Normal','Position', ...
              [xcenter-axwidth/2 ycenter-axheight/2 axwidth axheight])];
        axes(Axes(c))
        axis('off')
     
        hold on;                      % plot down left side of page first
        % set(h,'YLim',[ymin ymax]);    % set default plotting parameters
        % set(h,'XLim',[xmin xmax]);

        axislcolor = get(gca,'Xcolor');   %%CJH

        axis('off');
        if ISSPEC
          plot([xmin xmin],[0 ymax],'color',axislcolor); 
        else
          plot([0 0],[ymin ymax],'color',axislcolor); % draw vert axis at time 0  
        end  
        axis('off');
        plot([xmin xmax],[0 0],'color',axislcolor);  % draw horizontal axis 
                                                
        % secondx = 200;                             % draw second vert axis 
        % axis('off');plot([secondx secondx],[ymin ymax],'color',axislcolor); 
       %
       %%%%%%%%%%%%%%%%%%%%%%% Print channel names %%%%%%%%%%%%%%%%%%%%%%%%%%
       %
       NAME_OFFSET = -1;
       NAME_OFFSETY = -0.5;
       if channels~=0,                               % print channames
        if ISSPEC
          axis('off'),h=text(xmin-NAME_OFFSET*xdiff,ymax/2,[channames(c,:)]); 
            set(h,'HorizontalAlignment','right');    % print before traces
            set(h,'FontSize',CHANFONTSIZE);              % choose font size
        else % ~ISSPEC
          if ymin <= 0 & ymax >= 0,
            yht = 0;
          else
            yht = mean(SIGN*data(c,1+P*frames:1+P*frames+frames-1));
          end
          if ~ISRECT    % print before traces
             axis('off'),h=text(xmin-NAME_OFFSET*xdiff,yht-NAME_OFFSETY*ydiff,[channames(c,:)]); 
             set(h,'HorizontalAlignment','right');      
             set(h,'FontSize',CHANFONTSIZE);           % choose font size
          else % ISRECT
            if xmin<0
               xmn = 0;
            else
               xmn = xmin;
            end
            axis('off'),h=text(xmn,ymax,[channames(c,:)]); 
            set(h,'HorizontalAlignment','right');      
            set(h,'FontSize',TICKFONTSIZE);            % choose font size
          end % ISRECT
        end % ~ISSPEC
       end; % channels~=0
      end; % P=0
      %
      %%%%%%%%%%%%%%%%%%%%%%% Plot data traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      if ~ISSPEC % -/+ plot, normal case (e.g., not spectra), plot data trace
                                                    
        plot(x,SIGN*data(c,1+P*frames:1+P*frames+frames-1),colors(P+1),...
                       'linewidth',LINEWIDTH);   
        ymn = min(SIGN*[ymax ymin]);
        ymx = max(SIGN*[ymax ymin]);
        axis([xmin xmax ymn ymx]);          % set axis bounds
      else % ISSPEC
        plot(x,SIGN*data(c,1+P*frames:1+P*frames+frames-1),colors(P+1),...
                       'linewidth',LINEWIDTH);   
        ymaxm = ymax;
        if ymaxm/2. > ymax,
            ymaxm = ymaxm/2.;
        end;
        axis([xmin xmax ymin ymaxm]);      % set axis values
      end
      %
      %%%%%%%%%%%%%%%%%%%% plot vertical lines (optional) %%%%%%%%%%%%%%%%%
      %
      if ~isnan(vert)
       if ~ISSPEC % -/+ plot, normal case (e.g., not spectra), plot data trace
        ymean = (ymin+ymax)/2; 
        vmin = ymean-0.5*(ymean-ymin);
        vmax = ymean+0.5*(ymax-ymean);
        for v = vert
          plot([v v],[vmin vmax],'color',vertcolor); % draw vertical lines 
        end
       else
        for v = vert
          plot([v v],[0 ymax],'color',vertcolor); 
        end
       end
      end

     fprintf(' %d',c); % finished with channel plot
    end; % c, chans / subplot
  end; % P / epoch
  fprintf('\n');
  %
  %%%%%%%%%%%%%%%%%%%%% Make time and amp cal bar %%%%%%%%%%%%%%%%%%%%%%%%%
  %
  ax = axes('Units','Normal','Position', ...
                         [0.80 0.1 axwidth axheight]); % FIX!!!!
  axes(ax)
  axis('off');
  if ~ISSPEC,
    if xmin <=0
      p=plot([0 0],[ymn ymx],'color','k'); % draw vert axis at zero
    else
      p=plot([xmin xmin],[ymn ymx],'color','k'); % draw vert axis at zero
    end
    axis([xmin xmax ymn ymx]);        % set axis values
    hold on
    %set(p, 'Clipping','off');        % center text
  elseif ISSPEC
    ylo=0;
    plot([xmin xmin],[0 ymax],'color',axislcolor); 
    axis([xmin xmax ylo ymaxm]);      % set axis values
  end  
  p=plot([xmin xmax],[0 0],'color',axislcolor); % draw horizontal axis 
  axis([xmin xmax ymin ymax]);        % set axis values
  %
  %%%%%%%%%%%%%%%%%%%% plot vertical lines (optional) %%%%%%%%%%%%%%%%%
  %
  if ~isnan(vert)
   if ~ISSPEC % -/+ plot, normal case (e.g., not spectra), plot data trace
    for v = vert
      plot([v v],[vmin vmax],'color',vertcolor); % draw vertical lines 
    end
   else
    for v = vert
      plot([v v],[0 ymax],'color',vertcolor); 
    end
   end
  end
                                               
  % secondx = 200;                    % draw second vert axis 
  % axis('off');plot([secondx secondx],[ylo ymax],'color',axislcolor); 
  %
  %%%%%%%%%%%%%%%%%%%%% Plot negative-up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  if ~ISSPEC % not spectral data
                                                    
    signx = xmin-0.15*xdiff;
    axis('off');h=text(signx,SIGN*ymin,num2str(ymin,3)); % text ymin
    set(h,'FontSize',TICKFONTSIZE);               % choose font size
    set(h,'HorizontalAlignment','right','Clipping','off');

    axis('off');h=text(signx,SIGN*ymax,['+' num2str(ymax,3)]);  % text +ymax
    set(h,'FontSize',TICKFONTSIZE);         % choose font size
    set(h,'HorizontalAlignment','right','Clipping','off');

    ytick = -ymax-0.3*ydiff;
    tick = [int2str(xmin)]; h=text(xmin,ytick,tick);
    set(h,'FontSize',TICKFONTSIZE);         % choose font size
    set(h,'HorizontalAlignment','center',...
                        'Clipping','off');  % center text

    tick = [xlabel]; h=text(xmin+xdiff/2,ytick-0.5*ydiff,tick);
    set(h,'FontSize',TICKFONTSIZE);         % choose font size
    set(h,'HorizontalAlignment','center',...
                        'Clipping','off');  % center text

    tick = [int2str(xmax)]; h=text(xmax,ytick,tick);
    set(h,'FontSize',TICKFONTSIZE);         % choose font size
    set(h,'HorizontalAlignment','center',...
                        'Clipping','off');  % center text
    %
    %%%%%%%%%%%%%%%%%%%%% Plot positive-up [0,ymax] %%%%%%%%%%%%%%%%%%%%%%%%
    %
    else % ISSPEC
      ymin=0;
      signx = xmin-0.15*xdiff;

      axis('on');h=text(signx,-1*ymin,num2str(ymin,3));% text ymin
      set(h,'FontSize',TICKFONTSIZE);           % choose font size
      set(h,'HorizontalAlignment','right','Clipping','off');

      axis('on');h=text(signx,-1*ymax,['+' num2str(ymax,3)]); % text +ymax
      set(h,'FontSize',TICKFONTSIZE);           % choose font size
      set(h,'HorizontalAlignment','right','Clipping','off');

      ytick = -ymax-0.25*ydiff;

      tick = [int2str(xmin)]; h=text(xmin,ytick,tick);
      set(h,'FontSize',TICKFONTSIZE);         % choose font size
      set(h,'HorizontalAlignment','center',...
                          'Clipping','off');  % center text

      tick = [xlabel]; h=text(xmin+xdiff/2,ytick,tick);
      set(h,'FontSize',TICKFONTSIZE);         % choose font size
      set(h,'HorizontalAlignment','center',...
                          'Clipping','off');  % center text

      tick = [int2str(xmax)]; h=text(xmax,ytick,tick);
      set(h,'FontSize',TICKFONTSIZE);         % choose font size
      set(h,'HorizontalAlignment','center',...
                          'Clipping','off');  % center text
    end; % if ISSPEC

    axcopy(gcf, 'axis on; clear xlabel ylabel; xlabel(''''Time (ms)''''); ylabel(''''Voltage (\muV)'''');'); % turn on popup feature
%
%%%%%%%%%%%%%%%%%% Make printed figure fill page %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
orient tall
  % curfig = gcf;
  % h=figure(curfig);
  % set(h,'PaperPosition',[0.2 0.3 7.6 10]); % stretch out the plot on the page

