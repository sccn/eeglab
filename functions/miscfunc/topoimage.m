% topoimage() - plot concatenated multichannel time/frequency images 
%               in a topographic format
%               Uses a channel location file with the same format as topoplot() 
%               or else plots data on a rectangular grid of axes.
%               Click on individual images to examine separately.
%
% Usage:
%    >> topoimage(data,'chan_locs',ntimes,limits);
%    >> topoimage(data,[rows cols],ntimes,limits);
%    >> topoimage(data,'chan_locs',ntimes,limits,title,...
%                                 channels,axsize,colors,ydir,rmbase) 
%
% Inputs:
%   data       = data consisting of nchans images, each size (rows,ntimes*chans) 
%  'chan_locs' = file of channel locations as in >> topoplot example
%                Else [rows cols] matrix of locations. Example: [6 4]
%   ntimes     = columns per image 
%  [limits]    = [mintime maxtime minfreq maxfreq mincaxis maxcaxis]  
%                Give times in msec {default|0 (|both caxis 0) -> use data limits)
%  'title'     = plot title {0 -> none}
%   channels   = vector of channel numbers to plot & label {0 -> all}
%   axsize     = [x y] axis size {default [.08 .07]}
%  'colors'    = file of color codes, 3 chars per line  
%                 ( '.' = space) {0 -> default color order}
%   ydir       = y-axis polarity (pos-up = 1; neg-up = -1) {def -> pos-up}
%   rmbase     = if ~=0, remove the mean value for times<=0 for each freq {def -> no}
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 12-10-1999
%
% See also: topoplot(), timef() 

% Copyright (C) 12-10-99 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 1-16-00  debugged help msg and improved presentation -sm
% 3-16-00  added axcopy() -sm
% 8-07-00  added logimagesc() via 'LOGIT' -sm
% 9-02-00  added RMBASE option below, plus colorbar to key image -sm ???
% 1-25-02  reformated help & license, added link -ad 

function topoimage(data,loc_file,times,limits,plottitle,channels,axsize,colors,ydr,rmbas)

% Options:
% LOGIT = 1; % comment out for non-log imaging
% YVAL = 10;          % plot horizontal lines at 10 Hz (comment to omit)
% RMBASE = 0;  % remove <0 mean for each image row 

MAXCHANS = 256;
DEFAULT_AXWIDTH  = 0.08;
DEFAULT_AXHEIGHT = 0.07;
DEFAULT_SIGN = 1;   % Default - plot positive-up
LINEWIDTH = 2.0;
FONTSIZE = 14;      % font size to use for labels
CHANFONTSIZE = 10;  % font size to use for channel names
TICKFONTSIZE=10;    % font size to use for axis labels
TITLEFONTSIZE = 16;

PLOT_WIDTH = 0.75;  % width and height of plot array on figure!
PLOT_HEIGHT = 0.81;
ISRECT = 0;         % default

if nargin < 1,
    help topoimage
    return
end

if nargin < 4,
    help topoimage
    error('topoimage(): needs four arguments');
end

if times <0,
    help topoimage
    return
elseif times==1,
    fprintf('topoimage: cannot plot less than 2 times per image.\n');
    return
else
    freqs = 0;
end;

axcolor= get(0,'DefaultAxesXcolor'); % find what the default x-axis color is
plotfile = 'topoimage.ps';
ls_plotfile = 'ls -l topoimage.ps';

%
%%%%%%%%%%%%%%%%%%%%%%%%%% Substitute defaults for missing parameters %%%%%
%
SIGN = DEFAULT_SIGN;
if nargin < 10
   rmbas = 0;
end
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
  axwidth  = DEFAULT_AXWIDTH;
  axheight = DEFAULT_AXHEIGHT;
elseif size(axsize) == [1 1] & axsize(1) == 0
  axwidth  = DEFAULT_AXWIDTH;
  axheight = DEFAULT_AXHEIGHT;
elseif size(axsize) == [1 2]
  axwidth  = axsize(1);
  axheight = axsize(2);
  if axwidth > 1 | axwidth < 0 | axheight > 1 | axwidth < 0
    help topoimage
    return
  end
else
  help topoimage
  return
end

[freqs,framestotal]=size(data);             % data size
chans = framestotal/times;

fprintf('\nPlotting data using axis size [%g,%g]\n',axwidth,axheight);
if nargin < 6
   channels = 0;
end
if channels == 0
   channels = 1:chans;
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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%% Test parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  icadefs; % read MAXPLOTDATACHANS constant from icadefs.m

  if max(channels) > chans
    fprintf('topoimage(): max channel index > %d channels in data.\n',...
                       chans);
    return
  end
  if min(channels) < 1
    fprintf('topoimage(): min channel index (%g) < 1.\n',...
                       min(channels));
    return
  end;
  if length(channels)>MAXPLOTDATACHANS,
    fprintf('topoimage(): not set up to plot more than %d channels.\n',...
                       MAXPLOTDATACHANS);
    return
  end;
%
%%%%%%%%%%%%% Extend the size of the plotting area in the window %%%%%%%%%%%%
%
  curfig = gcf;
  h=figure(curfig);
  set(h,'Color',BACKCOLOR); % set the background color
  set(h,'PaperUnits','normalized'); % use percentages to avoid US/A4 difference
  set(h,'PaperPosition',[0.0235308 0.0272775 0.894169 0.909249]); % equivalent
  % orient portrait
  axis('normal');
%
%%%%%%%%%%%%%%%%%%%% Read the channel names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  if isstr(channels) == 0,
    % channames = zeros(MAXPLOTDATACHANS,4);
    % for c=1:length(channels),
    %     channames(c,:)= sprintf('%4d',channels(c));
    % end;
    channames = num2str(channels(:));                   %%CJH
  else,
    if ~isstr(channels)
       fprintf('topoimage(): channel file name must be a string.\n');
       return
    end
    chid = fopen(channels,'r');
    if chid <3,
        fprintf('topoimage(): cannot open file %s.\n',channels);
        return
    end;
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
data = matsel(data,times,0,0,channels);
chans = length(channels);
%
%%%%%%%%%%%%%%%%%%%%%%%%% Read the color names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  if colors ~=0,
    if ~isstr(colors)
       fprintf('topoimage(): color file name must be a string.\n');
       return
    end
    cid = fopen(colors,'r');
    % fprintf('cid = %d\n',cid);
    if cid <3,
        fprintf('topoimage: cannot open file %s.\n',colors);
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
    xmin=min(times);
    xmax=max(times);
    ymin=min(freqs);
    ymax=max(freqs);
  else
    if length(limits)~=6,
      fprintf( ...
       'topoimage: limits should be 0 or an array [xmin xmax ymin ymax zmin zmax].\n');
      return
    end;
    xmin = limits(1);
    xmax = limits(2);
    ymin = limits(3);
    ymax = limits(4);
    zmin = limits(5);
    zmax = limits(6);
  end;

  if xmax == 0 & xmin == 0,
    x = [0:1:times-1];
    xmin = min(x);
    xmax = max(x);
  else
    dx = (xmax-xmin)/(times-1);
    x=xmin*ones(1,times)+dx*(0:times-1); % compute x-values
    xmax = xmax*times/times;
  end;
  if xmax<=xmin,
      fprintf('topoimage() - xmax must be > xmin.\n')
      return
  end

  if ymax == 0 & ymin == 0,
      y=[1:1:freqs];
      ymax=freqs;
      ymin=1;
  else
    dy = (ymax-ymin)/(freqs-1);
    y=ymin*ones(1,freqs)+dy*(0:freqs-1); % compute y-values
    ymax = max(y);
  end
  if ymax<=ymin,
      fprintf('topoimage() - ymax must be > ymin.\n')
      return
  end

  if zmax == 0 & zmin == 0,
      zmax=max(max(data));
      zmin=min(min(data));
      fprintf('Color axis limits [%g,%g]\n',zmin,zmax);
  end
  if zmax<=zmin,
      fprintf('topoimage() - zmax must be > zmin.\n')
      return
  end

  xlabel = 'Time (ms)';
  ylabel = 'Hz';

%
%%%%%%%%%%%%%%%%%%%%%%%% Set up plotting environment %%%%%%%%%%%%%%%%%%%%%%%%%
%
  h = gcf;
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

  msg = [msg ' -> \n'];    % print starting info on screen . . .
  fprintf(...
    '\nlimits: [xmin,xmax,ymin,ymax] = [%4.1f %4.1f %4.2f %4.2f]\n',...
                xmin,xmax,ymin,ymax);

  set(h,'YLim',[ymin ymax]);            % set default plotting parameters
  set(h,'XLim',[xmin xmax]);
  set(h,'FontSize',FONTSIZE);            % choose font size

  set(h,'FontSize',FONTSIZE);            % choose font size
  set(h,'YLim',[ymin ymax]);            % set default plotting parameters
  set(h,'XLim',[xmin xmax]);

  axis('off')

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read chan_locs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if size(loc_file,2) == 2 % plot in a rectangular grid
   ISRECT = 1;
   ht = loc_file(1);
   wd = loc_file(2);
   if chans > ht*wd
      fprintf(...
         '\ntopoimage(): (d%) channels to be plotted > grid size [%d %d]\n\n',...
                           chans,ht,wd);
      return
   end
   hht = (ht-1)/2;
   hwd = (wd-1)/2;
   xvals = zeros(ht*wd,1);
   yvals = zeros(ht*wd,1);
   dist  = zeros(ht*wd,1);
   for i=1:wd
    for j=1:ht
      xvals(i+(j-1)*wd) = -hwd+(i-1);
      yvals(i+(j-1)*wd) =  hht-(j-1);
       dist(i+(j-1)*wd) =  sqrt(xvals(j+(i-1)*ht).^2+yvals(j+(i-1)*ht).^2);
    end
   end
   maxdist = max(dist);
   for i=1:wd
    for j=1:ht
      xvals(i+(j-1)*wd) = 0.499*xvals(i+(j-1)*wd)/maxdist;
      yvals(i+(j-1)*wd) = 0.499*yvals(i+(j-1)*wd)/maxdist;
    end
   end
   channames = repmat(' ',ht*wd,4);
   for i=1:ht*wd
     channum = num2str(i);
     channames(i,1:length(channum)) = channum;
   end
  
else % read chan_locs file
  fid = fopen(loc_file);
  if fid<1,
    fprintf('topoimage(): cannot open eloc_file "%s"\n',loc_file)
    return
  end
  A = fscanf(fid,'%d %f %f %s',[7 MAXCHANS]);
  fclose(fid);
  A = A';

 if length(channels) > size(A,1),
   error('topoimage(): data channels must be <= chan_locs channels')
 end

  channames = setstr(A(channels,4:7));
  idx = find(channames == '.');                       % some labels have dots
  channames(idx) = setstr(abs(' ')*ones(size(idx)));  % replace them with spaces

  Th = pi/180*A(channels,2);                          % convert degrees to rads
  Rd = A(channels,3);
  % ii = find(Rd <= 0.5); % interpolate on-head channels only
  % Th = Th(ii);
  % Rd = Rd(ii);

  [yvals,xvals] = pol2cart(Th,Rd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xvals = 0.5+PLOT_WIDTH*xvals;   % controls width of plot array on page!
yvals = 0.5+PLOT_HEIGHT*yvals;  % controls height of plot array on page!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  xdiff=xmax-xmin;
  rightmost = max(xvals);
  basetimes = find(x<=0);
P=0;
  Axes = [];
    fprintf('\ntrace %d: ',P+1);
    for I=1:chans,%%%%%%%%%% for each data channel %%%%%%%%%%%%%%%%%%%%%%%%%%
      if P>0
        axes(Axes(I))
        hold on;                      % plot down left side of page first
        axis('off')

      else % P <= 0

      %
      %%%%%%%%%%%%%%%%%%%%%%% Plot data images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
        xcenter = xvals(I);
        ycenter = yvals(I);
        Axes = [Axes axes('Units','Normal','Position', ...
              [xcenter-axwidth/2 ycenter-axheight/2 axwidth axheight])];
        axes(Axes(I))
        imageaxes = gca;
        axislcolor = get(gca,'Xcolor');   %%CJH
 
        dataimage = matsel(data,times,0,0,I);
        if rmbas~=0                  % rm baseline
          dataimage = dataimage ...
             - repmat(mean(matsel(data,times,basetimes,0,I)')',1,times); 
        end

        if exist('LOGIT')
          logimagesc(x,y,dataimage); % <---- plot logfreq image 
          if exist('YVAL')
             YVAL = log(YVAL);
          end
        else
          imagesc(x,y,dataimage);    % <---- plot image 
        end
        hold on

      curax = axis;
      xtk  = get(gca,'xtick'); % use these for cal axes below
      xtkl = get(gca,'xticklabel');
      ytk  = get(gca,'ytick');
      ytkl = get(gca,'yticklabel');

      set(gca,'tickdir','out');
      set(gca,'ticklength',[0.02 0.05]);
      set(gca,'xticklabel',[]);
      set(gca,'yticklabel',[]);
      set(gca,'ydir','normal');
      caxis([zmin zmax]);

        if exist('YVAL') & YVAL>=curax(3) & YVAL<=curax(4)
            hold on
            hp=plot([xmin xmax],[YVAL YVAL],'r-');%,'color',axislcolor);  
                                                       % draw horizontal axis 
            set(hp,'Linewidth',1.0)
        end

        if xmin<0 & xmax>0 
          hold on
          vl= plot([0 0],[curax(3) curax(4)],'color',axislcolor);  % draw vert axis 
          set(vl,'linewidth',2);
        end

       % if xcenter == rightmost
       %    colorbar
       %    rightmost = Inf;
       % end

       % secondx = 200;                              % draw second vert axis 
       % axis('off');plot([secondx secondx],[ymin ymax],'color',axislcolor); 
 
       %
       %%%%%%%%%%%%%%%%%%%%%%% Print channel names %%%%%%%%%%%%%%%%%%%%%%%%%%
       %
       NAME_OFFSET = 0.01;
       if channels~=0,                               % print channames
          if ~ISRECT    % print before topographically arrayed image 

             % axis('off');
             hold on
             h=text(xmin-NAME_OFFSET*xdiff,(curax(4)+curax(3))*0.5,[channames(I,:)]); 
             set(h,'HorizontalAlignment','right');      
             set(h,'FontSize',CHANFONTSIZE);                % choose font size

          else % print before rectangularly arrayed image 
            if xmin<0
               xmn = 0;
            else
               xmn = xmin;
            end
            % axis('off');
            h=text(xmin-NAME_OFFSET*xdiff,ymax,[channames(I,:)]); 
            set(h,'HorizontalAlignment','right');      
            set(h,'FontSize',TICKFONTSIZE);                % choose font size
          end
       end; % channels
      end; % P=0
                                                    
      % if xcenter == rightmost
      %    colorbar
      %    rightmost = Inf;
      % end

      % if xmin<0 & xmax>0
       % axes(imageaxes);
       % hold on; plot([0 0],[curax(3) curax(4)],'k','linewidth',2);
      % end

      drawnow
      fprintf(' %d',I);
    end; % %%%%%%%%%%%%%%% chan I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf('\n');

%
%%%%%%%%%%%%%%%%%%%%% Make time and freq cal axis %%%%%%%%%%%%%%%%%%%%%%%%%
%
      ax = axes('Units','Normal','Position', ...
                         [0.80 0.1 axwidth axheight]);
      axes(ax)
      axis('off');
      imagesc(x,y,zeros(size(dataimage))); hold on   % <---- plot green 
      caxis([zmin zmax]);
      set(gca,'ydir','normal');
      if xmin <=0
          py=plot([0 0],[curax(3) curax(4)],'color','k'); % draw vert axis at time zero
      else
          py=plot([xmin xmin],[curax(3) curax(4)],'color','k'); % vert axis at xmin
      end
      hold on
      if exist('YVAL') & YVAL>=curax(3)
          px=plot([xmin xmax],[YVAL YVAL],'color',axislcolor); 
                                                      % draw horiz axis at YVAL
      else
          px=plot([xmin xmax],[curax(3) curax(4)],'color',axislcolor); 
                                                      % draw horiz axis at ymin
      end
                                               
      axis(curax);
      set(gca,'xtick',xtk); % use these for cal axes 
      set(gca,'xticklabel',xtkl);
      set(gca,'ytick',ytk);
      set(gca,'yticklabel',ytkl);
      set(gca,'ticklength',[0.02 0.05]);
      set(gca,'tickdir','out');
      h = colorbar;
      cbp = get(h,'position');
      set(h,'position',[cbp(1) cbp(2) 2*cbp(3) cbp(4)]);
      caxis([zmin zmax]);
                                                    
      % secondx = 200;                    % draw second vert axis 
      % axis('off');plot([secondx secondx],[curax(3) curax(4)],'color',axislcolor); 

      %
      %%%%%%%%%%%%%%%%%%%%% Plot axis values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %

    if 0 % DETOUR
          signx = xmin-0.15*xdiff;
          axis('off');h=text(signx,SIGN*curax(3),num2str(curax(3),3));  
            set(h,'FontSize',TICKFONTSIZE);         % choose font size
            set(h,'HorizontalAlignment','right','Clipping','off');

          textx = xmin-0.6*xdiff;
          axis('off');h=text(textx,(curax(3)+curax(4))/2,ylabel);        % text Hz
            set(h,'Rotation',90);
            set(h,'FontSize',TICKFONTSIZE);         % choose font size
            set(h,'HorizontalAlignment','center','Clipping','off');

        % axis('off');h=text(signx,SIGN*ymax,['+' num2str(ymax,3)]);  % text +ymax
          axis('off');h=text(signx,SIGN*ymax,[    num2str(ymax,3)]);  % text  ymax
            set(h,'FontSize',TICKFONTSIZE);         % choose font size
            set(h,'HorizontalAlignment','right','Clipping','off');

          ytick = curax(3)-0.3*(curax(4)-curax(3));
          tick = [int2str(xmin)]; h=text(xmin,ytick,tick);            % text xmin
            set(h,'FontSize',TICKFONTSIZE);         % choose font size
            set(h,'HorizontalAlignment','center',...
                      'Clipping','off');  % center text

           h=text(xmin+xdiff/2,ytick-0.5*(curax(4)-curax(3)),xlabel);% text Times
             set(h,'FontSize',TICKFONTSIZE);         % choose font size
             set(h,'HorizontalAlignment','center',...
                       'Clipping','off');  % center text

          tick = [int2str(xmax)]; h=text(xmax,ytick,tick);            % text xmax
            set(h,'FontSize',TICKFONTSIZE);         % choose font size
            set(h,'HorizontalAlignment','center',...
                      'Clipping','off');  % center text
          axis on            
          set(ax,'xticklabel','');
          set(ax,'yticklabel','');
          set(ax,'ticklength',[0.02 0.05]);
          set(ax,'tickdir','out');
          caxis([zmin zmax]);
          hc=colorbar;
          cmapsize = size(colormap,1);
          set(hc,'ytick',[1 cmapsize]); % 

          minlabel = num2str(zmin,3);
          while (length(minlabel)<4)
             if ~contains(minlabel,'.')
                  minlabel = [minlabel '.'];
             else
                  minlabel = [minlabel '0'];
             end
          end
          maxlabel = num2str(zmax,3);
          if zmin<0 & zmax>0
             maxlabel = ['+' maxlabel];
          end
          while (length(maxlabel)<length(minlabel))
             if ~contains(maxlabel,'.')
                  maxlabel = [maxlabel '.'];
             else
                  maxlabel = [maxlabel '0'];
             end
          end
          while (length(maxlabel)>length(minlabel))
             if ~contains(minlabel,'.')
                  minlabel = [minlabel '.'];
             else
                  minlabel = [minlabel '0'];
             end
          end
          set(hc,'yticklabel',[minlabel;maxlabel]);
          set(hc,'Color',BACKCOLOR);
          set(hc,'Zcolor',BACKCOLOR);
      end % DETOUR

axcopy(gcf); % turn on pop-up axes
%
%%%%%%%%%%%%%%%%%% Make printed figure fill page %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  % orient tall
  % curfig = gcf;
  % h=figure(curfig);
  % set(h,'PaperPosition',[0.2 0.3 7.6 10]); % stretch out the plot on the page

function [returnval] = contains(strng,chr)
   returnval=0;
   for i=1:length(strng)
     if strng(i)==chr
         returnval=1;
         break
     end
   end
