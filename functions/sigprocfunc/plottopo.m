%  plottopo() - plot concatenated multichannel data epochs in a topographic or
%               rectangular array. Uses a channel location file with the same 
%               format as topoplot(), or else plots data on a rectangular grid. 
%               If data are all positive, they are assumed to be spectra.
% Usage:
%    >> plottopo(data,'chan_locs')
%    >> plottopo(data,[rows cols])
%    >> plottopo(data,'chan_locs',frames,limits,title,channels,...
%                                            axsize,colors,ydir,vert) 
% Inputs:
%   data       = data consisting of consecutive epochs of (chans,frames)
%                or (chans,frames,n)
%
% Optional inputs:
%  'chanlocs'  = [struct] channel structure or file plot ERPs at channel 
%                locations. See help readlocs() for data channel format.
%  'geom'      = [rows cols] plot ERP in grid (overwrite previous option).
%                Grid size for rectangular matrix. Example: [6 4].
%  'frames'    = time frames (points) per epoch {def|0 -> data length}
%  'limits'    = [xmin xmax ymin ymax]  (x's in ms or Hz) {def|0 
%                 (or both y's 0) -> use data limits)
%  'title'     = [string] plot title {def|'' -> none}
%  'chan'      = vector of channel numbers to plot {def|0 -> all}
%  'axsize'    = [x y] axis size {default [.07 .07]}
%  'legend'    = [cell array] cell array of string for the legend. Note
%                the last element can be an integer to set legend 
%                position.
%  'showleg'   = ['on'|'off'] show or hide legend.
%  'colors'    = [cell array] cell array of plot aspect. E.g. { 'k' 'k--' }
%                for plotting the first curve in black and the second one
%                in black dashed. Can also contain additional formating.
%                { { 'k' 'fontweight' 'bold' } 'k--' } same as above but
%                the first line is bolded.
%  'ydir'      = [1|-1] y-axis polarity (pos-up = 1; neg-up = -1) {def -> 1}
%  'vert'      = [vector] of times (in ms or Hz) to plot vertical lines 
%                {def none}
%  'regions'   = [float array] float array of size (n,2) each line defining 
%                a time region [low high] to be highlighted.
%  'axsize'    = [x y] axis size {default [.07 .07]}
%
% Author: Scott Makeig and Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 3-2-98 
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
% Revision 1.15  2003/03/16 02:59:30  arno
% debug legend
%
% Revision 1.14  2003/03/16 02:43:50  arno
% allowing to show legend
%
% Revision 1.13  2003/03/16 01:49:09  arno
% debug
% last
%
% Revision 1.12  2003/03/16 01:47:24  arno
% allowing setting of each curve color and aspect
%
% Revision 1.11  2003/03/16 01:22:16  arno
% converting to 'key' 'val' parameter sequence
%
% Revision 1.10  2003/03/05 16:33:29  arno
% default linewidth set to 1
%
% Revision 1.9  2003/03/05 16:27:53  arno
% plotting lines after data
%
% Revision 1.8  2003/03/05 02:25:23  arno
% removing warnings and extra CR
%
% Revision 1.7  2003/02/21 00:36:09  scott
% header edit -sm
%
% Revision 1.6  2002/07/22 22:57:54  arno
% swap 2 first colors
%
% Revision 1.5  2002/04/24 18:24:39  scott
% added vert to cal axis -sm
%
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

function plottopo(data, varargin);
    
%
%%%%%%%%%%%%%%%%%%%%% Graphics Settings - can be customized %%%%%%%%%%%%%%%%%%
%
LINEWIDTH     = 1.0;     % data line widths (can be non-integer)
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
    
if nargin < 1
    help plottopo
    return
end

if length(varargin) > 0
    if length(varargin) == 1 | ~isstr(varargin{1}) | isempty(varargin{1})
        options = { 'chanlocs' varargin{1} };
        if nargin > 2, options = { options{:} 'frames' varargin{2} }; end;
        if nargin > 3, options = { options{:} 'limits' varargin{3} }; end;
        if nargin > 5, options = { options{:} 'chans'  varargin{5} }; end;
        if nargin > 6, options = { options{:} 'axsize' varargin{6} }; end;
        if nargin > 7, options = { options{:} 'colors' varargin{7} }; end;
        if nargin > 8, options = { options{:} 'ydir'   varargin{8} }; end;
        if nargin > 9, options = { options{:} 'vert'   varargin{9} }; end;
        if nargin > 4 & ~isequal(varargin{4}, 0), options = { options{:} 'title'  varargin{4} }; end;
        
        %    , chan_locs,frames,limits,plottitle,channels,axsize,colors,ydr,vert)
    else
        options = varargin;
    end;
else
    options = varargin;
end;
g = finputcheck(options, { 'chanlocs'  ''    []          '';
                    'frames'    'integer'               [1 Inf]     size(data,2);
                    'chans'     'integer'               [1 Inf]     0;
                    'geom'      'integer'               [1 Inf]     [];
                    'limits'    'float'                 []          0;
                    'title'     'string'                []          '';
                    'axsize'    'float'                 [0 1]       [nan nan];
                    'regions'   'float'                 []          [];
                    'colors'    'cell'                  []          {};
                    'legend'    'cell'                  []          {};
                    'showleg'   'string'                {'on' 'off'} 'on';
                    'ydir'      'integer'               [-1 1]      DEFAULT_SIGN;
                    'vert'      'float'                 []          []});
if isstr(g), error(g); end;
data = reshape(data, size(data,1), size(data,2), size(data,3));    
if length(g.chans) == 1 & g.chans(1) ~= 0, error('can not plot a single ERP'); end;

[chans,framestotal]=size(data);           % data size

%
%%%%%%%%%%%%%%% Substitute defaults for missing parameters %%%%%%%%%%%%%%%%
%
  
axwidth  = g.axsize(1);
axheight = g.axsize(2);

if ~isempty(g.chans) & g.chans == 0
   channelnos = 1:size(data,1);
elseif ~isstr(g.chans)
   channelnos = g.chans;
else
   channelnos = 1:size(data,1);
end

nolegend = 0;
if isempty(g.legend), nolegend = 1; end;

limitset = 0;
if length(g.limits)>1
    limitset = 1;
end

if isempty(g.chanlocs) & isempty(g.geom)
  n = ceil(sqrt(length(channelnos)));
  g.geom = [n n];
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Test parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  icadefs; % read BACKCOLOR, MAXPLOTDATACHANS constant from icadefs.m
  if g.frames <=0,
    g.frames = framestotal;    % default
    datasets=1;
  elseif g.frames==1,
    fprintf('plottopo: cannot plot less than 2 frames per trace.\n');
    return
    datasets=1;
  else
    datasets = fix(framestotal/g.frames);        % number of traces to overplot
  end;

  if max(channelnos) > chans
    fprintf('plottopo(): max channel index > %d channels in data.\n',...
                       chans);
    return
  end
  if min(channelnos) < 1
    fprintf('plottopo(): min channel index (%g) < 1.\n',...
                       min(g.chans));
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

  if ~isempty(g.geom)
      if isnan(axheight) % if not specified
          axheight = gcapos(4)/(g.geom(1)+1);
          axwidth  = gcapos(3)/(g.geom(2)+1);
      end
      % if chan_locs(2) > 5
      %     axwidth = 0.66/(chan_locs(2)+1);
      % end
  else
      axheight = DEFAULT_AXHEIGHT;
      axwidth =  DEFAULT_AXWIDTH;
  end
    fprintf('Plotting data using axis size [%g,%g]\n',axwidth,axheight);

%
%%%%%%%%%%%%% Extend the size of the plotting area in the window %%%%%%%%%%%%
%
  curfig = gcf;
  h=figure(curfig);
  set(h,'PaperUnits','normalized'); % use percentages to avoid US/A4 difference
  set(h,'PaperPosition',[0.0235308 0.0272775 0.894169 0.909249]); % equivalent
  orient portrait
  axis('normal');

  set(gca,'Color',BACKCOLOR);               % set the background color
  
  axcolor= get(0,'DefaultAxesXcolor'); % find what the default x-axis color is
  vertcolor = 'b';
  plotfile = 'plottopo.ps';
  ls_plotfile = 'ls -l plottopo.ps';
    
%
%%%%%%%%%%%%%%%%%%%% Read the channel names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  if ~isstr(g.chans) 
    % channames = zeros(MAXPLOTDATACHANS,4);
    % for c=1:length(g.chans),
    %     channames(c,:)= sprintf('%4d',g.chans(c));
    % end;
    if length(g.chans) > 1 | g.chans(1) ~= 0
        channames = num2str(g.chans(:));                   %%CJH
    end;
  else % isstr(g.chans)
    if ~isstr(g.chans)
       fprintf('plottopo(): channel file name must be a string.\n');
       return
    end
    chid = fopen(g.chans,'r');
    if chid <3,
        fprintf('plottopo(): cannot open file %s.\n',g.chans);
        return
    else
        fprintf('plottopo(): opened file %s.\n',g.chans);
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
  if isstr(g.colors)
      cid = fopen(g.colors,'r');
      % fprintf('cid = %d\n',cid);
      if cid <3,
          fprintf('plottopo: cannot open file %s.\n',g.colors);
          return
      end;
      g.colors = fscanf(cid,'%s',[3 MAXPLOTDATAEPOCHS]);
      g.colors = g.colors';
      [r c] = size(g.colors);
      for i=1:r
          for j=1:c
              if g.colors(i,j)=='.',
                  g.colors(i,j)=' ';
              end;
          end;
      end;
      g.colors = cellstr(g.colors);
      for c=1:length(g.colors)   % make white traces black unless axis color is white
          if g.colors{c}(1)=='w' & axcolor~=[1 1 1]
              g.colors{c}(1)='k';
          end
      end
  else % use default color order (no yellow!)
      tmpcolors = { 'b' 'r' 'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm' 'r' 'b' ...
                   'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm'};
      g.colors = {g.colors{:} tmpcolors{:} tmpcolors{:}};  % make > 64 available
  end;
%
%%%%%%%%%%%%%%%%%%%%%%% Read and adjust limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  if g.limits==0,      % == 0 or [0 0 0 0]
    xmin=0;
    xmax=g.frames-1;
    ymin=min(min(data));
    ymax=max(max(data));
  else
    if length(g.limits)~=4,
      fprintf( ...
       'plottopo: limits should be 0 or an array [xmin xmax ymin ymax].\n');
      return
    end;
    xmin = g.limits(1);
    xmax = g.limits(2);
    ymin = g.limits(3);
    ymax = g.limits(4);
  end;

  if xmax == 0 & xmin == 0,
    x = (0:1:g.frames-1);
    xmin = 0;
    xmax = g.frames-1;
  else
    dx = (xmax-xmin)/(g.frames-1);
    x=xmin*ones(1,g.frames)+dx*(0:g.frames-1); % compute x-values
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
    g.ydir = 1;
    fprintf('Plotting positive up. Assuming data are spectra.\n');
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
  h=gca;title(g.title,'FontSize',TITLEFONTSIZE); % title plot and
  hold on
  msg = ['Plotting %d traces of %d frames with colors: '];

  for c=1:datasets
      if iscell(g.colors{c})
          msg = [msg  '''' g.colors{c}{1} ''' ' ];
      else
          msg = [msg  '''' g.colors{c} ''' ' ];
      end;
  end
  msg = [msg '\n'];    % print starting info on screen . . .
  fprintf('limits: [xmin,xmax,ymin,ymax] = [%4.1f %4.1f %4.2f %4.2f]\n',...
                xmin,xmax,ymin,ymax);
  fprintf(msg,datasets,g.frames);

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
if isempty(g.chanlocs) % plot in a rectangular grid
    ISRECT = 1;
    ht = g.geom(1);
    wd = g.geom(2);
    if chans > ht*wd
        fprintf('topoplot(): (%d) channels to be plotted > grid size [%d %d]\n',...
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
	[tmp channames Th Rd] = readlocs(g.chanlocs);
	channames = strvcat(channames{channelnos});
	Th = pi/180*Th(channelnos);                 % convert degrees to radians
	Rd = Rd(channelnos); 
	
    if length(channelnos) > length(Th),
        error('plottopo(): data channels must be <= ''chanlocs'' channels')
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
      fprintf('trace %d: ',P+1);

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

            %
            %%%%%%%%%%%%%%%%%%%%%%% Highlight regions %%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            for index=1:size(g.regions,1)
                tmpreg = g.regions(index,:);
                tmph = patch([tmpreg(1) tmpreg(2) tmpreg(2) tmpreg(1)], ...
                             [-100 -100 100 100], [1 1 0.9]); hold on;
                set(tmph, 'edgecolor', [1 1 0.9]);
            end;
            
            % secondx = 200;                             % draw second vert axis 
            % axis('off');plot([secondx secondx],[ymin ymax],'color',axislcolor); 
            %
            %%%%%%%%%%%%%%%%%%%%%%% Print channel names %%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            NAME_OFFSET = -1;
            NAME_OFFSETY = -0.5;
            if ISSPEC
                axis('off'),h=text(xmin-NAME_OFFSET*xdiff,ymax/2,[channames(c,:)]); 
                set(h,'HorizontalAlignment','right');    % print before traces
                set(h,'FontSize',CHANFONTSIZE);              % choose font size
            else % ~ISSPEC
                if ymin <= 0 & ymax >= 0,
                    yht = 0;
                else
                    yht = mean(g.ydir*data(c,1+P*g.frames:1+P*g.frames+g.frames-1));
                end
                if ~ISRECT    % print before traces
                    axis('off'),h=text(xmin-NAME_OFFSET*xdiff,yht-NAME_OFFSETY*ydiff,[channames(c,:)]); 
                    set(h,'HorizontalAlignment','right');      
                    set(h,'FontSize',CHANFONTSIZE);           % choose font size
                else % ISRECT
                    xmn = xdiff/2+xmin;
                    axis('off'),h=text(xmn,ymax+0.05*ymax,[channames(c,:)]); 
                    set(h,'HorizontalAlignment','right');      
                    set(h,'FontSize',CHANFONTSIZE);            % choose font size
                end % ISRECT
            end % ~ISSPEC
            
            
        end; % P=0 
        %
        %%%%%%%%%%%%%%%%%%%%%%% Plot data traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if ~iscell( g.colors{P+1} ), tmpcolor = { g.colors{P+1} 'linewidth' LINEWIDTH };
        else                         tmpcolor = g.colors{P+1};
        end;
        if ~ISSPEC % -/+ plot, normal case (e.g., not spectra), plot data trace           
            plot(x,g.ydir*data(c,1+P*g.frames:1+P*g.frames+g.frames-1), tmpcolor{:});   
            ymn = min(g.ydir*[ymax ymin]);
            ymx = max(g.ydir*[ymax ymin]);
            axis([xmin xmax ymn ymx]);          % set axis bounds
        else % ISSPEC
            plot(x,g.ydir*data(c,1+P*g.frames:1+P*g.frames+g.frames-1), tmpcolor{:});   
            ymaxm = ymax;
            if ymaxm/2. > ymax,
                ymaxm = ymaxm/2.;
            end;
            axis([xmin xmax ymin ymaxm]);      % set axis values
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%% Plot lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if P == datasets-1
            if ISSPEC
                plot([xmin xmin],[0 ymax],'color',axislcolor); 
            else
                plot([0 0],[ymin ymax],'color',axislcolor); % draw vert axis at time 0  
            end  
            axis('off');
            plot([xmin xmax],[0 0],'color',axislcolor);  % draw horizontal axis 
        end;
        %
        %%%%%%%%%%%%%%%%%%%% plot vertical lines (optional) %%%%%%%%%%%%%%%%%
        %
        if ~isnan(g.vert)
            if ~ISSPEC % -/+ plot, normal case (e.g., not spectra), plot data trace
                ymean = (ymin+ymax)/2; 
                vmin = ymean-0.5*(ymean-ymin);
                vmax = ymean+0.5*(ymax-ymean);
                for v = g.vert
                    plot([v v],[vmin vmax],'color',vertcolor); % draw vertical lines 
                end
            else
                for v = g.vert
                    plot([v v],[0 ymax],'color',vertcolor); 
                end
            end
        end
        
        fprintf(' %d',c); % finished with channel plot
    end; % c, chans / subplot
    % handle legend
    if nolegend, g.legend{P+1} = ['Data ' int2str(P) ]; end;
        
    fprintf('\n');
  end; % P / epoch

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
  if ~isnan(g.vert)
   if ~ISSPEC % -/+ plot, normal case (e.g., not spectra), plot data trace
    for v = g.vert
      plot([v v],[vmin vmax],'color',vertcolor); % draw vertical lines 
    end
   else
    for v = g.vert
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
    axis('off');h=text(signx,g.ydir*ymin,num2str(ymin,3)); % text ymin
    set(h,'FontSize',TICKFONTSIZE);               % choose font size
    set(h,'HorizontalAlignment','right','Clipping','off');

    axis('off');h=text(signx,g.ydir*ymax,['+' num2str(ymax,3)]);  % text +ymax
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

    if length(g.legend) > 1 & strcmpi(g.showleg, 'on')
        tmpleg = vararg2str(g.legend);
        quotes = find(tmpleg == '''');
        for index = length(quotes):-1:1
            tmpleg(quotes(index)+1:end+1) = tmpleg(quotes(index):end);
            tmpleg(quotes(index)) = '''';
        end;
        tmpleg = [ 'legend(' tmpleg ');' ];
    else tmpleg = '';
    end;
    com = [ 'axis on;' ...
            'clear xlabel ylabel;' tmpleg ...
            'xlabel(''''Time (ms)'''');' ...
            'ylabel(''''Voltage (\muV)'''');' ];
    axcopy(gcf, com); % turn on popup feature
%
%%%%%%%%%%%%%%%%%% Make printed figure fill page %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
orient tall
  % curfig = gcf;
  % h=figure(curfig);
  % set(h,'PaperPosition',[0.2 0.3 7.6 10]); % stretch out the plot on the page

