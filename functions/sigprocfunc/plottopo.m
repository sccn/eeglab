% plottopo() - plot concatenated multichannel data epochs in a topographic
% or
%              rectangular array. Uses a channel location file with the same 
%              format as topoplot(), or else plots data on a rectangular grid. 
%              If data are all positive, they are assumed to be spectra.
% Usage:
%    >> plottopo(data, 'key1', 'val1', 'key2', 'val2')
% Or
%    >> plottopo(data,'chan_locs',frames,limits,title,channels,...
%                      axsize,colors,ydir,vert) % old function call
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
%  'ylim'      = [ymin ymax] y axis limits. Overwrite option above.
%  'title'     = [string] plot title {def|'' -> none}
%  'chans'     = vector of channel numbers to plot {def|0 -> all}
%  'axsize'    = [x y] axis size {default [.05 .08]}
%  'legend'    = [cell array] cell array of string for the legend. Note
%                the last element can be an integer to set legend 
%                position.
%  'showleg'   = ['on'|'off'] show or hide legend.
%  'colors'    = [cell array] cell array of plot aspect. E.g. { 'k' 'k--' }
%                for plotting the first curve in black and the second one
%                in black dashed. Can also contain additional formating.
%                { { 'k' 'linewidth' 2 } 'k--' } same as above but
%                the first line is bolded.
%  'ydir'      = [1|-1] y-axis polarity (pos-up = 1; neg-up = -1) {def -> -1}
%  'vert'      = [vector] of times (in ms or Hz) to plot vertical lines 
%                {def none}
%  'hori'      = [vector] plot horizontal line at given ordinate values.
%  'regions'   = [cell array] cell array of size nchan. Each cell contains a
%                float array of size (2,n) each column defining a time region 
%                [low high] to be highlighted.
%  'plotfunc'  = [cell] use different function for plotting data. The format
%                is { funcname arg2 arg3 arg2 ... }. arg1 is taken from the
%                data.
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
% Revision 1.59  2008/04/19 21:14:32  arno
% do not draw axis while plotting
%
% Revision 1.58  2007/08/09 22:22:57  arno
% typo
%
% Revision 1.57  2007/08/07 19:06:16  arno
% fix bug 334 and 330; remove channelnos variable
%
% Revision 1.56  2007/08/07 19:01:02  arno
% fixed bug 321
%
% Revision 1.55  2007/08/07 18:59:32  arno
% fixed unused variables bug 319
%
% Revision 1.54  2007/08/07 18:58:10  arno
% bug 307
%
% Revision 1.53  2007/08/06 17:55:04  arno
% remove ISSPEC, clean code and fix bug 328
%
% Revision 1.52  2007/08/06 17:42:24  arno
% 'chans' can now be a string
%
% Revision 1.51  2007/08/06 17:38:36  arno
% fix bug 327
%
% Revision 1.50  2007/02/09 01:41:31  toby
% 3rd attempt, looks good so far
%
% Revision 1.49  2007/02/09 00:31:39  toby
% 2nd attempt to fix bug, tested well
%
% Revision 1.48  2007/02/08 23:42:35  toby
% correcting bug discovered by nightly script ICADEMO
%
% Revision 1.47  2007/02/05 16:19:12  arno
% Bug due to Stefan's changes, legend color now fixed
%
% Revision 1.46  2007/02/02 16:36:04  arno
% speed plot
%
% Revision 1.42  2006/11/15 20:41:40  arno
% Voltage -> Potential
%
% Revision 1.41  2006/07/10 21:00:21  arno
% fix one channel plotting
%
% Revision 1.40  2006/03/09 17:20:39  arno
% stefan's changes
%
% Revision 1.39  2005/11/16 22:01:41  toby
% Edited vertical patch bars to be transparent for indicating regions of interest
% so that the plots will be visible underneath. Another solution would have been to
% plot the patches first.
%
% Revision 1.39  2005/11/24 10:32:0  stefan
% added horizontal line input, modified defaults, added axis ticks, changed
% plot order (data trace now on top)
%
% Revision 1.38  2005/03/21 16:13:09  arno
% highlight regions of significance at the end
%
% Revision 1.37  2004/09/14 09:30:04  arno
% fix matlab 7
%
% Revision 1.36  2004/09/08 20:55:32  scott
% same
%
% Revision 1.35  2004/09/08 02:46:47  scott
% same
%
% Revision 1.34  2004/09/08 02:09:39  scott
% insert figure(curfig) before plotting commands to squelch Matlab 7.0.0 bug -sm
%
% Revision 1.33  2004/04/06 17:21:38  arno
% draw full vertical lines
%
% Revision 1.32  2004/02/10 16:56:27  arno
% error msg
%
% Revision 1.31  2004/01/29 16:44:32  arno
% fix icademo bug
%
% Revision 1.30  2004/01/29 00:58:04  arno
% allowing to plot channels with no coordinates
%
% Revision 1.29  2003/09/17 01:44:57  arno
% debuging for icademo
%
% Revision 1.28  2003/07/25 17:24:43  arno
% allow to plot more than 50 trials
%
% Revision 1.27  2003/07/16 00:58:35  arno
% debug legnd
%
% Revision 1.26  2003/07/16 00:38:14  arno
% fixing ydir
%
% Revision 1.25  2003/07/16 00:31:18  arno
% debug ydir
%
% Revision 1.24  2003/07/16 00:26:49  arno
% debug legend
%
% Revision 1.23  2003/07/16 00:23:22  arno
% debug ydir
%
% Revision 1.22  2003/07/15 18:37:29  arno
% debug color
%
% Revision 1.21  2003/07/15 17:11:22  arno
% allowing empty g.chans & header typo
%
% Revision 1.20  2003/05/09 23:26:06  arno
% debuging regions
%
% Revision 1.19  2003/03/17 23:38:08  arno
% programing ylim option
%
% Revision 1.18  2003/03/17 23:33:28  arno
% debuging regions of interest
%
% Revision 1.17  2003/03/17 23:05:05  arno
% debuging regions
%
% Revision 1.16  2003/03/17 22:04:41  arno
% allow to highlight regions of interest
%
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
%   hori        = [vector] of amplitudes (in uV or dB) to plot horizontal lines {def none}
%

function plottopo(data, varargin);
    
%
%%%%%%%%%%%%%%%%%%%%% Graphics Settings - can be customized %%%%%%%%%%%%%%%%%%
%
LINEWIDTH     = 0.7;     % data line widths (can be non-integer)
FONTSIZE      = 10;      % font size to use for labels
CHANFONTSIZE  = 7;       % font size to use for channel names
TICKFONTSIZE  = 8;       % font size to use for axis labels
TITLEFONTSIZE = 12;      % font size to use for the plot title
PLOT_WIDTH    = 0.95;     % 0.75, width and height of plot array on figure
PLOT_HEIGHT   = 0.88;    % 0.88
gcapos = get(gca,'Position'); axis off;
PLOT_WIDTH    = gcapos(3)*PLOT_WIDTH; % width and height of gca plot array on gca
PLOT_HEIGHT   = gcapos(4)*PLOT_HEIGHT;
curfig = gcf;            % learn the current graphic figure number
%
%%%%%%%%%%%%%%%%%%%% Default settings - use commandline to override %%%%%%%%%%%
%
DEFAULT_AXWIDTH  = 0.05; %
DEFAULT_AXHEIGHT = 0.08; % 
DEFAULT_SIGN = -1;                        % Default - plot positive-up
ISRECT = 0;                               % default
    
if nargin < 1
    help plottopo
    return
end

if length(varargin) > 0
    if length(varargin) == 1 | ~isstr(varargin{1}) | isempty(varargin{1}) | ...
        (length(varargin)>2 &  ~isstr(varargin{3}))
        options = { 'chanlocs' varargin{1} };
        if nargin > 2, options = { options{:} 'frames' varargin{2} }; end;
        if nargin > 3, options = { options{:} 'limits' varargin{3} }; end;
        if nargin > 5, options = { options{:} 'chans'  varargin{5} }; end;
        if nargin > 6, options = { options{:} 'axsize' varargin{6} }; end;
        if nargin > 7, options = { options{:} 'colors' varargin{7} }; end;
        if nargin > 8, options = { options{:} 'ydir'   varargin{8} }; end;
        if nargin > 9, options = { options{:} 'vert'   varargin{9} }; end;
        if nargin > 10,options = { options{:} 'hori'  varargin{10} }; end;
        if nargin > 4 & ~isequal(varargin{4}, 0), options = {options{:} 'title'  varargin{4} }; end;
        %    , chan_locs,frames,limits,plottitle,channels,axsize,colors,ydr,vert)
    else
        options = varargin;
    end;
else
    options = varargin;
end;
g = finputcheck(options, { 'chanlocs'  ''    []          '';
                    'frames'    'integer'               [1 Inf]     size(data,2);
                    'chans'     { 'integer' 'string' }  { [1 Inf] [] }    0;
                    'geom'      'integer'               [1 Inf]     [];
                    'limits'    'float'                 []          0;
                    'ylim'      'float'                 []          [];
                    'title'     'string'                []          '';
                    'plotfunc'  'cell'                  []          {};
                    'axsize'    'float'                 [0 1]       [nan nan];
                    'regions'   'cell'                  []          {};
                    'colors'    { 'cell' 'string' }     []          {};
                    'legend'    'cell'                  []          {};
                    'showleg'   'string'                {'on' 'off'} 'on';
                    'ydir'      'integer'               [-1 1]      DEFAULT_SIGN;
                    'vert'      'float'                 []          [];
                    'hori'      'float'                 []          []});
if isstr(g), error(g); end;
data = reshape(data, size(data,1), size(data,2), size(data,3));    
%if length(g.chans) == 1 & g.chans(1) ~= 0, error('can not plot a single ERP'); end;

[chans,framestotal]=size(data);           % data size

%
%%%%%%%%%%%%%%% Substitute defaults for missing parameters %%%%%%%%%%%%%%%%
%
  
axwidth  = g.axsize(1);
if length(g.axsize) < 2
    axheight = NaN;
else 
    axheight = g.axsize(2);
end;
if isempty(g.chans) | g.chans == 0
   g.chans = 1:size(data,1);
elseif ~isstr(g.chans)
   g.chans = g.chans;
end

nolegend = 0;
if isempty(g.legend), nolegend = 1; end;

if ~isempty(g.ylim)
    g.limits(3:4) = g.ylim;
end;
plotgrid = 0;
if isempty(g.chanlocs) % plot in a rectangular grid
    plotgrid = 1;
elseif ~isfield(g.chanlocs, 'theta')
    plotgrid = 1;
end;
if length(g.chans) < 4 & ~plotgrid
    disp('Not enough channels, does not use channel coordinate to plot axis');
    plotgrid = 1;
end;
if plotgrid & isempty(g.geom)
  n = ceil(sqrt(length(g.chans)));
  g.geom = [n ceil(length(g.chans)/n)];
end
if ~isempty(g.geom)
    plotgrid = 1;
end;
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

  if max(g.chans) > chans
    fprintf('plottopo(): max channel index > %d channels in data.\n',...
                       chans);
    return
  end
  if min(g.chans) < 1
    fprintf('plottopo(): min channel index (%g) < 1.\n',...
                       min(g.chans));
    return
  end;
  if length(g.chans)>MAXPLOTDATACHANS,
    fprintf('plottopo(): not set up to plot more than %d traces.\n',...
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
    vertcolor = 'k';
    horicolor = vertcolor;
    
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
    data = data(g.chans,:);
    chans = length(g.chans);
    %
    %%%%%%%%%%%%%%%%%%%%%%%%% Read the color names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if isstr(g.colors) % filename for backward compatibility but not documented
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
        g.colors = {g.colors{:} tmpcolors{:}};  % make > 64 available
    end;
    %
    %%%%%%%%%%%%%%%%%%%%%%% Read and adjust limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if g.limits==0,      % == 0 or [0 0 0 0]
        xmin=0;
        xmax=g.frames-1;
        % for abs max scaling:
        ymax=max(max(abs(data)));
        ymin=ymax*-1;
        % for data limits:
        %ymin=min(min(data));
        %ymax=max(max(data));
    else
        if length(g.limits)~=4,
            error('plottopo: limits should be 0 or an array [xmin xmax ymin ymax].\n');
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
        % for abs max scaling:
        ymax=max(max(abs(data)));
        ymin=ymax*-1;
        % for data limits:
        %ymin=min(min(data));
        %ymax=max(max(data));
    end
    if ymax<=ymin,
        fprintf('plottopo() - ymax must be > ymin.\n')
        return
    end

    xlabel = 'Time (ms)';
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
    figure(curfig); h=gca;title(g.title,'FontSize',TITLEFONTSIZE); % title plot 
    hold on
    msg = ['Plotting %d traces of %d frames with colors: '];

    for c=1:datasets
        cind = mod(c-1, length(g.colors))+1;
        if iscell(g.colors{cind})
            msg = [msg  '''' g.colors{cind}{1} ''' ' ];
        else
            msg = [msg  '''' g.colors{cind} ''' ' ];
        end;
    end
    msg = [msg '\n'];    % print starting info on screen . . .
    fprintf('limits: [xmin,xmax,ymin,ymax] = [%4.1f %4.1f %4.2f %4.2f]\n',...
            xmin,xmax,ymin,ymax);
    fprintf(msg,datasets,g.frames);

    set(h,'FontSize',FONTSIZE);           % choose font size
    set(h,'YLim',[ymin ymax]);            % set default plotting parameters
    set(h,'XLim',[xmin xmax]);

    axis('off')
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Read chan_locs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if plotgrid
        ISRECT = 1;
        ht = g.geom(1);
        wd = g.geom(2);
        if chans > ht*wd
            fprintf('plottopo(): (%d) channels to be plotted > grid size [%d %d]\n',...
                    chans,ht,wd);
            return
        end
        xvals = 0; yvals = 0;
        if ~exist('channames') 
            if isfield(g.chanlocs,'labels') && ~iscellstr({g.chanlocs.labels})
                channames = strvcat(g.chanlocs.labels);
            else
                channames = repmat(' ',ht*wd,4);
                for i=1:ht*wd
                    channum = num2str(i);
                    channames(i,1:length(channum)) = channum;
                end
            end
        end
        
    else % read chan_locs file
         % read the channel location file
         % ------------------------------
        if isstruct(g.chanlocs)
            nonemptychans = cellfun('isempty', { g.chanlocs.theta });
            nonemptychans = find(~nonemptychans);
            [tmp channames Th Rd] = readlocs(g.chanlocs(nonemptychans));
            channames = strvcat({ g.chanlocs.labels });
        else
            [tmp channames Th Rd] = readlocs(g.chanlocs);
            channames = strvcat(channames);
            nonemptychans = [1:length(channames)];
        end;
        Th = pi/180*Th;                 % convert degrees to radians
        Rd = Rd; 
        
        if length(g.chans) > length(g.chanlocs),
            error('plottopo(): data channels must be <= ''chanlocs'' channels')
        end
        
        [yvalstmp,xvalstmp] = pol2cart(Th,Rd); % translate from polar to cart. coordinates
        xvals(nonemptychans) = xvalstmp;
        yvals(nonemptychans) = yvalstmp;
        
        % find position for other channels
        % --------------------------------
        totalchans = length(g.chanlocs);
        emptychans = setdiff(1:totalchans, nonemptychans);
        totalchans = floor(sqrt(totalchans))+1;
        for index = 1:length(emptychans)
            xvals(emptychans(index)) = 0.7+0.2*floor((index-1)/totalchans);
            yvals(emptychans(index)) = -0.4+mod(index-1,totalchans)/totalchans;
        end;
        channames = channames(g.chans,:);
        xvals     = xvals(g.chans);
        yvals     = yvals(g.chans);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % xvals = 0.5+PLOT_WIDTH*xvals;   % controls width of  plot array on page!
    % yvals = 0.5+PLOT_HEIGHT*yvals;  % controls height of plot array on page!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if length(xvals) > 1
        if length(unique(xvals)) > 1
            xvals = (xvals-mean([max(xvals) min(xvals)]))/(max(xvals)-min(xvals)); % recenter
            xvals = gcapos(1)+gcapos(3)/2+PLOT_WIDTH*xvals;   % controls width of plot 
                                                              % array on current axes
        end;
    end;
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
                if plotgrid
                    Axes = [ Axes sbplot(g.geom(1), g.geom(2), c)];
                else
                    xcenter = xvals(c);
                    ycenter = yvals(c);
                    Axes = [Axes axes('Units','Normal','Position', ...
                                      [xcenter-axwidth/2 ycenter-axheight/2 axwidth axheight])];
                end;
                %axes(Axes(c))
                axis('off')
                
                hold on;                      % plot down left side of page first
                                              % set(h,'YLim',[ymin ymax]);    % set default plotting parameters
                                              % set(h,'XLim',[xmin xmax]);
                
                axislcolor = get(gca,'Xcolor');   %%CJH
                
                axis('off');

                % secondx = 200;                             % draw second vert axis 
                % axis('off');plot([secondx secondx],[ymin ymax],'color',axislcolor); 
                %
                %%%%%%%%%%%%%%%%%%%%%%% Print channel names %%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                NAME_OFFSET = -.25;
                NAME_OFFSETY = .2;
                if ymin <= 0 & ymax >= 0,
                    yht = 0;
                else
                    yht = mean(data(c,1+P*g.frames:1+P*g.frames+g.frames-1));
                end
                if ~ISRECT    % print before traces
                    xt = double(xmin-NAME_OFFSET*xdiff);
                    yt = double(yht-NAME_OFFSETY*ydiff); 
                    str = [channames(c,:)];
                    h=text(xt,yt,str);
                    set(h,'HorizontalAlignment','right');      
                    %set(h,'FontSize',CHANFONTSIZE);           % choose font size
                else % ISRECT
                    xmn = xdiff/2+xmin;
                    h=text(double(xmn),double(ymax+0.05*ymax),[channames(c,:)]); 
                    set(h,'HorizontalAlignment','right');      
                    %set(h,'FontSize',CHANFONTSIZE);            % choose font size
                end % ISRECT
                
                
                %
                %%%%%%%%%%%%%%%%%%%%%%% Highlight regions %%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                if ~isempty(g.regions)
                    for index=1:size(g.regions{c},2)
                        tmpreg = g.regions{c}(:,index);
                        if tmpreg(1) ~= tmpreg(2)
                            tmph = patch([tmpreg(1) tmpreg(2) tmpreg(2) tmpreg(1)], ...
                                         [-100 -100 100 100], [0.9 0.9 0.9]); hold on;
                            set(tmph, 'edgecolor', [0.9 0.9 0.9]); %,'facealpha',0.5,'edgealpha',0.5);
                        end;
                    end;
                end;
                
            end; % P=0 
            
            %
            %%%%%%%%%%%%%%%%%%%%%%% Plot data traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            Pind = mod(P+1-1, length(g.colors))+1;
            if ~iscell( g.colors{Pind} ), tmpcolor = { g.colors{Pind} 'linewidth' LINEWIDTH };
            else                          tmpcolor = g.colors{Pind};
            end;
            ymn = min([ymax ymin]);
            ymx = max([ymax ymin]);
            if isempty(g.plotfunc)
                if isstr(tmpcolor{1}) & length(tmpcolor) > 1
                    plot(x,data(c,1+P*g.frames:1+P*g.frames+g.frames-1), tmpcolor{1}, tmpcolor{2:end});   
                else
                    plot(x,data(c,1+P*g.frames:1+P*g.frames+g.frames-1), 'color', tmpcolor{:});   
                end; 
                if g.ydir == -1
                    set(gca, 'ydir', 'reverse');
                end;
                axis([xmin xmax ymn ymx]);          % set axis bounds
            elseif P == 1
                func = eval( [ '@' g.plotfunc{1} ] );
                feval(func, data(c,:), g.plotfunc{2:end});
            end;
            
            if P == datasets-1 % last pass
                               %
                               %%%%%%%%%%%%%%%%%%%%%%% Plot lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                               %
                plot([0 0],[ymin ymax],'color',axislcolor); % draw vert axis at time 0  
                axis('off');
                plot([xmin xmax],[0 0],'color',axislcolor);  % draw horizontal axis 

                %
                %%%%%%%%%%%%%%%%%%%% plot vertical lines (optional) %%%%%%%%%%%%%%%%%
                %
                
                if isempty(g.vert)
                    g.vert = [xmin xmax];
                    ymean = (ymin+ymax)/2; 
                    vmin = ymean-0.1*(ymean-ymin);
                    vmax = vmin*-1;  %ymean+0.2*(ymax-ymean);
                elseif ~isnan(g.vert)
                    ymean = (ymin+ymax)/2; 
                    vmin = ymean-0.1*(ymean-ymin);
                    vmax = vmin*-1;  %ymean+0.2*(ymax-ymean);
                    for v = g.vert
                        plot([v v],[vmin vmax],'color',vertcolor); % draw vertical lines 
                    end
                end
                
                %
                %%%%%%%%%%%%%%%%%%%% plot horizontal lines (optional) %%%%%%%%%%%%%%%
                %
                if isempty(g.hori)
                    g.hori = [ymin ymax]; 
                end
                if ~isnan(g.hori)
                    xmean = 0; 
                    hmin = xmean-0.2*(xmean-xmin);
                    hmax = hmin*-1; %xmean+0.3*(xmax-xmean);
                    for v = g.hori
                        plot([hmin hmax],[v v], 'color',horicolor); % draw horizontal lines 
                    end
                end

            end;
            
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
              [0.85 0.1 axwidth axheight]); % FIX!!!!
    axes(ax)
    axis('off');
    if xmin <=0
        figure(curfig);p=plot([0 0],[ymn ymx],'color','k'); % draw vert axis at zero
    else
        figure(curfig);p=plot([xmin xmin],[ymn ymx],'color','k'); % draw vert axis at zero
    end
    if g.ydir == -1
        set(gca, 'ydir', 'reverse');
    end;
    axis([xmin xmax ymn ymx]);        % set axis values
    hold on
    %set(p, 'Clipping','off');        % center text
    figure(curfig);p=plot([xmin xmax],[0 0],'color',axislcolor); % draw horizontal axis 
    axis([xmin xmax ymin ymax]);        % set axis values
                                        %
                                        %%%%%%%%%%%%%%%%%%%% plot vertical lines (optional) %%%%%%%%%%%%%%%%%
                                        %
    if ~isnan(g.vert)
        for v = g.vert
            figure(curfig);plot([v v],[vmin vmax],'color',vertcolor); % draw vertical lines 
        end
    end
    %
    %%%%%%%%%%%%%%%%%%%% plot horizontal lines (optional) %%%%%%%%%%%%%%%%%
    %
    if ~isnan(g.hori)
        xmean = 0; 
        hmin = xmean-0.2*(xmean-xmin);
        hmax = hmin*-1; %xmean+0.3*(xmax-xmean);
        for v = g.hori
            figure(curfig);plot([hmin hmax],[v v], 'color',horicolor); % draw horizontal lines 
        end
    end
    
    % secondx = 200;                    % draw second vert axis 
    % axis('off');plot([secondx secondx],[ylo ymax],'color',axislcolor); 
    %
    %%%%%%%%%%%%%%%%%%%%% Plot negative-up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    signx = xmin-0.15*xdiff;
    figure(curfig);axis('off');h=text(double(signx),double(ymin),num2str(ymin,3)); % text ymin
    set(h,'FontSize',TICKFONTSIZE);               % choose font size
    set(h,'HorizontalAlignment','right','Clipping','off');

    figure(curfig);axis('off');h=text(double(signx), double(ymax),['+' num2str(ymax,3)]);  % text +ymax
    set(h,'FontSize',TICKFONTSIZE);         % choose font size
    set(h,'HorizontalAlignment','right','Clipping','off');

    ytick = g.ydir*(-ymax-0.3*ydiff);
    figure(curfig);tick = [int2str(xmin)]; h=text(double(xmin),double(ytick),tick);
    set(h,'FontSize',TICKFONTSIZE);         % choose font size
    set(h,'HorizontalAlignment','center',...
          'Clipping','off');  % center text

    tick = [xlabel]; figure(curfig);h=text(double(xmin+xdiff/2),double(ytick-0.5*g.ydir*ydiff),tick);
    set(h,'FontSize',TICKFONTSIZE);         % choose font size
    set(h,'HorizontalAlignment','center',...
          'Clipping','off');  % center text

    tick = [int2str(xmax)]; figure(curfig);h=text(double(xmax),double(ytick),tick);
    set(h,'FontSize',TICKFONTSIZE);         % choose font size
    set(h,'HorizontalAlignment','center',...
          'Clipping','off');  % center text

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
            'ylabel(''''Potential (\muV)'''');' ];
    axcopy(gcf, com); % turn on popup feature
                      %
                      %%%%%%%%%%%%%%%%%% Make printed figure fill page %%%%%%%%%%%%%%%%%%%%%%%%%%%
                      %
    orient tall
    % curfig = gcf;
    % h=figure(curfig);
    % set(h,'PaperPosition',[0.2 0.3 7.6 10]); % stretch out the plot on the page

