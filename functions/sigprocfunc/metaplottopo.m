% metaplottopo() - plot concatenated multichannel data epochs in a topographic or
%               rectangular array. Uses a channel location file with the same
%               format as topoplot(), or else plots data on a rectangular grid.
%
% Usage:
%    >> axes = metaplottopo(data, 'key1', 'val1', 'key2', 'val2')
%
% Inputs:
%   data       = data consisting of consecutive epochs of (chans,frames)
%                or (chans,frames,n)
%
% Optional inputs:
%  'chanlocs'  = [struct] channel structure or file plot ERPs at channel
%                locations. See help readlocs() for data channel format.
%  'geom'      = [rows cols] plot ERP in grid (overwrite previous option).
%                Grid size for rectangular matrix. Example: [6 4].
%  'title'     = [string] general plot title {def|'' -> none}
%  'chans'     = vector of channel numbers to plot {def|0 -> all}
%  'axsize'    = [x y] axis size {default [.07 .07]}
%  'plotfunc'  = [string] plot function name. If none is entered, axes
%                are created and returned.
%  'plotargs'  = [cell] plotting function arguments. 
%  'datapos'   = [integer] position of data array in the function call. 
%                Default is 1.
%
% Output:
%  Axes        = [real] array of axes handles of the same length as the
%                number of plotted channels.
%  Channames   = [cell] cell array of channel name for each plot
%
% Author: Arnaud Delorme, Scott Makeig, CERCO, CNRS, 2007-
%
% See also: plottopo()

% Copyright (C) 2007, Arnaud Delorme, CERCO, arno@sccn.ucsd.edu
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

function [Axes, outchannames ]= metaplottopo(data, varargin);

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
MAXCHANS      = 256;     % can be increased
curfig = gcf;            % learn the current graphic figure number
%
%%%%%%%%%%%%%%%%%%%% Default settings - use commandline to override %%%%%%%%%%%
%
DEFAULT_AXWIDTH  = 0.05; %
DEFAULT_AXHEIGHT = 0.08; %
DEFAULT_SIGN = -1;                        % Default - plot positive-up
ISRECT = 0;                               % default
ISSPEC = 0;                               % Default - not spectral data

outchannames = {};
if nargin < 1
    help metaplottopo
    return
end

if iscell(data), nchans = size(data{1},1);
else             nchans = size(data,1);
end;
g = finputcheck(varargin, { 'chanlocs'  ''    []          '';
    'chans'     'integer'               [1 size(data,1)]  [1:nchans];
    'geom'      'integer'               [1 Inf]     [];
    'title'     'string'                []          '';
    'plotfunc'  'string'                []          '';
    'plotargs'  'cell'                  []          {};
    'datapos'   'integer'               []          1;
    'calbar'    'real'                  []          [];
    'axcopycom' 'string'                []          '';
    'axsize'    'float'                 [0 1]       [nan nan]}, 'metaplottopo' );
if isstr(g), error(g); end;
if length(g.chans) == 1 & g.chans(1) ~= 0, error('can not plot a single ERP'); end;

[chans,framestotal]=size(data);           % data size

%
%%%%%%%%%%%%%%% Substitute defaults for missing parameters %%%%%%%%%%%%%%%%
%

if isempty(g.chanlocs) & isempty(g.geom)
    n = ceil(sqrt(length(g.chans)));
    g.geom = [n n];
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Test parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
axwidth  = g.axsize(1);
axheight = g.axsize(2);
if ~isempty(g.geom)
    if isnan(axheight) % if not specified
        axheight = gcapos(4)/(g.geom(1)+1);
        axwidth  = gcapos(3)/(g.geom(2)+1);
    end
else
    axheight = DEFAULT_AXHEIGHT*(gcapos(4)*1.25);
    axwidth =  DEFAULT_AXWIDTH*(gcapos(3)*1.3);
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

axcolor= get(0,'DefaultAxesXcolor'); % find what the default x-axis color is
vertcolor = 'k';
horicolor = vertcolor;

%
%%%%%%%%%%%%%%%%%%%%%%%%% Plot and label specified channels %%%%%%%%%%%%%%%%%%
%

if ~isempty(data)
    if iscell(data)
        data{1} = data{1}(g.chans,:);
        data{2} = data{2}(g.chans,:);
    else
        data = data(g.chans,:);
    end;
end;

%
%%%%%%%%%%%%%%%%%%%%%%%%%% Print plot info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
figure(curfig); h=gca;title(g.title); % title plot
hold on
axis('off')

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Read chan_locs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if isempty(g.chanlocs) % plot in a rectangular grid
    ISRECT = 1;
    ht = g.geom(1);
    wd = g.geom(2);
    if length(g.chans) > ht*wd
        fprintf('metaplottopo(): (%d) channels to be plotted > grid size [%d %d]\n',...
            length(g.chans),ht,wd);
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
        error('metaplottopo(): data channels must be <= ''chanlocs'' channels')
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
    xvals = (xvals-mean([max(xvals) min(xvals)]))/(max(xvals)-min(xvals)); % recenter
    xvals = gcapos(1)+gcapos(3)/2+PLOT_WIDTH*xvals;   % controls width of plot
    % array on current axes
end;
yvals = gcapos(2)+gcapos(4)/2+PLOT_HEIGHT*yvals;  % controls height of plot
% array on current axes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

Axes = [];
fprintf('Plotting all channel...');
for c=1:length(g.chans), %%%%%%%% for each data channel %%%%%%%%%%%%%%%%%%%%%%%%%%

    xcenter = xvals(c); if isnan(xcenter), xcenter = 0.5; end; 
    ycenter = yvals(c); if isnan(ycenter), ycenter = 0.5; end;
    Axes = [Axes axes('Units','Normal','Position', ...
        [xcenter-axwidth/2 ycenter-axheight/2 axwidth axheight])];
    hold on;
    axis('off');
    %axes('Units','Normal','Position', [xcenter-axwidth/2 ycenter-axheight/2 axwidth axheight])
    %axis('off');

    %
    %%%%%%%%%%%%%%%%%%%%%%% Plot data traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if ~isempty( g.plotfunc )
        %figure(curfig);
        eval( [ 'func = @' g.plotfunc ';' ] );
        if iscell(data), tmp = { g.plotargs{1:g.datapos(1)-1} data{1}(c,:) g.plotargs{g.datapos(1):g.datapos(2)-1} data{2}(c,:) g.plotargs{g.datapos(2):end}};
        else             tmp = { g.plotargs{1:g.datapos-1}    data(c,:)    g.plotargs{g.datapos:end} };
        end;
        tmp = { tmp{:} 'title' channames(c,:) 'plotmode' 'topo'};
        feval(func, tmp{:});
    end;
    outchannames{c} = deblank(channames(c,:));
end; % c, chans / subplot

fprintf('\n');

%
%%%%%%%%%%%%%%%%%%%%% Make time and amp cal bar %%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~isempty(g.calbar)
    ax = axes('Units','Normal','Position', [0.85 0.1 axwidth axheight]); % FIX!!!!
    axes(ax)
    axis('off');
    if   g.calbar(3) < g.calbar(4), g.ydir = 1;
    else g.calbar(5) = g.calbar(3); g.calbar(3) = []; g.ydir = -1;
    end;
    [xmin xmax ymin ymax] = deal(g.calbar(1), g.calbar(2),g.calbar(3), g.calbar(4));
    figure(curfig);p=plot([0 0],[ymin ymax],'color','k'); % draw vert axis at zero
    if g.ydir == -1
        set(gca, 'ydir', 'reverse');
    end;
    hold on
    figure(curfig);p=plot([xmin xmax],[0 0],'color','k'); % draw horizontal axis
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    axis off;
    
    % draw text limits
    % ----------------
    xdiff = xmax - xmin;
    ydiff = ymax - ymin;
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

    tick = 'Time (ms)'; figure(curfig);h=text(double(xmin+xdiff/2),double(ytick-0.5*g.ydir*ydiff),tick);
    set(h,'FontSize',TICKFONTSIZE);         % choose font size
    set(h,'HorizontalAlignment','center',...
                        'Clipping','off');  % center text

    tick = [int2str(xmax)]; figure(curfig);h=text(double(xmax),double(ytick),tick);
    set(h,'FontSize',TICKFONTSIZE);         % choose font size
    set(h,'HorizontalAlignment','center',...
                        'Clipping','off');  % center text
    
end;

%        'set(gcbf, ''''unit'''', ''''pixel'''');' ...
%        'tmp = get(gcbf, ''''position'''');' ...
%        'tmp2 = get(0, ''''screensize'''');' ...
%        'if tmp(2)+500 > tmp2(4), tmp(2) = tmp2(4)-500; end;' ...
%        'set(gcbf, ''''position'''', [ tmp(1) tmp(2) 560   420]);' ...
g.axcopycom = [ 'axis on;' ...
        'tmp = get(gca, ''''userdata'''');' ...
        'if ~isempty(tmp), xlabel(tmp{1});' ...
        'ylabel(tmp{2});' ...
        'legend(tmp{3}{:}); end; clear tmp tmp2;' ];
axcopy(gcf, g.axcopycom); % turn on popup feature
icadefs;
set(gca,'Color',BACKCOLOR);               % set the background color
