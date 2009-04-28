function Fig = eeg_plot(time,erp,scale,interval,width,titles)

% eeg_plot - create a report style plot of ERP data
%
% Usage: eeg_plot(time,erp,scale,interval,width,titles)
% 
% time  -  an Mx1 column vector of times (in msec) for each row of erp
% 
% erp   -  an MxN matrix, M rows of time points, N channels
% 
% scale -  a 1x4 row vector = [xmin,xmax,ymin,ymax], the
%          defaults are: min(time), max(time)
%                        -/+ max(max(abs(erp)))
% 
% interval - a 1x2 row vector = [xint, yint], to specify the
%            x,y tick mark intervals (defaults are 100 msec and 1 uV)
% 
% width -  1 = journal page width (default, 15.8 cm)
%          0 = journal column width (7.4 cm)
% 
%          plot height is (2/3) * width
% 
% titles   - a text struct with (default) fields:
% 
%            titles.title  - figure title ('')
%            titles.xlabel - x label ('msec')
%            titles.xalign - text alignment ('center')
%            titles.xxpos  - position of x label relative to x axis ([])
%            titles.xypos  - position of x label relative to y axis ([])
%            titles.ylabel - y label ('\muV', tex for uV)
%            titles.yalign - text alignment ('center')
%            titles.yypos  - position of y label relative to y axis ([])
%            titles.yxpos  - position of y label relative to x axis ([])
% 
% The figure may appear small on screen, but it is formatted to
% print in either the width of a journal page (15.8 cm) or the
% width of a journal column (7.5 cm).
% To adjust it to your preferred dimensions, hack the code, which is
% well documented.  If you want different colours or all black lines,
% hack the code at the plot command. Alternatively, use
% a graphics package to convert specific colours or convert all
% colours to greyscale.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $

% Licence:  GNU GPL, no express or implied warranties
% History:  07/2000, Darren.Weber_at_radiology.ucsf.edu
%                    Adapted from (c) Drea Thomas, June 1995
%                                  MathWorks newsletter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Requirements:
% Grids 
% Legend on top of diagram 

if ~exist('interval','var'), interval = []; end



% ----- Set the figure dimensions and scaling -----

if ~exist('width','var'),
    % Width of a journal page = 15.8 cm
    width = 15.8;
    fontsize = 8;
    ylabeloffset = 0.1;
elseif isempty(width),
    % Width of a journal page = 15.8 cm
    width = 15.8;
    fontsize = 8;
    ylabeloffset = 0.1;
elseif width,
    % Width of a journal page = 15.8 cm
    width = 15.8;
    fontsize = 8;
    ylabeloffset = 0.1;
else
    % Width of a journal column = 7.4 cm
    width = 7.4;
    fontsize = 6;
    ylabeloffset = 0.1;
end

height = (2/3) * width;

% If you want to control exactly how a plot is going to look on paper, 
% it is a good idea to WYSIWYG (what you see is what you get) the figure 
% window. That is, set the figure window 'Units' property to physical 
% units, and set the figure window size the same as the 'PaperPosition' 
% size. For our case,
Fig = figure;
set(Fig,'PaperUnits','centimeters','PaperPosition',[1 1 width height]) 
set(Fig,'Units','centimeters','position',get(Fig,'PaperPosition'))
movegui(Fig,'center');
% Note that the lower left-hand corner of the graph is offset by 1 cm 
% in X and Y because most printers can't physically print to the edge of 
% a piece of paper. If you specified a [0 0] offset, the plot would be 
% cut off.



% ----- Set default fonts -----

% Journals require a letter size no smaller than approx. 2 mm after
% reduction of the figure.
% Font sizes are usually measured in points. Since there are about 
% 2.8 points per mm, the minimum fontsize should be 6 or 7.  Let's 
% set the default fontsize for the figure window to 8.
set(Fig,'DefaultAxesFontSize',fontsize)
set(Fig,'DefaultAxesFontName','Arial')
% See Drea's Desk, April 1994 for a discussion of default properties
% (ftp://ftp.mathworks.com/pub/doc/tmw-digest/apr94).



% ----- plot the data -----

%plot(time,erp,'k'); % black lines

if size(erp,2) == 1,
    plot(time,erp(:,1),'r-', 'linewidth',1.0);
elseif size(erp,2) == 2,
    plot(time,erp(:,1),'r-', 'linewidth',1.0); hold on;
    plot(time,erp(:,2),'b--','linewidth',0.8);
else
    plot(time,erp,'linewidth',0.8); % msec
    %plot(time./1000,erp,'linewidth',0.5); % convert to sec
end



% ----- Set the axis properties -----

make_scale = 0;
if ~exist('scale','var'), make_scale = 1; end
if isempty(scale), make_scale = 1; end
if make_scale,
    %xmin = min(time) / 1000; % convert to sec
    %xmax = max(time) / 1000;
    
    xmin = min(time); % msec
    xmax = max(time);
    
    elecmax = max(abs(erp));
    ymax = max(elecmax);
    ymin = -1 * ymax;
    
    scale = [xmin,xmax,ymin,ymax];
else
    xmin = scale(1);
    xmax = scale(2);
    ymin = scale(3);
    ymax = scale(4);
end
if ~exist('interval','var'),
    interval(1) = 100; % msec
    interval(2) = (ymax - ymin) / 10;
    %interval = [0.1 1]; % sec
end

axis(scale);

Ax = gca;
set(Ax,'YDir','reverse');
set(Ax,'Box','off');

set(Ax,'XGrid','on');
if ~isempty(interval),
    set(Ax,'XTickLabelMode','manual');
    set(Ax,'XTickMode','manual');
    set(Ax,'XTick',     [xmin:interval(1):xmax]);
    set(Ax,'XTickLabel',[xmin:interval(1):xmax]);
end
set(Ax,'YGrid','on');
if ~isempty(interval),
    set(Ax,'YTickLabelMode','manual');
    set(Ax,'YTickMode','manual');
    set(Ax,'YTick'     ,[ymin:interval(2):ymax]);
    set(Ax,'YTickLabel',[ymin:interval(2):ymax]);
end


% ----- Reposition the axis -----

% The default axis size is too large (the labels get cut off). 
% Hence, modify the size and location of the axis within the 
% figure window,
pos = get(Ax,'position');      % This is in normalized coordinates
pos(1) = pos(1)+.01;           % Shift axes left by a factor of .01
pos(2) = pos(2)+.02;           % Shift axes up by a factor of .02
pos(3) = pos(3) + pos(3)*.05;  % Stretch x axis by a factor of .05
set(Ax,'position',pos);



% ----- Add titles -----

make_title = 0;
if ~exist('titles','var'), make_title = 1; end
if isempty(titles), make_title = 1; end
if make_title,
	titles.title  = '';
	titles.xlabel = 'msec';
	titles.xalign = 'center';
	titles.xxpos  = [];
	titles.xypos  = [];
	titles.ylabel = '\muV';
	titles.yalign = 'center';
    titles.yrot   = 0;
	titles.yypos  = [];
	titles.yxpos  = [];
end

if ~isfield(titles,'title'),  titles.title  = ''; end
if ~isfield(titles,'xlabel'), titles.xlabel = 'msec'; end
if ~isfield(titles,'xalign'), titles.xalign = 'center'; end
if ~isfield(titles,'xxpos'),  titles.xxpos  = []; end
if ~isfield(titles,'xypos'),  titles.xypos  = []; end

if ~isfield(titles,'ylabel'), titles.ylabel = '\muV'; end
if ~isfield(titles,'yalign'), titles.yalign = 'center'; end
if ~isfield(titles,'yrot'),   titles.yrot = 0; end
if ~isfield(titles,'yypos'),  titles.yypos  = []; end
if ~isfield(titles,'yxpos'),  titles.yxpos  = []; end

title(titles.title);
Hx = xlabel(titles.xlabel,'Rotation',0,'HorizontalAlignment',titles.xalign);
Hy = ylabel(titles.ylabel,'Rotation',titles.yrot,'HorizontalAlignment',titles.yalign);
Yextent = get(Hy,'extent');
% Yextent is the position and size of text. A four-element read-only 
% vector that defines the size and position of the text string.
% [left,bottom,width,height]
% If the Units property is set to data (the default), left and bottom 
% are the x and y coordinates of the lower-left corner of the text 
% Extent rectangle.  For all other values of Units, left and bottom 
% are the distance from the lower-left corner of the axes position 
% rectangle to the lower-left corner of the text Extent rectangle. 
% width and height are the dimensions of the Extent rectangle. All 
% measurements are in units specified by the Units property.
if titles.xxpos,
    xxoffset = titles.xxpos;
else
    pos = get(Hx,'Position');
    xxoffset = pos(1);
end
if titles.xypos,
    xyoffset = titles.xypos;
else
    pos = get(Hx,'Position');
    xyoffset = pos(2);
end
set(Hx,'Position',[xxoffset xyoffset]);

if titles.yypos,  yyoffset = titles.yypos;
else              yyoffset = Yextent(4) / 4;
end
if titles.yxpos,  yxoffset = titles.yxpos;
else
    pos = get(Hy,'Position');
    yxoffset = pos(1) + ( pos(1) * 0.1 );
end
set(Hy,'Position',[yxoffset yyoffset]);



% ----- Now for the legend -----

% h = legend('Sin','Cos'); 
% axes(h); 
% refresh;

% Legend is an axis. Unless it is the current axes, it is 
% drawn "behind" the other axis (i.e., your plot). So, just 
% before you print, make it the current axis and it will 
% appear "on top" (you wouldn't want to always make it the 
% current axis because then the next PLOT command would go 
% into the legend). 
