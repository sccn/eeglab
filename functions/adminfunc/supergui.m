% supergui() - a comprehensive gui automatic builder. This function help
%              to create GUI very fast without bothering about the 
%              positions of the elements. After creating a geometry, 
%              elements just place themselves into the predefined 
%              locations. It is especially usefull for figure where you
%              intend to put text button and descriptions.
%
% Usage:
%   >> [handlers, width, height ] = ...
%             supergui( 'key1', 'val1', 'key2', 'val2', ... );
% 
% Inputs:
%   'fig'       - figure handler, if not given, create a new figure.
%   'geom'      - cell array of cell array of integer vector. Each cell
%               array defines the coordinate of a given input in the following
%               manner: { nb_row nb_col [x_topcorner y_topcorner]
%               [x_bottomcorner y_bottomcorner] };
%   'geomhoriz' - integer vector or cell array of numerical vectors describing the 
%               geometry of the elements in the figure. 
%               - if integer vector, vector length is the number of rows and vector 
%               values are the number of 'uilist' elements in each row.
%               For example, [2 3 2] means that the
%               figures will have 3 rows, with 2 elements in the first
%               and last row and 3 elements in the second row.
%               - if cell array, each vector describes the relative widths
%               of items in each row. For example, { [2 8] [1 2 3] } which means
%               that figures will have 2 rows, the first one with 2
%               elements of relative width 2 and 8 (20% and 80%). The
%               second row will have 3 elements of relative size 1, 2 
%               and 3 (1/6 2/6 and 3/6).
%   'geomvert' - describting geometry for the rows. For instance
%               [1 2 1] means that the second row will be twice the height
%               of the other ones. If [], all the lines have the same height.
%   'uilist'   - list of uicontrol lists describing elements properties
%               { { ui1 }, { ui2 }... }, { 'uiX' } being GUI matlab 
%               uicontrol arguments such as { 'style', 'radiobutton', 
%               'String', 'hello' }. See Matlab function uicontrol() for details.
%   'borders'  - [left right top bottom] GUI internal borders in normalized
%               units (0 to 1). Default values are 
%   'title'    - optional figure title
%   'userdata' - optional userdata input for the figure
%   'inseth'   - horizontal space between elements. Default is 2% 
%                of window size.
%   'insetv'   - vertical space between elements. Default is 2% 
%                of window height.
%   'spacing'  - [horiz vert] spacing in normalized units. Default 
%   'spacingtype' - ['absolute'|'proportional'] abolute means that the 
%                spacing values are fixed. Proportional means that they
%                depend on the number of element in a line.
%   'minwidth' - [integer] minimal width in pixels. Default is none.
%
% Hint:
%    use 'print -mfile filemane' to save a matlab file of the figure.
%
% Output:
%    handlers  - all the handler of the elements (in the same order as the
%                uilist input).
%    height    - adviced widht for the figure (so the text look nice).   
%    height    - adviced height for the figure (so the text look nice).   
%
% Example:
%    figure;   
%    supergui( 'geomhoriz', { 1 1 }, 'uilist', { ...
%        { 'style', 'radiobutton', 'string', 'radio' }, ...
%        { 'style', 'pushbutton' , 'string', 'push' } } );
%      
% Author: Arnaud Delorme, CNL / Salk Institute, La Jolla, 2001-
%
% See also: eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.59  2009/12/16 02:02:31  arno
% spacing for Mac`
%
% Revision 1.58  2009/11/11 03:27:18  arno
% checkbox color under Linux
%
% Revision 1.57  2009/11/11 00:28:52  arno
% New GUI format
%
% Revision 1.56  2009/04/21 19:24:36  arno
% Fix Mac GUI ratio problem
%
% Revision 1.55  2007/01/08 05:21:00  toby
% help edit
%
% Revision 1.54  2007/01/06 08:48:27  toby
% Help clarified concerning 'geomhoriz'
%
% Revision 1.53  2006/10/18 02:34:21  toby
% spelling
%
% Revision 1.52  2006/09/12 16:43:24  arno
% allow vertshift
%
% Revision 1.51  2006/02/23 22:32:43  arno
% nothing
%
% Revision 1.50  2006/02/10 22:33:13  arno
% nothing
%
% Revision 1.49  2006/02/09 23:12:35  arno
% nothing
%
% Revision 1.48  2006/01/27 21:32:02  arno
% hieght for pop-upmenus and buttons
%
% Revision 1.47  2006/01/24 19:40:23  arno
% nothing
%
% Revision 1.46  2006/01/19 19:42:40  arno
% color of popupmenus
%
% Revision 1.45  2005/11/09 23:23:05  arno
% nothing
%
% Revision 1.44  2005/11/09 23:06:29  arno
% header
%
% Revision 1.43  2005/11/09 22:58:33  arno
% changing inset default; rewrote function input
%
% Revision 1.42  2005/09/27 21:55:53  arno
% change multiplicative factor
%
% Revision 1.41  2005/09/27 21:55:09  arno
% remove *1.15 aspect ratio multuplication for Windows
%
% Revision 1.40  2004/11/10 16:51:13  arno
% debug last + header
%
% Revision 1.39  2004/11/10 16:47:37  arno
% do not take listbox size into account
%
% Revision 1.38  2004/11/05 19:23:04  arno
% same
%
% Revision 1.37  2004/11/05 19:20:41  arno
% pushbutton case
%
% Revision 1.36  2003/02/21 16:50:34  arno
% create uicontrol in current figure
%
% Revision 1.35  2002/11/13 00:54:49  arno
% replace gcf by fig
%
% Revision 1.34  2002/10/23 15:06:40  arno
% isppc -> computer
%
% Revision 1.33  2002/10/15 16:25:15  arno
% magnify edit boxes windows
%
% Revision 1.32  2002/10/15 14:42:15  arno
% button text aspect
%
% Revision 1.31  2002/08/20 22:33:34  arno
% debug for Mac
%
% Revision 1.30  2002/08/19 19:12:16  arno
% debugging last
%
% Revision 1.29  2002/08/19 19:10:13  arno
% figure bigger for MAC
%
% Revision 1.28  2002/08/17 02:38:42  arno
% debugging
%
% Revision 1.27  2002/08/14 18:17:23  arno
% new supergui arg
%
% Revision 1.26  2002/08/14 18:07:20  arno
% changing default checkbox backcolor
%
% Revision 1.25  2002/08/14 16:32:37  arno
% yfact takes into account button now
%
% Revision 1.24  2002/08/13 23:04:51  arno
% debug pop_merge
%
% Revision 1.23  2002/08/13 18:59:49  arno
% update automatic INSETY
%
% Revision 1.22  2002/08/13 18:56:01  arno
% adjustments
%
% Revision 1.21  2002/08/13 18:50:54  arno
% same
%
% Revision 1.20  2002/08/13 18:50:31  arno
% vert align
%
% Revision 1.19  2002/08/13 17:40:41  arno
% same
%
% Revision 1.18  2002/08/13 17:38:13  arno
% ordinate adjustment
%
% Revision 1.17  2002/08/13 17:28:56  arno
% new geometry
%
% Revision 1.16  2002/08/13 00:20:42  arno
% same
%
% Revision 1.15  2002/08/13 00:20:22  arno
% update position
%
% Revision 1.14  2002/08/12 21:42:59  arno
% ignore pushbutton extent
%
% Revision 1.13  2002/08/12 16:00:57  arno
% same
%
% Revision 1.12  2002/08/12 16:00:13  arno
% do not adapt size for edit windows
%
% Revision 1.11  2002/08/12 15:57:04  arno
% size calculation
%
% Revision 1.10  2002/08/12 14:50:15  arno
% color
%
% Revision 1.9  2002/08/12 14:47:40  arno
% color
%
% Revision 1.8  2002/08/12 14:30:55  arno
% background
%
% Revision 1.7  2002/08/12 01:17:42  arno
% update colors
%
% Revision 1.6  2002/08/12 00:41:41  arno
% updating default colors
%
% Revision 1.5  2002/07/18 17:18:31  arno
% offscreen correction
%
% Revision 1.4  2002/07/18 17:13:05  arno
% same
%
% Revision 1.3  2002/07/18 17:11:19  arno
% correct out-of screen problem
%
% Revision 1.2  2002/07/18 17:07:40  arno
% no modif
%
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%

function [handlers, outheight, allhandlers] = supergui( varargin);

% handlers cell format
% allhandlers linear format
handlers = {};
outheight = 0;

if nargin < 2
	help supergui;
	return;
end;

% decoding input and backward compatibility
% -----------------------------------------
if isstr(varargin{1})
    options = varargin;
else
    options = { 'fig'      varargin{1} 'geomhoriz' varargin{2} ...
                'geomvert' varargin{3} 'uilist'    varargin(4:end) }; 
end
g = finputcheck(options, { 'geomhoriz' 'cell'   []      {};
                           'fig'       'real'   []      0;
                           'geom'      'cell'   []      {};
                           'uilist'    'cell'   []      {};
                           'title'     'string' []      '';
                           'userdata'  ''       []      [];
                           'geomvert'  'real'   []      [];
                           'horizontalalignment'  'string'   { 'left' 'right' 'center' } 'left';
                           'minwidth'  'real'   []      10;
                           'borders'   'real'   []      [0.05 0.04 0.07 0.06];
                           'spacing'   'real'   []      [0.02 0.01];
                           'inseth'    'real'   []      0.02; % x border absolute (5% of width)
                           'insetv'    'real'   []      0.02 }, 'supergui');
if isstr(g), error(g); end
if ~isempty(g.geomhoriz)
    maxcount = sum(cellfun('length', g.geomhoriz));
    if maxcount ~= length(g.uilist)
        warning('Wrong size for ''geomhoriz'' input');
    end;
    if ~isempty(g.geomvert)
        if length(g.geomvert) ~= length(g.geomhoriz)
            warning('Wrong size for ''geomvert'' input');
        end;
    end;
    g.insetv = g.insetv/length(g.geomhoriz);
end;
if ~isempty(g.geom)
    if length(g.geom) ~= length(g.uilist)
        warning('Wrong size for ''geom'' input');
    end;
    maxcount = length(g.geom);
end;

% create new figure
% -----------------
if g.fig == 0
	g.fig = figure('visible','off');
end

% converting the geometry formats
% -------------------------------
if ~isempty(g.geomhoriz) & ~iscell( g.geomhoriz )
	oldgeom = g.geomhoriz;
	g.geomhoriz = {};
	for row = 1:length(oldgeom)
		g.geomhoriz = { g.geomhoriz{:} ones(1, oldgeom(row)) };
	end;
end
if isempty(g.geomvert)
	g.geomvert = ones(1, length(g.geomhoriz));
end

% converting to the new format
% ----------------------------
if isempty(g.geom)
    count = 1;
    incy  = 0;
    sumvert  = sum(g.geomvert);
    maxhoriz = 1;
    for row = 1:length(g.geomhoriz)
        incx = 0;
        maxhoriz = max(maxhoriz, length(g.geomhoriz{row}));
        ratio = length(g.geomhoriz{row})/sum(g.geomhoriz{row});
        for column = 1:length(g.geomhoriz{row})
            g.geom{count} = { length(g.geomhoriz{row}) sumvert [incx incy] [g.geomhoriz{row}(column)*ratio g.geomvert(row)] };
            incx = incx+g.geomhoriz{row}(column)*ratio;
            count = count+1;
        end;
        incy = incy+g.geomvert(row);
    end;
    g.borders(1:2) = g.borders(1:2)/maxhoriz*5;
    g.borders(3:4) = g.borders(3:4)/sumvert*10;
    g.spacing(1)   = g.spacing(1)/maxhoriz*5;
    g.spacing(2)   = g.spacing(2)/sumvert*10;
end;

% get axis coordinates
% --------------------
set(g.fig, 'menubar', 'none', 'numbertitle', 'off');		
pos = [0 0 1 1]; % plot relative to current axes
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)]; % allow to use normalized position [0 100] for x and y
axis('off');

% creating guis
% -------------
row    = 1; % count the elements
column = 1; % count the elements
factmultx = 0;
factmulty = 0; %zeros(length(g.geomhoriz));
for counter = 1:maxcount

	% init
	clear rowhandle;
    gm = g.geom{counter};
    [posx posy width height] = getcoord(gm{1}, gm{2}, gm{3}, gm{4}, g.borders, g.spacing);
    
    try
        currentelem = g.uilist{ counter };
    catch
        fprintf('Warning: not all boxes were filled\n');
        return;
    end;		
    if ~isempty(currentelem)
        % decode metadata
        % ---------------
        if strcmpi(currentelem{1}, 'link2lines'), 
            currentelem(1) = []; 
            allhandlers(counter) = uicontrol(g.fig, 'unit', 'normalized', 'position', ...
                                          [posx-width/2 posy+3.6*height/2 width/2 0.005].*s+q, 'style', 'pushbutton', 'string', '');
            allhandlers(counter) = uicontrol(g.fig, 'unit', 'normalized', 'position', ...
                                          [posx-width/2 posy+0.7*height/2 width/2 0.005].*s+q, 'style', 'pushbutton', 'string', '');
            allhandlers(counter) = uicontrol(g.fig, 'unit', 'normalized', 'position', ...
                                          [posx posy+0.7*height/2 0.005 3/2*height].*s+q, 'style', 'pushbutton', 'string', '');
            allhandlers(counter) = uicontrol(g.fig, 'unit', 'normalized', 'position', ...
                                          [posx posy+2.1*height/2 width/2 0.005].*s+q, 'style', 'pushbutton', 'string', '');
            allhandlers(counter) = 0;
        else
            if strcmpi(currentelem{1}, 'width'),
                 curwidth = currentelem{2};
                 currentelem(1:2) = [];
            else curwidth = 0;
            end;
            if strcmpi(currentelem{1}, 'align'),
                 align = currentelem{2};
                 currentelem(1:2) = [];
            else align = 'right';
            end;
            if strcmpi(currentelem{1}, 'stickto'),
                 stickto = currentelem{2};
                 currentelem(1:2) = [];
            else stickto = 'none';
            end;
            if strcmpi(currentelem{1}, 'vertshift'), currentelem(1) = []; addvert = -height/2; 
            else                                                          addvert = 0;   
            end;
            allhandlers(counter) = uicontrol(g.fig, 'unit', 'normalized', 'position', ...
                                          [posx posy+addvert width height].*s+q, currentelem{:});

            % this simply compute a factor so that all uicontrol will be visible
            % ------------------------------------------------------------------
            style = get( allhandlers(counter), 'style');			
            set( allhandlers(counter), 'units', 'pixels');		
            curpos = get(allhandlers(counter), 'position');
            curext = get(allhandlers(counter), 'extent');
            if curwidth ~= 0
                curwidth = curwidth/((factmultx-1)/1.85+1);
                if strcmpi(align, 'right')
                    curpos(1) = curpos(1)+curpos(3)-curwidth;
                elseif strcmpi(align, 'center')
                    curpos(1) = curpos(1)+curpos(3)/2-curwidth/2;
                end;
                set(allhandlers(counter), 'position', [ curpos(1) curpos(2) curwidth curpos(4) ]);
                if strcmpi(stickto, 'on')
                    set( allhandlers(counter-1), 'units', 'pixels');
                    curpos2 = get(allhandlers(counter-1), 'position');
                    set(allhandlers(counter-1), 'position', [ curpos(1)-curpos2(3)-10 curpos2(2) curpos2(3) curpos2(4) ]);
                    set( allhandlers(counter-1), 'units', 'normalized');
                end;
                curext(3) = curwidth;
            end;
            set( allhandlers(counter), 'units', 'normalized');			

            if ~strcmp(style, 'edit') & ~strcmp(style, 'pushbutton')
                %tmp = curext(3)/curpos(3);
                %if tmp > 3*factmultx && factmultx > 0, adsfasd; end;
                factmultx = max(factmultx, curext(3)/curpos(3));
            end;
            if  ~strcmp(style, 'listbox')
                factmulty = max(factmulty, curext(4)/curpos(4));
            end;

            % Uniformize button text aspect (first letter must be upercase)
            % -----------------------------
            if strcmp(style, 'pushbutton')
                tmptext = get(allhandlers(counter), 'string');
                if length(tmptext) > 1
                    if upper(tmptext(1)) ~= tmptext(1) | lower(tmptext(2)) ~= tmptext(2)
                        tmptext = lower(tmptext);
                        try, tmptext(1) = upper(tmptext(1)); catch, end;
                    end;
                end;
                set(allhandlers(counter), 'string', tmptext);
            end;
        end;
    else 
        allhandlers(counter) = 0;
    end;
end;

% adjustments
% -----------
factmultx = factmultx*1.02;% because some text was still hidden
%factmultx = factmultx*1.2; 
if factmultx < 0.1
	factmultx = 0.1;
end;

% for MAC (magnify figures that have edit fields)
% -------
warning off;
try, 
	if strcmpi(computer, 'MAC') | strcmpi(computer, 'MACI')
        factmulty = factmulty*1.5;
	elseif ~isunix % windows
        factmulty = factmulty*1.08;
    end;
catch, end;
factmulty = factmulty*0.9; % global shinking
warning on;	

% scale and replace the figure in the screen
% -----------------------------------------
pos = get(g.fig, 'position');
if factmulty > 1
	pos(2) = max(0,pos(2)+pos(4)-pos(4)*factmulty);
end;
pos(1) = pos(1)+pos(3)*(1-factmultx)/2;
pos(3) = max(pos(3)*factmultx, g.minwidth);
pos(4) = pos(4)*factmulty;
set(g.fig, 'position', pos);

% vertical alignment to bottom for text
% ---------------------------------------
for index = 1:length(allhandlers)
	if allhandlers(index) ~= 0
		if strcmp(get(allhandlers(index), 'style'), 'text')
            set(allhandlers(index), 'unit', 'pixel');
			curpos = get(allhandlers(index), 'position');
			curext = get(allhandlers(index), 'extent');
			set(allhandlers(index), 'position', [curpos(1) curpos(2)-4 curpos(3) curext(4)]);
            set(allhandlers(index), 'unit', 'normalized');
		end;
	end;
end;

% setting defaults colors
%------------------------
try, icadefs;
catch,
	GUIBACKCOLOR        =  [.8 .8 .8];     
	GUIPOPBUTTONCOLOR   = [.8 .8 .8];    
	GUITEXTCOLOR        = [0 0 0];
end;

hh = findobj(allhandlers, 'parent', g.fig, 'style', 'text');
%set(hh, 'BackgroundColor', get(g.fig, 'color'), 'horizontalalignment', 'left');
set(hh, 'Backgroundcolor', GUIBACKCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
set(g.fig, 'color',GUIBACKCOLOR );
set(hh, 'horizontalalignment', g.horizontalalignment);

hh = findobj(allhandlers, 'style', 'edit');
set(hh, 'BackgroundColor', [1 1 1]); %, 'horizontalalignment', 'right');

hh =findobj(allhandlers, 'parent', g.fig, 'style', 'pushbutton');
if ~strcmpi(computer, 'MAC') & ~strcmpi(computer, 'MACI') % this puts the wrong background on macs
    set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
    set(hh, 'foregroundcolor', GUITEXTCOLOR);
end;
hh =findobj(allhandlers, 'parent', g.fig, 'style', 'popupmenu');
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'parent', g.fig, 'style', 'checkbox');
set(hh, 'backgroundcolor', GUIBACKCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'parent', g.fig, 'style', 'listbox');
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'parent', g.fig, 'style', 'radio');
set(hh, 'foregroundcolor', GUITEXTCOLOR);
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(g.fig, 'visible', 'on');

% set userdata and title
% ----------------------
if ~isempty(g.userdata), set(g.fig, 'userdata', g.userdata); end;
if ~isempty(g.title   ), set(g.fig, 'name',     g.title   ); end;

return;

function [posx posy width height] = getcoord(geom1, geom2, coord1, sz, borders, spacing);
    
    coord2 = coord1+sz;
    borders(1:2) = borders(1:2)-spacing(1);
    borders(3:4) = borders(3:4)-spacing(2);
    
    % absolute positions
    posx   = coord1(1)/geom1;
    posy   = coord1(2)/geom2;
    posx2  = coord2(1)/geom1;
    posy2  = coord2(2)/geom2;
    width  = posx2-posx;
    height = posy2-posy;

    % add spacing
    posx   = posx+spacing(1)/2;
    width  = max(posx2-posx-spacing(1), 0.001);
    height = max(posy2-posy-spacing(2), 0.001);
    posy   = max(0, 1-posy2)+spacing(2)/2;
    
    % add border
    posx   = posx*(1-borders(1)-borders(2))+borders(1);
    posy   = posy*(1-borders(3)-borders(4))+borders(4);
    width  = width*( 1-borders(1)-borders(2));
    height = height*(1-borders(3)-borders(4));
    
function [posx posy width height] = getcoordold(geom1, geom2, coord1, sz);
    
    coord2 = coord1+sz;
    horiz_space  = 0.05/geom1;
    vert_space   = 0.05/geom2;
    horiz_border = min(0.1, 1/geom1)-horiz_space;
    vert_border  = min(0.2, 1.5/geom2)-vert_space;
    
    % absolute positions
    posx   = coord1(1)/geom1;
    posy   = coord1(2)/geom2;
    posx2  = coord2(1)/geom1;
    posy2  = coord2(2)/geom2;
    width  = posx2-posx;
    height = posy2-posy;

    % add spacing
    posx   = posx+horiz_space/2;
    width  = max(posx2-posx-horiz_space, 0.001);
    height = max(posy2-posy- vert_space, 0.001);
    posy   = max(0, 1-posy2)+vert_space/2;
    
    % add border
    posx   = posx*(1-horiz_border)+horiz_border/2;
    posy   = posy*(1- vert_border)+vert_border/2;
    width  = width*(1-horiz_border);
    height = height*(1-vert_border);
 
    %     posx   = coord1(1)/geom1+horiz_border*1/geom1/2;
%     posy   = 1-(coord1(2)/geom2+vert_border*1/geom2/2)-1/geom2;
%     
%     posx2  = coord2(1)/geom1+horiz_border*1/geom1/2;
%     posy2  = 1-(coord2(2)/geom2+vert_border*1/geom2/2)-1/geom2;
%     
%     width  = posx2-posx;
%     height = posy-posy2;
    
    %h = axes('unit', 'normalized', 'position', [ posx posy width height ]);
    %h = axes('unit', 'normalized', 'position', [ coordx/geom1 1-coordy/geom2-1/geom2 1/geom1 1/geom2 ]);