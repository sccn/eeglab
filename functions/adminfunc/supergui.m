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
%   'geomhoriz' - cell array describing the geometry of the elements
%               in the figure. For instance, [2 3 2] means that the
%               figures will have 3 rows, with 2 elements in the first
%               and last row and 3 elements in the second row.
%               An other syntax is { [2 8] [1 2 3] } which means
%               that figures will have 2 rows, the first one with 2
%               elements of relative width 2 and 8 (20% and 80%). The
%               second row will have 3 elements of relative size 1, 2 
%               and 3.
%   'geomvert' - describting geometry for the rows. For instance
%               [1 2 1] means that the second row will be twice the height
%               of the other ones. If [], all the lines have the same height.
%   'uilist'   - list of uicontrol lists describing elements properties
%               { { ui1 }, { ui2 }... }, { 'uiX' } being GUI matlab 
%               uicontrol arguments such as { 'style', 'radiobutton', 
%               'String', 'hello' }. See supergui() for details.
%   'title'    - optional figure title
%   'userdata' - optional userdata input for the figure
%   'inseth'   - horizontal space between elements. Default is 2% 
%                of window size.
%   'insetv'   - vertical space between elements. Default is 2% 
%                of window height.
%
% Hint:
%    use 'print -mfile filemane' to save a matlab file of the figure.
%
% Output:
%    handlers  - all the handler of the elements (in the same ordre as the
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
end;
g = finputcheck(options, { 'geomhoriz' 'cell'   []      [];
                           'fig'       'real'   []      0;
                           'uilist'    'cell'   []      {};
                           'title'     'string' []      '';
                           'userdata'  ''       []      [];
                           'geomvert'  'real'   []      [];
                           'inseth'    'real'   []      0.02; % x border absolute (5% of width)
                           'insetv'    'real'   []      0.02 }, 'supergui');
if isstr(g), error(g); end;
g.insetv = g.insetv/length(g.geomhoriz);

% create new figure
% -----------------
if g.fig == 0
	g.fig = figure('visible','off');
end;

% converting the geometry formats
% -------------------------------
if ~iscell( g.geomhoriz )
	oldgeom = g.geomhoriz;
	g.geomhoriz = {};
	for row = 1:length(oldgeom)
		g.geomhoriz = { g.geomhoriz{:} ones(1, oldgeom(row)) };
	end;
end;
if isempty(g.geomvert)
	g.geomvert = ones(1, length(g.geomhoriz));
end;

% setting relative width in percent
% ---------------------------------
for row = 1:length(g.geomhoriz)
	tmprow = g.geomhoriz{row};
	sumrow = sum(g.geomhoriz{row});
	g.geomhoriz{row} = 1.05*g.geomhoriz{row}/sumrow;
	g.geomhoriz{row} = g.geomhoriz{row} - g.inseth*(length(tmprow)-1)/length(tmprow);
end;

% setting relative height in percent
% ---------------------------------
sumcol = sum(g.geomvert);
ind1   = find(g.geomvert == 1); ind1 = ind1(1);
g.geomvert  = (1.03+0.003*sumcol)*g.geomvert/sumcol;
g.geomvert  = g.geomvert - g.insetv*(length(g.geomvert)-1)/length(g.geomvert);

% get axis coordinates
% --------------------
set(g.fig, 'menubar', 'none', 'numbertitle', 'off');		
pos = get(gca,'position'); % plot relative to current axes
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)]; % allow to use normalized position [0 100] for x and y
axis('off');

% creating guis
% -------------
counter = 1; % count the elements
outwidth = 0;
outheight = 0;
%height = 1.05/(length(g.geomhoriz)+1)*(1-g.insetv);
%posy = 1 - height - 1/length(g.geomhoriz)*g.insetv;
factmultx = 0;
factmulty = 0; %zeros(length(g.geomhoriz));
posy = 0.98+(0.003*sumcol)/2+g.insetv;
for row = 1:length(g.geomhoriz)

	% init
    posx = -0.05;
	clear rowhandle;
	tmprow = g.geomhoriz{row};
    height = g.geomvert(row);
	posy = posy - height - g.insetv;
	
	for column = 1:length(tmprow)

		width  = tmprow(column);
		try
			currentelem = g.uilist{ counter };
		catch
			fprintf('Warning: not all boxes were filled\n');
			return;
		end;		
		if ~isempty(currentelem)
            if ~strcmp(currentelem{2}, 'popupmenu') & ~strcmp(currentelem{2}, 'pushbutton')
                rowhandle(column) = uicontrol(g.fig, 'unit', 'normalized', 'position', ...
						                      [posx posy width height].*s+q, currentelem{:});
            else % force height to be unitary
                rowhandle(column) = uicontrol(g.fig, 'unit', 'normalized', 'position', ...
						[posx posy+height-(height+g.geomvert(ind1))/2 width g.geomvert(ind1)].*s+q, currentelem{:});
            end;
            
			% this simply compute a factor so that all uicontrol will be visible
			% ------------------------------------------------------------------
			style = get( rowhandle(column), 'style');			
			set( rowhandle(column), 'units', 'pixels');			
			curpos = get(rowhandle(column), 'position');
			curext = get(rowhandle(column), 'extent');
			if ~strcmp(style, 'edit') & ~strcmp(style, 'pushbutton')
				factmultx = max(factmultx, curext(3)/curpos(3));
			end;
            if  ~strcmp(style, 'listbox')
                factmulty = max(factmulty, curext(4)/curpos(4));
            end;
			set( rowhandle(column), 'units', 'normalized');			

			% Uniformize button text aspect (first letter must be Capital)
            % -----------------------------
            if strcmp(style, 'pushbutton')
                tmptext = get(rowhandle(column), 'string');
                if length(tmptext) > 1
                    if upper(tmptext(1)) ~= tmptext(1) | lower(tmptext(2)) ~= tmptext(2)
                        tmptext = lower(tmptext);
                        try, tmptext(1) = upper(tmptext(1)); catch, end;
                    end;
                end;
                set(rowhandle(column), 'string', tmptext);
            end;
        else 
			rowhandle(column) = 0;
		end;
		
		handlers{ row } = rowhandle;
		allhandlers(counter) = rowhandle(column);
		
		posx   = posx + width + g.inseth;
		counter = counter+1;
	end;
	%posy      = posy - height - 1/length(g.geomhoriz)*g.insetv; %compensate for inset 
end;

% adjustments
% -----------
factmultx = factmultx+0.1; % because some text was still hidden
if factmultx < 0.3
	factmultx = 0.3;
end;

% for MAC (magnify figures that have edit fields)
% -------
warning off;
try, 
	if strcmp(computer, 'MAC')
		hh = findobj(allhandlers, 'style', 'edit');
		if ~isempty(hh)
			factmulty = factmulty;
		end;
	elseif ~isunix % windows
		hh = findobj(allhandlers, 'style', 'edit');
		if ~isempty(hh)
			factmulty = factmulty*1.08;
		end;
    end;
catch, end;
warning on;	

% scale and replace the figure in the screen
% -----------------------------------------
pos = get(g.fig, 'position');
if factmulty > 1
	pos(2) = max(0,pos(2)+pos(4)-pos(4)*factmulty);
end;
pos(1) = pos(1)+pos(3)*(1-factmultx)/2;
pos(3) = pos(3)*factmultx;
pos(4) = pos(4)*factmulty;
set(g.fig, 'position', pos);

% vertical alignment to bottom for text
% ---------------------------------------
for index = 1:length(allhandlers)
	if allhandlers(index) ~= 0
		if strcmp(get(allhandlers(index), 'style'), 'text')
			curpos = get(allhandlers(index), 'position');
			curext = get(allhandlers(index), 'extent');
			set(allhandlers(index), 'position',[curpos(1) curpos(2) curpos(3) curext(4)]);
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
set(hh, 'horizontalalignment', 'left');

hh = findobj(allhandlers, 'style', 'edit');
set(hh, 'BackgroundColor', [1 1 1]); %, 'horizontalalignment', 'right');

hh =findobj(allhandlers, 'parent', g.fig, 'style', 'pushbutton');
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'parent', g.fig, 'style', 'popupmenu');
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'parent', g.fig, 'style', 'checkbox');
if isunix
	set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
	set(hh, 'foregroundcolor', GUITEXTCOLOR);	
else 
	set(hh, 'backgroundcolor', GUIBACKCOLOR);
	set(hh, 'foregroundcolor', GUITEXTCOLOR);
end;
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
