% supergui() - a comprehensive gui automatic builder. This function help
%              to create GUI very fast without bothering about the 
%              positions of the elements. After creating a geometry, 
%              elements just place themselves into the predefined 
%              locations. It is especially usefull for figure where you
%              intend to put text button and descriptions.
%
% Usage:
%   >> [handlers, width, height ] = ...
%             supergui( geomx, geomy, { arguments1 }, { arguments2 }... );
% 
% Inputs:
%   geomx   - cell array describing the geometry of the elements
%             in the figure. For instance, [2 3 2] means that the
%             figures will have 3 rows, with 2 elements in the first
%             and last row and 3 elements in the second row.
%             An other syntax is { [2 8] [1 2 3] } which means
%             that figures will have 2 rows, the first one with 2
%             elements of relative width 2 and 8 (20% and 80%). The
%             second row will have 3 elements of relative size 1, 2 
%             and 3.
%   geomy  - [array] describting geometry for the rows. For instance
%            [1 2 1] means that the second row will be twice the height
%            of the other ones. If [], all the lines have the same height.
%  {argument} - GUI matlab element arguments. Ex { 'style', 
%               'radiobutton', 'String', 'hello' }.
%
% Hint:
%    use 'print -mfile filemane' to save a matlab file of the figure.
%
% Output:
%    handlers  - all the handler of the elements (in the same form as the
%                geometry cell input.
%    height    - adviced widht for the figure (so the text look nice).   
%    height    - adviced height for the figure (so the text look nice).   
%
% Example:
%    figure;   
%    supergui( [1 1], [], { 'style', 'radiobutton', 'string', 'radio' }, ...
%        { 'style', 'pushbutton', 'string', 'push' });
%      
% Author: Arnaud Delorme, CNL / Salk Institute, La Jolla, 2001
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

function [handlers, outheight, allhandlers] = supergui( geomx, geomy, varargin);

% handlers cell format
% allhandlers linear format

INSETX = 0.05; % x border absolute (5% of width)
INSETY = 0.01;  % y border relative (50% of heigth)  

if nargin < 2
	help supergui;
	return;
end;

% converting the geometry formats
% -------------------------------
if ~iscell( geomx )
	oldgeom = geomx;
	geomx = {};
	for row = 1:length(oldgeom)
		geomx = { geomx{:} ones(1, oldgeom(row)) };
	end;
end;
if isempty(geomy)
	geomy = ones(1, length(geomx));
end;

% setting relative width in percent
% ---------------------------------
for row = 1:length(geomx)
	tmprow = geomx{row};
	sumrow = sum(geomx{row});
	geomx{row} = 1.05*geomx{row}/sumrow;
	geomx{row} = geomx{row} - INSETX*(length(tmprow)-1)/length(tmprow);
end;

% setting relative height in percent
% ---------------------------------
sumcol = sum(geomy);
geomy  = (1.03+0.003*sumcol)*geomy/sumcol;
geomy  = geomy - INSETY*(length(geomy)-1)/length(geomy);

% $$$ % counting rows
% $$$ % -------------
% $$$ nbrows = 0;
% $$$ counter = 1; % count the elements
% $$$ for row = 1:length(geomx)
% $$$ 	nbrowtmp = 1;
% $$$ 	for column = 1:length(tmprow)
% $$$ 		currentelem = varargin{ counter };
% $$$ 		if ~isempty(currentelem)
% $$$ 			try, 
% $$$ 				if size(currentelem,1) > 1
% $$$ 			nbrowtmp = 2;
% $$$ 		end;
% $$$ 		counter = counter+1;
% $$$ 	end;
% $$$ 	nbrows = nbrow+nbrowtmp;
% $$$ end;

% get axis coordinates
% --------------------
set(gcf, 'menubar', 'none', 'numbertitle', 'off');		
pos = get(gca,'position'); % plot relative to current axes
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)]; % allow to use normalized position [0 100] for x and y
axis('off');

% creating guis
% -------------
counter = 1; % count the elements
outwidth = 0;
outheight = 0;
%height = 1.05/(length(geomx)+1)*(1-INSETY);
%posy = 1 - height - 1/length(geomx)*INSETY;
factmultx = 0;
factmulty = 0; %zeros(length(geomx));
posy = 0.98+(0.003*sumcol)/2+INSETY;
for row = 1:length(geomx)

	% init
    posx = -0.05;
	clear rowhandle;
	tmprow = geomx{row};
    height = geomy(row);
	posy = posy - height - INSETY;
	
	for column = 1:length(tmprow)

		width  = tmprow(column);
		try
			currentelem = varargin{ counter };
		catch
			fprintf('Warning: not all boxes were filled\n');
			return;
		end;		
		if ~isempty(currentelem)
			rowhandle(column) = uicontrol( 'unit', 'normalized', 'position', ...
						                      [posx posy width height].*s+q, currentelem{:});
						
			% this simply compute a factor so that all uicontrol will be visible
			% ------------------------------------------------------------------
			style = get( rowhandle(column), 'style');			
			set( rowhandle(column), 'units', 'pixels');			
			curpos = get(rowhandle(column), 'position');
			curext = get(rowhandle(column), 'extent');
			if ~strcmp(style, 'edit') & ~strcmp(style, 'pushbutton')
				factmultx = max(factmultx, curext(3)/curpos(3));
			end;
			if ~strcmp(style, 'pushbutton')
				factmulty = max(factmulty, curext(4)/curpos(4));
			end;
			set( rowhandle(column), 'units', 'normalized');			
        else 
			rowhandle(column) = 0;
		end;
		
		handlers{ row } = rowhandle;
		allhandlers(counter) = rowhandle(column);
		
		posx   = posx + width + INSETX;
		counter = counter+1;
	end;
	%posy      = posy - height - 1/length(geomx)*INSETY; %compensate for inset 
end;

% recompute ordinates based on extet info
% ---------------------------------------

%scale and replace the figure in the screen
pos = get(gcf, 'position');
if factmulty > 1
	pos(2) = max(0,pos(2)+pos(4)-pos(4)*factmulty);
end;
pos(1) = pos(1)+pos(3)*(1-factmultx)/2;
pos(3) = pos(3)*factmultx;
pos(4) = pos(4)*factmulty;
set(gcf, 'position', pos);


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

hh = findobj(allhandlers, 'parent', gcf, 'style', 'text');
%set(hh, 'BackgroundColor', get(gcf, 'color'), 'horizontalalignment', 'left');
set(hh, 'Backgroundcolor', GUIBACKCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
set(gcf, 'color', get(hh(1), 'BackgroundColor'));
set(hh, 'horizontalalignment', 'left');

hh = findobj(allhandlers, 'style', 'edit');
set(hh, 'BackgroundColor', [1 1 1]); %, 'horizontalalignment', 'right');

hh =findobj(allhandlers, 'parent', gcf, 'style', 'pushbutton');
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'parent', gcf, 'style', 'checkbox');
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'parent', gcf, 'style', 'listbox');
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'parent', gcf, 'style', 'radio');
set(hh, 'foregroundcolor', GUITEXTCOLOR);
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);

return;
