% supergui() - a comprehensive gui automatic builder. This function help
%              to create GUI very fast without bothering about the 
%              positions of the elements. After creating a geometry, 
%              elements just place themselves into the predefined 
%              locations. It is especially useful for figure where you
%              intend to put text button and descriptions.
%
% Usage:
%   >> [handles, height, allhandles ] = ...
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
%               of window size.
%   'insetv'   - vertical space between elements. Default is 2% 
%               of window height.
%   'spacing'  - [horiz vert] spacing in normalized units. Default 
%   'spacingtype' - ['absolute'|'proportional'] absolute means that the 
%               spacing values are fixed. Proportional means that they
%               depend on the number of element in a line.
%   'minwidth' - [integer] minimal width in pixels. Default is none.
%   'screenpos' - [x y] position of the right top corner of the graphic
%               interface. 'center' may also be used to center the GUI on
%               the screen.
%   'adjustbuttonwidth' - ['on'|'off'] adjust button width in the GUI.
%               Default is 'off'.
%
% Hint:
%    use 'print -mfile filemane' to save a matlab file of the figure.
%
% Output:
%    handles  - all the handles of the elements (in the same order as the
%                uilist input).
%    height    - advised height for the figure (so the text look nice).
%    allhandles - all the handles in object format
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

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [handlers, outheight, allhandlers] = supergui( varargin);

% handlers cell format
% allhandlers linear format
handlers = {};
outheight = 0;

if nargin < 2
	help supergui;
	return;
end

% get version and
% set additional parameters
% -------------------------
v = version;
indDot  = find(v == '.');
versnum = str2num(v(1:indDot(2)-1));
if versnum >= 7.14
     addParamFont = { 'fontsize' 12 };
else addParamFont = { };
end

warning off MATLAB:hg:uicontrol:ParameterValuesMustBeValid

% decoding input and backward compatibility
% -----------------------------------------
if ischar(varargin{1})
    options = varargin;
else
    options = { 'fig'      varargin{1} 'geomhoriz' varargin{2} ...
                'geomvert' varargin{3} 'uilist'    varargin(4:end) }; 
end
g = finputcheck(options, { 'geomhoriz' 'cell'   []      {};
                           'fig'       ''       []      0;
                           'geom'      'cell'   []      {};
                           'uilist'    'cell'   []      {};
                           'title'     'string' []      '';
                           'userdata'  ''       []      [];
                           'adjustbuttonwidth' 'string' { 'on' 'off' } 'off';
                           'geomvert'  'real'   []      [];
                           'screenpos' { 'real' 'string' } [] [];
                           'horizontalalignment'  'string'   { 'left','right','center' } 'left';
                           'minwidth'  'real'   []      10;
                           'borders'   'real'   []      [0.05 0.04 0.07 0.06];
                           'spacing'   'real'   []      [0.02 0.01];
                           'inseth'    'real'   []      0.02; % x border absolute (5% of width)
                           'insetv'    'real'   []      0.02 }, 'supergui');
if ischar(g), error(g); end
if ~isempty(g.geomhoriz)
    maxcount = sum(cellfun('length', g.geomhoriz));
    if maxcount ~= length(g.uilist)
        warning('Wrong size for ''geomhoriz'' input');
    end
    if ~isempty(g.geomvert)
        if length(g.geomvert) ~= length(g.geomhoriz)
            warning('Wrong size for ''geomvert'' input');
        end
    end
    g.insetv = g.insetv/length(g.geomhoriz);
end
if ~isempty(g.geom)
    if length(g.geom) ~= length(g.uilist)
        warning('Wrong size for ''geom'' input');
    end
    maxcount = length(g.geom);
end

% create new figure
% -----------------
if g.fig == 0
	g.fig = figure('visible','off');
end

% converting the geometry formats
% -------------------------------
if ~isempty(g.geomhoriz) && ~iscell( g.geomhoriz )
	oldgeom = g.geomhoriz;
	g.geomhoriz = {};
	for row = 1:length(oldgeom)
		g.geomhoriz = { g.geomhoriz{:} ones(1, oldgeom(row)) };
	end
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
        end
        incy = incy+g.geomvert(row);
    end
    g.borders(1:2) = g.borders(1:2)/maxhoriz*5;
    g.borders(3:4) = g.borders(3:4)/sumvert*10;
    g.spacing(1)   = g.spacing(1)/maxhoriz*5;
    g.spacing(2)   = g.spacing(2)/sumvert*10;
end

% disp new geometry
% -----------------
if 0
    fprintf('{ ...\n');
    for index = 1:length(g.geom)
        fprintf('{ %g %g [%g %g] [%g %g] } ...\n', g.geom{index}{1}, g.geom{index}{2}, ...
            g.geom{index}{3}(1), g.geom{index}{3}(2), g.geom{index}{4}(1), g.geom{index}{3}(2));
    end
    fprintf('};\n');
end

% get axis coordinates
% --------------------
try 
    set(g.fig, 'menubar', 'none', 'numbertitle', 'off');
catch
end
    
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
transformTextUIList = [];
for counter = 1:maxcount

	% init
	clear rowhandle;
    gm = g.geom{counter};
    [posx posy width height] = getcoord(gm{1}, gm{2}, gm{3}, gm{4}, g.borders, g.spacing);
%     str1 = sprintf('%d, %d, [%1.2f,%1.2f], [%1.2f,%1.2f]', gm{1}, gm{2}, gm{3}(1), gm{3}(2), gm{4}(1), gm{4}(2));
%     str2 = sprintf('%.4f, ', [posx posy width height]);
%     fprintf('%s -> %s -> %s\n', str1, str2, g.uilist{ counter }{4});
    
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
            hf1 = 3.6/2-0.3;
            hf2 = 0.7/2-0.3;
            allhandlers{counter} = uicontrol(g.fig, 'unit', 'normalized', 'position', ...
                                          [posx-width/2 posy+hf1*height width/2 0.005].*s+q, 'style', 'pushbutton', 'string', '');
            allhandlers{counter} = uicontrol(g.fig, 'unit', 'normalized', 'position', ...
                                          [posx-width/2 posy+hf2*height width/2 0.005].*s+q, 'style', 'pushbutton', 'string', '');
            allhandlers{counter} = uicontrol(g.fig, 'unit', 'normalized', 'position', ...
                                          [posx posy+hf2*height 0.005 (hf1-hf2+0.1)*height].*s+q, 'style', 'pushbutton', 'string', '');
            allhandlers{counter} = uicontrol(g.fig, 'unit', 'normalized', 'position', ...
                                          [posx posy+(hf1+hf2)/2*height width/2 0.005].*s+q, 'style', 'pushbutton', 'string', '');
            allhandlers{counter} = nan;
        else
            if strcmpi(currentelem{1}, 'width')
                 curwidth = currentelem{2};
                 currentelem(1:2) = [];
            else curwidth = 0;
            end
            if strcmpi(currentelem{1}, 'align')
                 align = currentelem{2};
                 currentelem(1:2) = [];
            else align = 'right';
            end
            transformText = false;
            if strcmpi(currentelem{2}, 'text2')
                transformTextUIList = [ transformTextUIList counter ];
                currentelem{2} = 'text';
            end
            if strcmpi(currentelem{1}, 'stickto')
                 stickto = currentelem{2};
                 currentelem(1:2) = [];
            else stickto = 'none';
            end
            if strcmpi(currentelem{1}, 'vertshift'), currentelem(1) = []; addvert = -height/2; 
            else                                                          addvert = 0;   
            end
            if strcmpi(currentelem{1}, 'vertexpand'), heightfactor = currentelem{2}; addvert = -(heightfactor-1)*height; currentelem(1:2) = []; 
            else                                      heightfactor = 1;   
            end
            
            % position adjustment depending on GUI type
            if ischar(currentelem{2}) && strcmpi(currentelem{2}, 'popupmenu')
                posy = posy-height/5;
            end
            if ischar(currentelem{2}) && strcmpi(currentelem{2}, 'text')
                posy = posy+height/5;
            end
                
            if strcmpi(currentelem{1}, 'function'),
                % property grid argument
                panel = uipanel('Title','','FontSize',12,'BackgroundColor','white','Position',[posx posy+addvert width height*heightfactor].*s+q);
                allhandlers{counter} = arg_guipanel(panel, currentelem{:});
            elseif strcmpi(currentelem{1}, 'panel'),
                % property grid argument
                uipanel(currentelem{2:end},'FontSize',12,'Position',[posx posy+addvert width height*heightfactor].*s+q);
                allhandlers{counter} = nan;
            elseif strcmpi(currentelem{1}, 'uitable'),
                uitable(g.fig, currentelem{2:end}, 'unit', 'normalized', 'Position',[posx posy+addvert width height*heightfactor].*s+q);
                allhandlers{counter} = nan;
            else
                if ~isempty(strmatch('fontsize', currentelem(1:2:end)))
                    allhandlers{counter} = uicontrol(g.fig, 'unit', 'normalized', 'position', ...
                        [posx posy+addvert width height*heightfactor].*s+q, currentelem{:});
                else
                    allhandlers{counter} = uicontrol(g.fig, 'unit', 'normalized', 'position', ...
                        [posx posy+addvert width height*heightfactor].*s+q, currentelem{:}, addParamFont{:});
                end
                
                % this simply compute a factor so that all uicontrol will be visible
                % ------------------------------------------------------------------
                style = get( allhandlers{counter}, 'style');
                set( allhandlers{counter}, 'units', 'pixels');
                curpos = get(allhandlers{counter}, 'position');
                curext = get(allhandlers{counter}, 'extent');
                if curwidth ~= 0
                    curwidth = curwidth/((factmultx-1)/1.85+1);
                    if strcmpi(align, 'right')
                        curpos(1) = curpos(1)+curpos(3)-curwidth;
                    elseif strcmpi(align, 'center')
                        curpos(1) = curpos(1)+curpos(3)/2-curwidth/2;
                    end
                    set(allhandlers{counter}, 'position', [ curpos(1) curpos(2) curwidth curpos(4) ]);
                    if strcmpi(stickto, 'on')
                        set( allhandlers{counter-1}, 'units', 'pixels');
                        curpos2 = get(allhandlers{counter-1}, 'position');
                        set(allhandlers{counter-1}, 'position', [ curpos(1)-curpos2(3)-10 curpos2(2) curpos2(3) curpos2(4) ]);
                        set( allhandlers{counter-1}, 'units', 'normalized');
                    end
                    curext(3) = curwidth;
                end
                set( allhandlers{counter}, 'units', 'normalized');			
            end

            if ~strcmp(style, 'edit') && (~strcmp(style, 'pushbutton') || strcmpi(g.adjustbuttonwidth, 'on'))
                if strcmpi( currentelem{3}, 'string') && (~iscell(currentelem{4}) || ...
                        (~isempty(currentelem{4}) && isempty(findstr('html', currentelem{4}{1}))))
                    %tmp = curext(3)/curpos(3);
                    %if tmp > 3*factmultx && factmultx > 0, adsfasd; end
                    if exist('OCTAVE_VERSION','builtin') ~= 0
                      curtxt = get(allhandlers{counter}, 'string');
                      if ischar(curtxt)
                        nLines = length(find(curtxt == 10))+1;
                        curext(3) = curext(3)/nLines;
                      end
                    end
                    
                    if curext(3)/curpos(3) < 3
                        factmultx = max(factmultx, curext(3)/curpos(3));
                    end
                    tmptxt = get(allhandlers{counter}, 'string');
%                   % text to get which UI modify horizontal size
%                    if iscell(tmptxt) fprintf('%s -> %f\n', 'cell', curext(3)/curpos(3));
%                    else              fprintf('%s -> %f\n', tmptxt, curext(3)/curpos(3));
%                    end                      
                    if strcmp(style, 'pushbutton'), factmultx = factmultx*1.1; end
                end
            end
            if  ~strcmp(style, 'listbox') &&  ~strcmp(style, 'popupmenu')
                factmulty = max(factmulty, curext(4)/curpos(4));
            end
            
            % Uniformize button text aspect (first letter must be uppercase)
            % -----------------------------
            if strcmp(style, 'pushbutton') && ishandle(allhandlers{counter})
                tmptext = get(allhandlers{counter}, 'string');
                if length(tmptext) > 1
                    if upper(tmptext(1)) ~= tmptext(1) || lower(tmptext(2)) ~= tmptext(2) && ~strcmpi(tmptext, 'STATS')
                        tmptext = lower(tmptext);
                        try, tmptext(1) = upper(tmptext(1)); catch, end
                    end
                end
                set(allhandlers{counter}, 'string', tmptext);
            end
        end
    else 
        allhandlers{counter} = nan;
    end
end

% adjustments
% -----------
factmultx = factmultx*1.1;% because some text was still hidden
if factmultx < 0.1
	factmultx = 0.1;
end

% for MAC (magnify figures that have edit fields)
% -------
warning off;
try, 
    comp = computer;
	if length(comp) > 2 && strcmpi(comp(1:3), 'MAC')
        factmulty = factmulty*1.5;
	elseif ~isunix % windows
        factmulty = factmulty*1.08;
    end
catch, end
if exist('OCTAVE_VERSION', 'builtin') == 0
    factmulty = factmulty*0.9; % global shinking
else
    factmulty = factmulty*1.5;
end
warning on;	

% scale and replace the figure in the screen
% -----------------------------------------
pos = get(g.fig, 'position');
if factmulty > 1
	pos(2) = max(0,pos(2)+pos(4)-pos(4)*factmulty);
end
pos(1) = max(0,pos(1)+pos(3)*(1-factmultx)/2);
pos(3) = max(pos(3)*factmultx, g.minwidth);
pos(4) = pos(4)*factmulty;
set(g.fig, 'position', pos);

% vertical alignment to bottom for text (isnumeric by ishanlde was changed here)
% ---------------------------------------
for index = 1:length(allhandlers)
	if ishandle(allhandlers{index})
		if strcmp(get(allhandlers{index}, 'style'), 'text')
            set(allhandlers{index}, 'unit', 'pixel');
			curpos = get(allhandlers{index}, 'position');
			curext = get(allhandlers{index}, 'extent');
			set(allhandlers{index}, 'position', [curpos(1) curpos(2)-4 curpos(3) curext(4)]);
            set(allhandlers{index}, 'unit', 'normalized');
		end
	end
end

% setting defaults colors
%------------------------
try, icadefs;
catch,
	GUIBACKCOLOR        =  [.8 .8 .8];     
	GUIPOPBUTTONCOLOR   = [.8 .8 .8];    
	GUITEXTCOLOR        = [0 0 0];
end

numobjects = cellfun(@ishandle, allhandlers); % (isnumeric by ishanlde was changed here)
allhandlersnum = [ allhandlers{numobjects} ];
hh = findobj(allhandlersnum,'flat', 'parent', g.fig, 'style', 'text');
%set(hh, 'BackgroundColor', get(g.fig, 'color'), 'horizontalalignment', 'left');
set(hh, 'Backgroundcolor', GUIBACKCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
try
set(g.fig, 'color',GUIBACKCOLOR );
catch
end
set(hh, 'horizontalalignment', g.horizontalalignment);

hh = findobj(allhandlersnum,'flat', 'style', 'edit');
set(hh, 'BackgroundColor', [1 1 1]); %, 'horizontalalignment', 'right');

hh =findobj(allhandlersnum,'flat', 'parent', g.fig, 'style', 'pushbutton');
comp = computer;
if length(comp) < 3 || ~strcmpi(comp(1:3), 'MAC') % this puts the wrong background on macs
    set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
    set(hh, 'foregroundcolor', GUITEXTCOLOR);
end
hh =findobj(allhandlersnum,'flat', 'parent', g.fig, 'style', 'popupmenu');
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlersnum,'flat', 'parent', g.fig, 'style', 'checkbox');
set(hh, 'backgroundcolor', GUIBACKCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlersnum,'flat', 'parent', g.fig, 'style', 'listbox');
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlersnum,'flat', 'parent', g.fig, 'style', 'radio');
set(hh, 'foregroundcolor', GUITEXTCOLOR);
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(g.fig, 'visible', 'on');

% screen position
% ---------------
if ~isempty(g.screenpos)
    pos = get(g.fig, 'position');
    if isnumeric(g.screenpos)
        set(g.fig, 'position', [ g.screenpos pos(3) pos(4)]);
    else
        screenSize = get(0, 'screensize');
        pos(1) = (screenSize(3)-pos(3))/2;
        pos(2) = (screenSize(4)-pos(4))/2+pos(4);
        set(g.fig, 'position', pos);
    end
end;    

% set userdata and title
% ----------------------
if ~isempty(g.userdata), set(g.fig, 'userdata', g.userdata); end
if ~isempty(g.title   ), set(g.fig, 'name',     g.title   ); end

% transform selected text UI
% --------------------------
a = axes('position', [0 0 1 1]);
for iList = transformTextUIList
    set(allhandlers{iList}, 'unit', 'normalized');
    sTmp = get(allhandlers{iList},'string');
    pTmp = get(allhandlers{iList},'position');
    eTmp = get(allhandlers{iList},'extent');
    tTmp = get(allhandlers{iList},'tag');
    % Remove the UICONTROL
    delete(allhandlers{iList});
    % Replace it with a TEXT object
    allhandlers{iList} = text(pTmp(1),pTmp(2),sTmp, 'interpreter','latex', 'tag', tTmp);
end
set(a, 'visible', 'off');

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
