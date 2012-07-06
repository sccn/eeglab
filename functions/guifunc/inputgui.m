% inputgui() - A comprehensive gui automatic builder. This function helps
%              to create GUI very quickly without bothering about the 
%              positions of the elements. After creating a geometry, 
%              elements just place themselves in the predefined 
%              locations. It is especially useful for figures in which 
%              you intend to put text buttons and descriptions.
%
% Usage:
%   >> [ outparam ] = inputgui( 'key1', 'val1', 'key2', 'val2', ... );
%   >> [ outparam userdat strhalt outstruct] = ...
%             inputgui( 'key1', 'val1', 'key2', 'val2', ... );
% 
% Inputs:
%   'geom'       - cell array of cell array of integer vector. Each cell
%                  array defines the coordinate of a given input in the 
%                  following manner: { nb_row nb_col [x_topcorner y_topcorner]
%                  [x_bottomcorner y_bottomcorner] };
%   'geometry'   - cell array describing horizontal geometry. This corresponds 
%                  to the supergui function input 'geomhoriz'
%   'geomvert'   - vertical geometry argument, this argument is passed on to
%                  the supergui function
%   'uilist'     - list of uicontrol lists describing elements properties
%                  { { ui1 }, { ui2 }... }, { 'uiX' } being GUI matlab 
%                  uicontrol arguments such as { 'style', 'radiobutton', 
%                  'String', 'hello' }. See Matlab function uicontrol() for details.
%   'helpcom'    - optional help command 
%   'title'      - optional figure title
%   'userdata'   - optional userdata input for the figure
%   'mode'       - ['normal'|'noclose'|'plot' fignumber]. Either wait for
%                  user to press OK or CANCEL ('normal'), return without
%                  closing window input ('noclose'), only draw the gui ('plot')
%                  or process an existing window which number is given as 
%                  input (fignumber). Default is 'normal'.
%   'eval'       - [string] command to evaluate at the end of the creation 
%                  of the GUI but before waiting for user input. 
%   'screenpos'  - see supergui.m help message.
%   'skipline'   - ['on'|'off'] skip a row before the "OK" and "Cancel"
%                  button. Default is 'on'.
%
% Output:
%   outparam   - list of outputs. The function scans all lines and
%                add up an output for each interactive uicontrol, i.e
%                edit box, radio button, checkbox and listbox.
%   userdat    - 'userdata' value of the figure.
%   strhalt    - the function returns when the 'userdata' field of the
%                button with the tag 'ok' is modified. This returns the
%                new value of this field.
%   outstruct  - returns outputs as a structure (only tagged ui controls
%                are considered). The field name of the structure is
%                the tag of the ui and contain the ui value or string.
%
% Note: the function also adds three buttons at the bottom of each 
%       interactive windows: 'CANCEL', 'HELP' (if callback command
%       is provided) and 'OK'.
%
% Example:
%   res = inputgui('geometry', { 1 1 }, 'uilist', ...
%                         { { 'style' 'text' 'string' 'Enter a value' } ...
%                           { 'style' 'edit' 'string' '' } });
%
%   res = inputgui('geom', { {2 1 [0 0] [1 1]} {2 1 [1 0] [1 1]} }, 'uilist', ...
%                         { { 'style' 'text' 'string' 'Enter a value' } ...
%                           { 'style' 'edit' 'string' '' } });
%
% Author: Arnaud Delorme, CNL / Salk Institute, La Jolla, 1 Feb 2002
%
% See also: supergui(), eeglab()

% Copyright (C) Arnaud Delorme, CNL/Salk Institute, 27 Jan 2002, arno@salk.edu
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

function [result, userdat, strhalt, resstruct] = inputgui( varargin);

if nargin < 2
   help inputgui;
   return;
end;	

% decoding input and backward compatibility
% -----------------------------------------
if isstr(varargin{1})
    options = varargin;
else
    options = { 'geometry' 'uilist' 'helpcom' 'title' 'userdata' 'mode' 'geomvert' };
    options = { options{1:length(varargin)}; varargin{:} };
    options = options(:)';
end;

% checking inputs
% ---------------
g = finputcheck(options, { 'geom'     'cell'                []      {}; ...
                           'geometry' {'cell','integer'}    []      []; ...
                           'uilist'   'cell'                []      {}; ...
                           'helpcom'  { 'string','cell' }   { [] [] }      ''; ...
                           'title'    'string'              []      ''; ...
                           'eval'     'string'              []      ''; ...
                           'skipline' 'string'              { 'on' 'off' } 'on'; ...
                           'userdata' ''                    []      []; ...
                           'getresult' 'real'               []      []; ...
                           'screenpos' ''                   []      []; ...
                           'mode'     ''                    []      'normal'; ...
                           'geomvert' 'real'                []       [] ...
                          }, 'inputgui');
if isstr(g), error(g); end;

if isempty(g.getresult)
    if isstr(g.mode)
        fig = figure('visible', 'off');
        set(fig, 'name', g.title);
        set(fig, 'userdata', g.userdata);
        if ~iscell( g.geometry )
            oldgeom = g.geometry;
            g.geometry = {};
            for row = 1:length(oldgeom)
                g.geometry = { g.geometry{:} ones(1, oldgeom(row)) };
            end;
        end
        if strcmpi(g.skipline, 'on')
             g.geometry = { g.geometry{:} [1] [1 1 1] }; % add button to geometry
        else g.geometry = { g.geometry{:} [1 1 1] }; % add button to geometry
        end;
        if ~isempty(g.geom)
            for ind = 1:length(g.geom)
                g.geom{ind}{2} = g.geom{ind}{2}+2;
            end;
            g.geom = { g.geom{:}, ...
                      {1 g.geom{1}{2} [0 g.geom{1}{2}-2] [1 1] }, ... 
                      {3 g.geom{1}{2} [0 g.geom{1}{2}-1] [1 1] }, ... 
                      {3 g.geom{1}{2} [1 g.geom{1}{2}-1] [1 1] }, ...
                      {3 g.geom{1}{2} [2 g.geom{1}{2}-1] [1 1] } };
        end;

        % add the three buttons (CANCEL HELP OK) at the bottom of the GUI
        % ---------------------------------------------------------------
        if strcmpi(g.skipline, 'on'),  g.uilist = { g.uilist{:}, {} }; end;
        options = { 'width' 80 'stickto' 'on' };
        if ~isempty(g.helpcom)
            if ~iscell(g.helpcom) | isempty(g.geom)
                g.uilist = { g.uilist{:}, { 'width' 80 'align' 'left' 'Style', 'pushbutton', 'string', 'Help', 'tag', 'help', 'callback', g.helpcom } };
            else
                g.uilist = { g.uilist{:}, { 'width' 80 'align' 'left' 'Style', 'pushbutton', 'string', 'Help gui', 'callback', g.helpcom{1} } };
                g.uilist = { g.uilist{:}, { 'width' 80 'align' 'left' 'Style', 'pushbutton', 'string', 'More help', 'callback', g.helpcom{2} } };
                g.geometry{end} = [1 1 1 1];
            end;
        else
            g.uilist = { g.uilist{:}, {} };
        end;
        g.uilist = { g.uilist{:}, { 'width' 80 'align' 'right' 'Style', 'pushbutton', 'string', 'Cancel', 'tag' 'cancel' 'callback', 'close gcbf' } };
        g.uilist = { g.uilist{:}, { 'width' 80 'align' 'right' 'stickto' 'on' 'Style', 'pushbutton', 'tag', 'ok', 'string', 'OK', 'callback', 'set(gcbo, ''userdata'', ''retuninginputui'');' } };
        if ~isempty(g.geom)
            [tmp tmp2 allobj] = supergui( 'fig', fig, 'minwidth', 200, 'geom', g.geom, 'uilist', g.uilist, 'screenpos', g.screenpos );
        elseif isempty(g.geomvert)
            [tmp tmp2 allobj] = supergui( 'fig', fig, 'minwidth', 200, 'geomhoriz', g.geometry, 'uilist', g.uilist, 'screenpos', g.screenpos );
        else
            if strcmpi(g.skipline, 'on'),  g.geomvert = [g.geomvert(:)' 1]; end;
            [tmp tmp2 allobj] = supergui( 'fig', fig, 'minwidth', 200, 'geomhoriz', g.geometry, 'uilist', g.uilist, 'screenpos', g.screenpos, 'geomvert', [g.geomvert(:)' 1] );
        end;
    else 
        fig = g.mode;
        set(findobj('parent', fig, 'tag', 'ok'), 'userdata', []);
        allobj = findobj('parent',fig);
        allobj = allobj(end:-1:1);
    end;

    % evaluate command before waiting?
    % --------------------------------
    if ~isempty(g.eval), eval(g.eval); end;

    % create figure and wait for return
    % ---------------------------------
    if isstr(g.mode) & (strcmpi(g.mode, 'plot') | strcmpi(g.mode, 'return') )
        if strcmpi(g.mode, 'plot')
           return; % only plot and returns
        end;
    else 
        waitfor( findobj('parent', fig, 'tag', 'ok'), 'userdata');
    end;
else
    fig = g.getresult;
    allobj = findobj('parent',fig);
    allobj = allobj(end:-1:1);
end;

result = {};
userdat = [];
strhalt = '';
resstruct = [];
try, findobj(fig); % figure still exist ?
catch, return; end;
strhalt = get(findobj('parent', fig, 'tag', 'ok'), 'userdata');

% output parameters
% -----------------
counter = 1;
for index=1:length(allobj)
   try,
      objstyle = get(allobj( index ), 'style');
      switch lower( objstyle )
      case { 'listbox', 'checkbox', 'radiobutton' 'popupmenu' 'radio' }
         result{counter} = get( allobj( index ), 'value');
         counter = counter+1;
      case 'edit' 
         result{counter} = get( allobj( index ), 'string');
         counter = counter+1;
      end;
   catch, end;
end;   
userdat = get(fig, 'userdata');
if nargout >= 4
	resstruct = myguihandles(fig);
end;

if isempty(g.getresult) && isstr(g.mode) && ( strcmp(g.mode, 'normal') || strcmp(g.mode, 'return') )
	close(fig);
end;
drawnow; % for windows

% function for gui res
% --------------------
function g = myguihandles(fig)
	g = [];
	h = findobj('parent', fig);
	for index = 1:length(h)
		if ~isempty(get(h(index), 'tag'))
			try, 
				switch get(h(index), 'style')
				 case 'edit', g = setfield(g, get(h(index), 'tag'), get(h(index), 'string'));
				 case { 'value' 'radio' 'checkbox' 'listbox' 'popupmenu' 'radiobutton'  }, ...
					  g = setfield(g, get(h(index), 'tag'), get(h(index), 'value'));
				end;
			catch, end;
		end;
	end;
