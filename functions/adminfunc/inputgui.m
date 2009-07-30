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
%   'geometry'   - this corresponds to supergui() key-val input 'geomhoriz'
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
%   'geomvert'   - vertical geometry argument, this argument is passed on to
%                  supergui()
%   'eval'       - [string] command to evaluate at the end of the creation 
%                  of the GUI but before waiting for user input. 
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
% Author: Arnaud Delorme, CNL / Salk Institute, La Jolla, 1 Feb 2002
%
% See also: supergui(), eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.38  2009/02/09 11:03:20  arno
% Fix CVS problem for Windows
%
% Revision 1.37  2007/03/05 18:55:00  arno
% fixed helpcom for string and cell array
%
% Revision 1.36  2007/01/08 05:21:34  toby
% help edit
%
% Revision 1.35  2007/01/06 08:45:37  toby
% Geometry input altered for consistency with supergui.m help
%
% Revision 1.34  2007/01/06 06:11:11  toby
% Help text update
%
% Revision 1.33  2006/11/15 21:16:30  arno
% typo in header
%
% Revision 1.32  2006/02/10 23:36:56  arno
% assign popupmenu too
%
% Revision 1.31  2006/01/13 00:18:33  arno
% nothing
%
% Revision 1.30  2006/01/05 21:15:03  arno
% adding popupmenu
%
% Revision 1.29  2005/11/09 23:07:01  arno
% nothing
%
% Revision 1.28  2005/11/09 22:44:56  arno
% fixing calling format and header
%
% Revision 1.27  2005/11/09 22:30:32  arno
% new calling format
%
% Revision 1.26  2005/03/03 17:33:57  arno
% *** empty log message ***
%
% Revision 1.25  2005/03/02 19:57:11  hilit
% added a listbox object to the returned output structure
%
% Revision 1.24  2004/11/10 17:04:04  arno
% nothing
%
% Revision 1.23  2003/03/12 02:44:36  arno
% tow buttons
%
% Revision 1.22  2003/03/12 02:41:19  arno
% adding help gui button
%
% Revision 1.21  2003/02/21 16:49:22  arno
% nothing
%
% Revision 1.20  2002/12/24 01:26:19  arno
% adding 'plot' option
%
% Revision 1.19  2002/11/13 17:06:44  scott
% hj
% help msg
%
% Revision 1.18  2002/11/13 00:53:43  arno
% replace gcf by fig
%
% Revision 1.17  2002/11/13 00:49:31  arno
% tag the figure
% /
%
% Revision 1.16  2002/10/15 17:23:21  arno
% debug drawnow
%
% Revision 1.15  2002/10/15 16:53:37  arno
% add drawnow for windows
%
% Revision 1.14  2002/08/23 17:41:13  arno
% implementing new option
%
% Revision 1.13  2002/08/14 18:17:11  arno
% new supergui arg
%
% Revision 1.12  2002/08/13 21:44:55  arno
% debug
%
% Revision 1.11  2002/08/13 18:21:24  arno
% passes on geomvert
%
% Revision 1.10  2002/08/13 17:29:07  arno
% new supergui call
%
% Revision 1.9  2002/08/12 18:24:39  arno
% empty help
%
% Revision 1.8  2002/07/18 17:14:48  arno
% no modif
%
% Revision 1.7  2002/04/27 01:20:08  arno
% debugging mode
%
% Revision 1.6  2002/04/27 00:49:02  arno
% adding extra output parameter
%
% Revision 1.5  2002/04/27 00:17:25  arno
% debugging function call
%
% Revision 1.4  2002/04/26 23:32:04  arno
% updated mode argument value and processing
%
% Revision 1.3  2002/04/26 23:24:37  arno
% adding mode
%
% Revision 1.2  2002/04/09 03:44:35  arno
% adding input userdata
%
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%
% 02/15/02 add userdat option -ad
% 02/16/02 add figure title option -ad

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
g = finputcheck(options, { 'geometry' {'cell','integer'}    []      []; ...
                           'uilist'   'cell'                []      {}; ...
                           'helpcom'  { 'string' 'cell' }   { [] [] }      ''; ...
                           'title'    'string'              []      ''; ...
                           'eval'     'string'              []      ''; ...
                           'userdata' ''                    []      []; ...
                           'mode'     ''                    []      'normal'; ...
                           'geomvert' 'real'                []       [] ...
                          }, 'inputgui');
if isstr(g), error(g); end;

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
	g.geometry = { g.geometry{:} [1] [1 1 1] }; % add button to geometry
	
	% add the three buttons (CANCEL HELP OK) at the bottom of the GUI
	% ---------------------------------------------------------------
	g.uilist = { g.uilist{:}, {}, { 'Style', 'pushbutton', 'string', 'Cancel', 'callback', 'close gcbf' } };
	if ~isempty(g.helpcom)
        if ~iscell(g.helpcom)
            g.uilist = { g.uilist{:}, { 'Style', 'pushbutton', 'string', 'Help', 'callback', g.helpcom } };
        else
            g.uilist = { g.uilist{:}, { 'Style', 'pushbutton', 'string', 'Help gui', 'callback', g.helpcom{1} } };
            g.uilist = { g.uilist{:}, { 'Style', 'pushbutton', 'string', 'More help', 'callback', g.helpcom{2} } };
            g.geometry{end} = [1 1 1 1];
        end;
	else
		g.uilist = { g.uilist{:}, {} };
	end;   
	g.uilist = { g.uilist{:}, { 'Style', 'pushbutton', 'tag', 'ok', 'string', 'OK', 'callback', 'set(gcbo, ''userdata'', ''retuninginputui'');' } };
	if isempty(g.geomvert)
		[tmp tmp2 allobj] = supergui( fig, g.geometry, [], g.uilist{:} );
	else
		[tmp tmp2 allobj] = supergui( fig, g.geometry, [g.geomvert(:)' 1 1], g.uilist{:} );
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
      case { 'listbox', 'checkbox', 'radiobutton' 'popupmenu' }
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

if isstr(g.mode) & ( strcmp(g.mode, 'normal') | strcmp(g.mode, 'return') )
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
				 case { 'value' 'radio' 'checkbox' 'listbox' 'popupmenu' }, ...
					  g = setfield(g, get(h(index), 'tag'), get(h(index), 'value'));
				end;
			catch, end;
		end;
	end;
