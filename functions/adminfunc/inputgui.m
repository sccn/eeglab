% inputgui() - a comprehensive gui automatic builder. This function help
%              to create GUI very fast without bothering about the 
%              positions of the elements. After creating a geometry, 
%              elements just place themselves into the predefined 
%              locations. It is especially usefull for figure where you
%              intend to put text button and descriptions.
%
% Usage:
%   >> [ outparam ] = inputgui( geometry, listui );
%   >> [ outparam userdat strhalt] = ...
%             inputgui( geometry, listui, help, title, userdat, mode );
% 
% Inputs:
%   geometry   - see supergui()
%   listui     - list of uicontrol lists describing elements properties
%                { { ui1 }, { ui2 }... }, { 'uiX' } being GUI matlab 
%                uicontrol arguments such as { 'style', 'radiobutton', 
%                'String', 'hello' }. See supergui() for details.
%   help       - optional help command 
%   title      - optional figure title
%   userdat    - optional userdata input for the figure
%   mode       - ['normal'|'noclose'|fignumber], either wait for user to press
%                OK or CANCEL ('normal'), return without closing window
%                input ('noclose'), or process an existing 
%                window which number is given as input (fignumber). 
%                Default is 'normal'.
%
% Output:
%   outparam   - list of outputs. The function scan all lines and
%                add up an output for each interactive uicontrol, i.e
%                edit box, radio button, checkbox and listbox.
%   userdat    - 'userdata' value of the figure.
%   strhalt    - the funciton returns when the 'userdata' field of the
%                button with the tag 'ok' is modified. THis return the
%                new value of this field.
%
% Note: the function also add three buttons at the bottom of each 
%       interactive windows: 'CANCEL', 'HELP' (if callback command
%       is provided) and 'OK'.
%
% Author: Arnaud Delorme, CNL / Salk Institute, La Jolla, 1 Feb 2002
%
% See also: supergui(), eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 27 Jan 2002, arno@salk.edu
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

function [result, userdat, strhalt] = inputgui( geometry, listui, helpcom, mytitle, userdat, mode);

if nargin < 2
   help inputgui;
   return;
end;	
if exist('mode') ~= 1
	mode = 'normal';
end;

if isstr(mode)
	fig = figure;
	if exist('mytitle') == 1, set(gcf, 'name', mytitle); end;
	if exist('userdat') == 1, set(gcf, 'userdata', userdat); end; 
	geometry = { geometry{:} [1] [1 1 1] }; % add button to geometry
	
	% add the three buttons
	% ---------------------
	listui = { listui{:}, {}, { 'Style', 'pushbutton', 'string', 'Cancel', 'callback', 'close gcbf' } };
	if exist('helpcom') == 1 & ~isempty(helpcom)
		listui = { listui{:}, { 'Style', 'pushbutton', 'string', 'Help', 'callback', helpcom } };
	else
		listui = { listui{:}, {} };
	end;   
	listui = { listui{:}, { 'Style', 'pushbutton', 'tag', 'ok', 'string', 'OK', 'callback', 'set(gcbo, ''userdata'', ''retuninginputui'');' } };
	[tmp tmp2 allobj] = supergui( geometry, listui{:} );
else 
	fig = mode;
	set(findobj('parent', fig, 'tag', 'ok'), 'userdata', []);
	allobj = findobj('parent',fig);
	allobj = allobj(end:-1:1);
end;

% create figure and wait for return
% ---------------------------------
waitfor( findobj('parent', fig, 'tag', 'ok'), 'userdata');

result = {};
userdat = [];
strhalt = '';
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
      case { 'listbox', 'checkbox', 'radiobutton' }
         result{counter} = get( allobj( index ), 'value');
         counter = counter+1;
      case 'edit' 
         result{counter} = get( allobj( index ), 'string');
         counter = counter+1;
      end;
   catch, end;
end;   
userdat = get(fig, 'userdata');
if isstr(mode) & strcmp(mode, 'normal')
	close(fig);
end;