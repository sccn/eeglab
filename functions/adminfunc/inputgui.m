% inputgui() - a comprehensive gui automatic builder. This function help
%              to create GUI very fast without bothering about the 
%              positions of the elements. After creating a geometry, 
%              elements just place themselves into the predefined 
%              locations. It is especially usefull for figure where you
%              intend to put text button and descriptions.
%
% Usage:
%   >> [ outparam userdat ] = ...
%             inputgui( geometry, listui, help, title );
% 
% Inputs:
%   geometry   - see supergui()
%   listui     - list of uicontrol lists describing elements properties
%                { { ui1 }, { ui2 }... }, { 'uiX' } being GUI matlab 
%                uicontrol arguments such as { 'style', 'radiobutton', 
%                'String', 'hello' }. See supergui() for details.
%   help       - optional help command 
%   title      - optional figure title
%
% Output:
%   outparam   - list of outputs. The function scan all lines and
%                add up an output for each interactive uicontrol, i.e
%                edit box, radio button, checkbox and listbox.
%   userdat    - 'userdata' value of the figure.
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
% 02/15/02 add userdat option -ad
% 02/16/02 add figure title option -ad

function [result, userdat] = inputgui( geometry, listui, helpcom, mytitle);

if nargin < 2
   help inputgui;
   return;
end;	

fig = figure;
if exist('mytitle') == 1, set(gcf, 'name', mytitle); end;
geometry = { geometry{:} [1] [1 1 1] }; % add button to geometry

% add the three buttons
% ---------------------
listui = { listui{:}, {}, { 'Style', 'pushbutton', 'string', 'Cancel', 'callback', 'close gcbf' } };
if exist('helpcom') == 1
	listui = { listui{:}, { 'Style', 'pushbutton', 'string', 'Help', 'callback', helpcom } };
else
   listui = { listui{:}, {} };
end;   
listui = { listui{:}, { 'Style', 'pushbutton', 'tag', 'ok', 'string', 'OK', 'callback', 'set(gcbo, ''userdata'', ''retuninginputui'');' } };

% create figure and wait for return
% ---------------------------------
result = {};
[tmp tmp2 allobj] = supergui( geometry, listui{:} );

waitfor( findobj('parent', fig, 'tag', 'ok'), 'userdata');
%cont = 1;
%objok = findobj('parent', fig, 'tag', 'ok');
%while cont
%    pause(0.02);
%    try, cont = isempty(get(objok, 'userdata'));
%    catch, cont = 0; end;
%end;    

userdat = [];
try, findobj(fig); % figure still exist ?
catch, return; end;

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
close(fig);
