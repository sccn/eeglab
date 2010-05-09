% moveaxes() - move, resize, or copy Matlab axes using the mouse
%
% Usage: >> moveaxes
%        >> moveaxes(fig)
%        >> moveaxes off 
%
% Note: clicking the left mouse button selects an axis
%       dragging the left mouse button resizes a selected axis
%       dragging the right mouse button copies a selected axis
%       clicking the middle mouse button deselects a selected axis
%
% Authors: Colin Humphries & Scott Makeig, CNL / Salk Institute, La Jolla

% Copyright (C) Colin Humphries & Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 01-25-02 reformated help & license -ad 

function moveaxes(fig)

if nargin<1
  fig = gcf;
end

if ~isstr(fig)
   tmp = findobj('tag','moveaxes');
   if ~isempty(tmp)   % turn off previous moveaxes Tags 
      ax = get(tmp,'children');
      set(ax,'ButtonDownFcn','','Selected','Off');
      set(tmp,'Tag','','UserData',[]);
   end

   ax = findobj('parent',fig,'type','axes');
   for a=1:length(ax)                    % make all axes visible
     axvis = get(ax(a),'visible');
     set(ax(a),'visible','on','UserData',axvis);
   end

   set(ax,'ButtonDownFcn','selectmoveresize');
   set(fig,'UserData',ax,'Tag','moveaxes');

elseif strcmp(fig,'off')
   fig=findobj('Tag','moveaxes');
   ax = get(fig,'UserData');
   for a=1:length(ax)                    % reset original axis visibility
     axvis= get(ax(a),'UserData')
     set(ax(a),'visible',axvis);
   end
   set(ax,'ButtonDownFcn','','Selected','off');
   set(fig,'Tag','','UserData',[]);
end
   
   
