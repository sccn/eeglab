% axcopy() - move, resize, or copy Matlab axes using the mouse
%
% Usage:  >> selaxiscopy
%         >> selaxiscopy(fig)
%         >> selaxiscopy('noticks')
%         >> selaxiscopy(fig, command)
%
% Notes:
%   1) Clicking the left mouse button on an axis and copy
%   the objects in the axis to a new (popup) figure. 
%   Option 'noticks' does not make x and y tickloabelmodes 'auto'
%   in the popup.
%   2) The additional command option (a string that will be evaluated
%   when the figure pops-up allows to customize display).
%
% Authors: Tzyy-Ping Jung & Scott Makeig, SCCN/INC/UCSD, La Jolla, 2000 

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2000 Tzyy-Ping Jung & Scott Makeig, SCCN/INC/UCSD, 3-16-00 
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

% requires copyaxes.m
% 4-2-00 added 'noticks' -sm
% 01-25-02 reformated help & license -ad 
% 02-16-02 debuged exist('fig') -> exist('fig') == 1 -ad 
% 02-16-02 added the ignore parent size option -ad 
% 03-11-02 remove size option and added command option -ad 

function axcopy(fig, command)

if (exist('fig') == 1) & strcmp(fig,'noticks')
   noticks = 1;
   if nargin> 1
     shift
   else
     fig = [];
   end
end
if ~(exist('fig') ==1) | isempty(fig) | fig == 0 
   fig = gcf;
end

hndl= findobj('parent',fig,'type','axes');
offidx=[];
for a=1:length(hndl)                    % make all axes visible
    if exist('command') ~= 1
        set(findobj('parent',hndl(a)),'ButtonDownFcn','copyaxis');
    else
        set(findobj('parent',hndl(a)),'ButtonDownFcn',['copyaxis(''' command ''')']);
    end;        
end
figure(fig);
if exist('command') ~= 1
    set(hndl(a),'ButtonDownFcn','copyaxis');
else
    set(hndl,'ButtonDownFcn',['copyaxis(''' command ''')']);
end;        
%set(hndl,'ButtonDownFcn','copyaxis');
%if ~exist('noticks')
  %axis on
  %set(gca,'xticklabelmode','auto')
  %set(gca,'yticklabelmode','auto')
%end
