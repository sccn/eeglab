% axcopy() - Copy a Matlab figure axis and its graphic objects to a new pop-up window 
%            using the left mouse button.
%
% Usage:  >> axcopy
%         >> axcopy(fig)
%         >> axcopy('noticks')
%         >> axcopy(fig, command)
%
% Notes:
%   1) Clicking the left mouse button on a Matlab figure axis copies the graphic objects 
%      in the axis to a new (pop-up) figure window. 
%   2) Option 'noticks' does not make x and y tickloabelmodes 'auto' in the pop-up.
%   2) The command option is an optional string that is evaluated in the new window
%      when it pops up. This allows the user to customize the pop-up display.
%   3) Deleting the pop-up window containing the copied axis leaves the selected axis
%      as the current graphic axis (gca).

% Authors: Tzyy-Ping Jung & Scott Makeig, SCCN/INC/UCSD, La Jolla, 2000 

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

if ~strcmpi(get(fig, 'type'), 'axes')
    hndl= findobj('parent',fig,'type','axes');
else
    hndl=fig;
end;
offidx=[];
if exist('command') ~= 1
    comstr = 'copyaxis';
else
   command_dbl = double(command);
   comstr = double(['copyaxis(''' char(command_dbl) ''')']); 
end;
for a=1:length(hndl)                    % make all axes visible
    allobjs = findobj('parent',hndl(a));
    for index = 1:length(allobjs)
        if isempty(get(allobjs(index), 'ButtonDownFcn'))
            set(allobjs(index), 'ButtonDownFcn', char(comstr));
        end;
    end;
end
if ~strcmpi(get(fig, 'type'), 'axes')
    figure(fig);
else
    figure(get(fig, 'parent'));
end;
if exist('command') ~= 1
    set(hndl(a),'ButtonDownFcn','copyaxis');
else
    % set(hndl,'ButtonDownFcn',['copyaxis(''' command ''')']);
    comstr = double(['copyaxis(''' char(command_dbl) ''')']);
    set(hndl,'ButtonDownFcn',char(comstr));
end;        
%set(hndl,'ButtonDownFcn','copyaxis');
%if ~exist('noticks')
  %axis on
  %set(gca,'xticklabelmode','auto')
  %set(gca,'yticklabelmode','auto')
%end
