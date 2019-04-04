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

% requires copyaxes.m
% 4-2-00 added 'noticks' -sm
% 01-25-02 reformated help & license -ad 
% 02-16-02 debuged exist('fig') -> exist('fig') == 1 -ad 
% 02-16-02 added the ignore parent size option -ad 
% 03-11-02 remove size option and added command option -ad 

function axcopy(fig, command)

if (exist('fig') == 1) && strcmp(fig,'noticks')
   noticks = 1;
   if nargin> 1
     shift
   else
     fig = [];
   end
end
if ~(exist('fig') ==1) || isempty(fig) || fig == 0 
   fig = gcf;
end

if ~strcmpi(get(fig, 'type'), 'axes')
    %hndl= findobj('parent',fig,'type','axes');
    %--
    childs = get(fig,'Children');
    hndl   = childs(find(strcmpi(get(childs,'Type'),'axes')));
else
    hndl=fig;
end
offidx=[];
if exist('command') ~= 1
    comstr = 'copyaxis';
else
   command_dbl = double(command);
   comstr = double(['copyaxis(''' char(command_dbl) ''')']); 
end
for a=1:length(hndl)                    % make all axes visible
    %allobjs = findobj('parent',hndl(a));
    %--
    allobjs = get(hndl(a),'Children');
    for index = 1:length(allobjs)
        if isempty(get(allobjs(index), 'ButtonDownFcn'))
            set(allobjs(index), 'ButtonDownFcn', char(comstr));
        end
    end
end
if ~strcmpi(get(fig, 'type'), 'axes')
    figure(fig);
else
    figure(get(fig, 'parent'));
end
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
