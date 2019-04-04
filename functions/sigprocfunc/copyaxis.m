% copyaxis() - helper function for axcopy()
%
% Usage: >> copyaxis();
%        >> copyaxis( command );
% 
% Note: The optional command option (a string that will be evaluated
%       when the figure is created allows to customize display).
%
% Author: Tzyy-Ping Jung, SCCN/INC/UCSD, La Jolla, 2000 
%
% See also: axcopy()

% Copyright (C) 3-16-2000 Tzyy-Ping Jung, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 01-25-02 reformated help & license -ad 
% 02-16-02 added the ignore parent size option -ad 
% 03-11-02 remove size option and added command option -ad 

function copyaxis(command)

OFF_STD = 0.25; % std dev of figure offset
MIN_OFF = 0.15; % minimum offset for new figure
BORDER  = 0.04;  % screen edge tolerance

fig=gcf;
sel=gca;
%
% Get position for new figure
%
set(fig,'Units','normalized');
place=get(fig,'Position');
axiscolor=get(fig,'Color');
cmap=colormap;
newxy = (OFF_STD*randn(1,2))+place(1,1:2);
newx=newxy(1);newy=newxy(2);
if abs(newx-place(1,1))<MIN_OFF, newx=place(1,1)+sign(newx-place(1,1))*MIN_OFF;end
if abs(newy-place(1,1))<MIN_OFF, newy=place(1,1)+sign(newy-place(1,1))*MIN_OFF;end
if newx<BORDER, newx=BORDER; end
if newy<BORDER, newy=BORDER; end
if newx+place(3)>1-BORDER, newx=1-BORDER-place(3); end
if newy+place(4)>1-BORDER, newy=1-BORDER-place(4); end

newfig=figure('Units','Normalized','Position',[newx,newy,place(1,3:4)]);

%
% Copy object to new figure
%
set(newfig,'Color',axiscolor);
copyobj(sel,newfig);set(gca,'Position',[0.130 0.110 0.775 0.815]);
colormap(cmap);
%
% Increase font size
%
set(findobj('parent',newfig,'type','axes'),'FontSize',14);
set(get(gca,'XLabel'),'FontSize',16)
set(get(gca,'YLabel'),'FontSize',16)
set(get(gca,'Title'),'Fontsize',16);
%
% Add xtick and ytick labels if missing
%
if strcmp(get(gca,'Box'),'on')
   set(gca,'xticklabelmode','auto')
   set(gca,'xtickmode','auto')
   set(gca,'yticklabelmode','auto')
   set(gca,'ytickmode','auto')
end
%
% Turn on zoom in the new figure
%
zoom on;
%
% Turn off the axis copy for the new figure
%
hndl= findobj('parent',newfig,'type','axes');
set(hndl,'ButtonDownFcn','','Selected','off');
for a=1:length(hndl) % make all axes visible
  set(findobj('parent',hndl(a)),'ButtonDownFcn','','Selected','off');
end
%
% Execute additional command if present
%
if exist('command') == 1
   eval(command);
end;   
