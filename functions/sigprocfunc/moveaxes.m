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

function moveaxes(fig)

if nargin<1
  fig = gcf;
end

if ~ischar(fig)
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
   
   
