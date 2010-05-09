% cbar() - Display full or partial color bar
%
% Usage:
%    >> cbar % create a vertical cbar on the right side of a figure
%    >> cbar(type) % specify direction as 'vert' or 'horiz'
%    >> cbar(type,colors) % specify which colormap colors to plot
%  else
%    >> cbar(axhandle) % specify the axes to draw cbar in
%
%    >> h = cbar(type|axhandle,colors, minmax, grad)
%
% Inputs:
%  type      - ['vert'|'horiz'] direction of the cbar {default: 'vert')
%              ELSE axhandle = handle of axes to draw the cbar
%  colors    - vector of colormap indices to display, or integer to truncate upper 
%              limit by.
%              (int n -> display colors [1:end-n]) {default: 0}
%  minmax    - [min, max] range of values to label on colorbar 
%  grad      - [integer] number of tick labels. {default: 5}.
%
% Example:
%         >> colormap('default') % default colormap is 64-color 'jet'
%         >> cbar('vert',33:64); % plot a vertical cbar colored green->red 
%                                % useful for showing >0 (warm) and 0 (green) 
%                                % values only in a green=0 plot
%
% Author: Colin Humphries, Arnaud Delorme, CNL / Salk Institute, Feb. 1998-
%
% See also: colorbar()

% Copyright (C) Colin Humphries, CNL / Salk Institute, Feb. 1998 
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

% 12-13-98 added minmax arg -Scott Makeig
% 01-25-02 reformated help & license, added links -ad 

function [handle]=cbar(arg,colors,minmax, grad)

if nargin < 2
  colors = 0;
end
posscale = 'off';
if nargin < 1
  arg = 'vert';
  ax = [];
else
  if isempty(arg)
    arg = 0;
  end
  if arg(1) == 0
    ax = [];
    arg = 'vert';
  elseif strcmpi(arg, 'pos')
    ax = [];
    arg = 'vert';
    posscale = 'on';
  else      
    if ischar(arg)
      ax = [];
    else
      ax = arg;
      arg = [];
    end
  end
end

if nargin>2
  if size(minmax,1) ~= 1 | size(minmax,2) ~= 2
    help cbar
    fprintf('cbar() : minmax arg must be [min,max]\n');
    return
  end
end
if nargin < 4
    grad = 5;
end;

%obj = findobj('tag','cbar','parent',gcf);
%if ~isempty(obj) & ~isempty(arg)
%  arg = [];
%  ax = obj;
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose colorbar position
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (length(colors) == 1) & (colors == 0)
  t = caxis;
else
  t = [0 1];
end
if ~isempty(arg)
  if strcmp(arg,'vert')  
    cax = gca;
    pos = get(cax,'Position');
    stripe = 0.04; 
    edge = 0.01;
    space = .02;

%    set(cax,'Position',[pos(1) pos(2) pos(3)*(1-stripe-edge-space) pos(4)])
%    rect = [pos(1)+(1-stripe-edge)*pos(3) pos(2) stripe*pos(3) pos(4)];

    set(cax,'Position',[pos(1) pos(2) pos(3) pos(4)])
    rect = [pos(1)+pos(3)+space pos(2) stripe*pos(3) pos(4)];
    ax = axes('Position', rect);
  elseif strcmp(arg,'horiz')
    cax = gca;
    pos = get(cax,'Position');
    stripe = 0.075; 
    space = .1;  
    set(cax,'Position',...
        [pos(1) pos(2)+(stripe+space)*pos(4) pos(3) (1-stripe-space)*pos(4)])
    rect = [pos(1) pos(2) pos(3) stripe*pos(4)];
    ax = axes('Position', rect);
  end
else
  pos = get(ax,'Position');
  if pos(3) > pos(4)
    arg = 'horiz';
  else
    arg = 'vert';
  end
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw colorbar using image()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

map = colormap;
n = size(map,1);

if length(colors) == 1
  if strcmp(arg,'vert')
      if strcmpi(posscale, 'on')
          image([0 1],[0 t(2)],[ceil(n/2):n-colors]');
      else
          image([0 1],t,[1:n-colors]');
      end;
      set(ax,'xticklabelmode','manual')
      set(ax,'xticklabel',[],'YAxisLocation','right')
      
  else
    image(t,[0 1],[1:n-colors]);
    set(ax,'yticklabelmode','manual')
    set(ax,'yticklabel',[],'YAxisLocation','right')
  end
  set(ax,'Ydir','normal','YAxisLocation','right')

else % length > 1

  if max(colors) > n
    error('Color vector excedes size of colormap')
  end
  if strcmp(arg,'vert')
    image([0 1],t,[colors]');
    set(ax,'xticklabelmode','manual')
    set(ax,'xticklabel',[])
  else
    image([0 1],t,[colors]);
    set(ax,'yticklabelmode','manual')
    set(ax,'yticklabel',[],'YAxisLocation','right')
  end  
  set(ax,'Ydir','normal','YAxisLocation','right')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust cbar ticklabels
%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin > 2 
  if strcmp(arg,'vert')
      Cax = get(ax,'Ylim');
  else
      Cax = get(ax,'Xlim');
  end;
  CBTicks = [Cax(1):(Cax(2)-Cax(1))/(grad-1):Cax(2)]; % caxis tick positions
  CBLabels = [minmax(1):(minmax(2)-minmax(1))/(grad-1):minmax(2)]; % tick labels
  
  dec = floor(log10(max(abs(minmax)))); % decade of largest abs value
  CBLabels = ([minmax]* [ linspace(1,0, grad);linspace(0, 1, grad)]);
  %[1.0 .75 .50 .25 0.0; 0.0 .25 .50 .75 1.0]);
  if dec<1
    CBLabels = round(CBLabels*10^(1-dec))*10^(dec-1);
  elseif dec == 1
    CBLabels = round(CBLabels*10^(2-dec))*10^(dec-2);
  else
    CBLabels = round(CBLabels);
  end
% minmax
% CBTicks
% CBLabels

  if strcmp(arg,'vert')
    set(ax,'Ytick',CBTicks);
    set(ax,'Yticklabel',CBLabels);
  else
    set(ax,'Xtick',CBTicks);
    set(ax,'Xticklabel',CBLabels);
  end
end
handle = ax;

%%%%%%%%%%%%%%%%%%
% Adjust cbar tag
%%%%%%%%%%%%%%%%%%

set(ax,'tag','cbar')
