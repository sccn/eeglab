% ploterp() - plot a selected multichannel data epoch on a common timebase
%
% Usage: >> ploterp(data,frames,epoch,[limits],'title',[plotchans]);
%
% Inputs:
%  data       = EEG/ERP data epoch (chans,frames)
%  frames     = frames per epoch {default: data length}
%  epoch      = epoch to plot {default: 1}
%  [limits]   = [xmin xmax ymin ymax]  (x's in ms) 
%                   {def|0 or both y's 0 -> data limits}
% 'title'     = plot title {default|0 -> none}
%  plotchans  = data channels to plot {default|0 -> all}
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 6-11-98 
%
% See also: timtopo()

% Copyright (C) 6-11-98 from plotdata() Scott Makeig, SCCN/INC/UCSD,
% scott@sccn.ucsd.edu
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

% 01-25-02 reformated help & license, added links -ad 

function plot_handl = ploterp(data,frames,epoch,limits,titl,plotchans)

LABELFONT = 16;
TICKFONT = 14;
TITLEFONT = 16;

if nargin < 1
   help ploterp
   return
end

[chans,framestot] = size(data);
icadefs;   

if nargin < 6
   plotchans = 0;
end
if plotchans==0
   plotchans = 1:chans;
end

if nargin < 5,
   titl = '';     % DEFAULT TITLE
end

if nargin < 4,
    limits = 0;
end

if nargin < 3
  epoch = 0;
end
if epoch == 0
  epoch = 1;
end

if nargin<2
 frames = 0;
end
if frames == 0
 frames = size(data,2);
end
if floor(framestot/frames)*frames ~= framestot
  fprintf('ploterp(): frames argument does not divide data length.\n');
  return
end
if epoch*frames > framestot
  fprintf('ploterp(): data does not contain %d epochs of %d frames.\n',epoch,frames);
  return
end

%
%%%%%%%%%%%%%%%%%%%%%%% Read and adjust limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  if limits==0,      % == 0 or [0 0 0 0]
    xmin=0;
    xmax=frames-1;
    ymin=min(min(data));
    ymax=max(max(data));
  else
    if length(limits)~=4,
      fprintf( ...
       'ploterp: limits should be 0 or an array [xmin xmax ymin ymax].\n');
      return
    end;
    xmin = limits(1);
    xmax = limits(2);
    ymin = limits(3);
    ymax = limits(4);
  end;

  if xmax == 0 & xmin == 0,
    x = (0:1:frames-1);
    xmin = 0;
    xmax = frames-1;
  else
    dx = (xmax-xmin)/(frames-1);
    x=xmin*ones(1,frames)+dx*(0:frames-1); % compute x-values
  end;
  if xmax<=xmin,
      fprintf('ploterp() - xmax must be > xmin.\n')
      return
  end

  if ymax == 0 & ymin == 0,
      ymax=max(max(data));
      ymin=min(min(data));
  end
  if ymax<=ymin,
      fprintf('ploterp() - ymax must be > ymin.\n')
      return
  end

sampint = (xmax-xmin)/(frames-1); % sampling interval = 1000/srate;
x = xmin:sampint:xmax;   % make vector of x-values

%
%%%%%%%%%%%%%%%%%%%%%%% Plot the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
cla % clear the current figure axes
set(gca,'YColor',BACKCOLOR); % set the background color
set(gca,'Color',BACKCOLOR);

set(gca,'GridLineStyle',':')
set(gca,'Xgrid','off')
set(gca,'Ygrid','on')
set(gca,'Color',BACKCOLOR,'FontSize',TICKFONT,'FontWeight','bold');
plot_handl=plot(x,data(plotchans,(epoch-1)*frames+1:epoch*frames));    % plot the data
title(titl,'fontsize',TITLEFONT,'FontWeight','bold');

l= xlabel('Time (ms)');
set(l,'FontSize',LABELFONT,'FontWeight','bold');
l=ylabel('Potential (uV)');
set(l,'FontSize',LABELFONT,'FontWeight','bold');
axis([xmin xmax ymin ymax]);
