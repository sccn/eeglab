% tree() - Make a hierarchical (tree-diagram) component plot. Use
%          successive calls to this function to build the full plot.
%
% Usage:
%  >> tree(data,channel,frames,xvals,offset,scaleflag,xlimits)
%
% Inputs:
%   data      = (ncomponents,frames*epochs) data waveforms 
%   channel   = data channel/component activation (row) to plot from
%   offset    = y-value offset for this waveform
%   frames    = time frames per epoch (0-> length(data))
%   xvals     = x-indices to plot (e.g. [40:70]) ( 0-> all frames)
%   scaleflag = (1 = scale data: [min,max] -> [-0.5,0.5])
%                   (default 0 = use actual y-values)
%   xlimits   = [actual_xmin xmax]
%
% Authors: Tzyy-Ping Jung & Scott Makeig, SCCN/INC/UCSD, La Jolla, 11-30-96 

% Copyright (C) 11-30-96 Tzyy-Ping Jung & Scott Makeig, SCCN/INC/UCSD, 
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

% 01-25-02 reformated help & license -ad 

function tree(data,channel,frames,xvals,offset,scaleflag,xlimits)

  [chans,totalframes]=size(data);

  errorval = 0;
  if nargin<6,
	scaleflag = 0;
  end
  if nargin==3,
	frames = totalframes;
  end
  if nargin==4,
    xvals = [1:frames];
  end
  if nargin == 6,
  	xmin = xlimits(1);
  	xmax = xlimits(2);
  else
  	xmin = xlimits(1);
  	xmax = xlimits(2);
  end

	if xvals(length(xvals))<xvals(1),
		fprintf('tree(): xvals must be in ascending order!\n');
        errorval = 1;
	end
  if errorval == 0,
    mindata=min(min(data));
    maxdata=max(max(data));
    if scaleflag == 1,
  	  set(0,'DefaultAxesYLim',[0 chans])	% plotting parameters
    else
    %	set(0,'DefaultAxesYLim',[mindata maxdata])	% plotting parameters
    end
    set(0,'DefaultAxesXLim',[0 frames])
    set(0,'DefaultAxesFontSize',18);
    set(0,'DefaultTextFontSize',18);
  
    curfig = gcf;
    h=figure(curfig);
    set(h,'PaperPosition',[0.2 0.3 7.6 10]); % make figure full page
  
    % colors =['r';'b';'w';'g';'c';'m';'r';'b';'w';'g';'c';'m';'r';'b';'w';'g';'c';'m';'r';'b';'w';'g';'c';'m'];

 colors =['r';'b';'g';'c';'m';'w';'r';'b';'g';'c';'m';'w';'r';'b';'g';'c';'m';'w';'r';'b';'g';'c';'m';'w'];
    xdata = [xmin:(xmax-xmin)/(frames-1):xmax];
  
    msg = [int2str(channel)];
    text(xdata(xvals(length(xvals)))+10,offset,msg);
  
    for i=1:fix(totalframes/frames) % = epochs
      x=(i-1)*frames;
      hold on;
      tmp=data(channel,x+xvals);
	  if scaleflag == 1,
    	  tmp=(tmp-mindata)/(maxdata-mindata) -0.5;
	  end
      plot(xdata(xvals),tmp+offset,colors(i))
    end
    tmp = zeros(1,length(xvals));
    plot(xdata(xvals),tmp+offset,'w');
  end % if errorval...
