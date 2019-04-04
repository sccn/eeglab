% changeunits() - Takes one or more points in one axes and gives its position 
%                 in another axes. Useful for drawing lines between 
%                 sbplots (see sbplot()).
%
%  Usage: >> newpoint(s) = changeunits(point(s),curaxhandle,newaxhandle)
%
%  Inputs:
%     point(s)    - two-column vector of current [x y] data point locations 
%     curaxhandle - the point's current axes
%     newaxhandle - the new axes handle. If figure handle or absent, returns 
%                   the coordinates of the points in the whole figure 
%
%  Output:  
%     newpoint(s) - the coordinates of the same point(s) in the new axes
%
%  Example: 
%  >> figure
%  >> small1 = sbplot(4,4,10); % small axes in lower left
%  >> plot(0.3,0.4,'ro'); % mark point [0.3 0.4] in small1 axes
%  >> axis([0 1 0 1]); % set axes limits
%    
%  >> small2 = sbplot(4,4,7);  % small axes in upper right
%  >> plot(0.6,0.7,'ro'); % mark point [0.6 0.7] in small2 axes
%  >> axis([0 1 0 1]); % set axes limits
%    
%  >> large = sbplot(1,1,1);   % normal whole figure axes
%  >> % Now draw line from point [0.3 0.4] in small1 axes
%  >> %                 to point [0.6 0.7] in small2 axes
%  >> from = changeunits([0.3 0.4],small1,large); % point in small1 axes
%  >> to   = changeunits([0.6 0.7],small2,large); % point in small2 axes
%  >> plot([from(1) to(1)],[from(2) to(2)])
%  >> axis([0 1 0 1]); % set large axes limits
%  >> axis off % finally, hide large axes 
%
% Author: Colin Humphries, Salk Institute, La Jolla CA, Jan. 1998

% Copyright (C) 1998 Colin Humphries, Salk Institute, La Jolla CA, Jan. 1998
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

% 02-01-01 added default 3rd arg, multiple points, example -Scott Makeig
% 01-25-02 reformated help & license -ad 

function [newpnts] = changeunits(pnts,ax1,ax2)

curax = gca;
curfig = gcf;

if nargin<2
  help changeunits
  error('Must have at least 2 arguments');
end
if size(pnts,2) ~= 2
  help changeunits
  error('Points must be 2-column');
end

if ~strcmp(get(ax1,'type'),'axes')
  error('Second argument must be an axes handle')
end

if nargin<3
  figh = get(ax1,'Parent');
  figure(figh);
  ax2 = axes('Units','Normal','Position',[0 0 1 1],'Visible','Off');  
  % whole figure axes in ax1 figure
end


if strcmp(get(ax2,'type'),'axes')

   figh = get(ax1,'Parent');
   if figh ~= get(ax2,'Parent')
     error('Axes must be in the same figure.')
   end
   
   units1 = get(ax1,'Units');
   units2 = get(ax2,'Units');
   set(ax1,'Units','normalized')
   set(ax2,'Units','normalized')
   
   axpos1 = get(ax1,'Position');
   axpos2 = get(ax2,'Position');
   xlim1 = get(ax1,'Xlim');
   xlim2 = get(ax2,'Xlim');
   ylim1 = get(ax1,'Ylim');
   ylim2 = get(ax2,'Ylim');

   l1 = [xlim1(1) ylim1(1)];
   l2 = [xlim1(2) ylim1(2)];
   p1 = axpos1([1,2]);
   p2 = axpos1([3,4]);

   ll1 = [xlim2(1) ylim2(1)];
   ll2 = [xlim2(2) ylim2(2)];
   pp1 = axpos2([1,2]);
   pp2 = axpos2([3,4]);

   newpnts = zeros(size(pnts));
   for p = 1:size(pnts,1) 
     figpnt = (((pnts(p,:)-l1)./(l2-l1)).*p2) + p1;
     newpnts(p,:) = (((figpnt-pp1)./pp2).*(ll2-ll1)) + ll1;
   end

   set(ax1,'Units',units1)
   set(ax2,'Units',units2)
  
elseif strcmp(get(ax2,'type'),'figure')

  axpos1 = get(ax1,'Position');
  xlim1 = get(ax1,'Xlim'); % value limits in ax1
  ylim1 = get(ax1,'Ylim');

  l1 = [xlim1(1) ylim1(1)];
  l2 = [xlim1(2) ylim1(2)];
  p1 = axpos1([1,2]);
  p2 = axpos1([3,4]);

  newpnts = zeros(size(pnts));
  for p = 1:size(pnts,1) % unnecessary loop?
    newpnts(p,:) = (((pnts(p,:)-l1)./(l2-l1)).*p2) + p1;
  end

end 

if nargin<3
  delete(ax2)
end

figure(curfig); % restore gcf
axes(curax);    % restore gca
