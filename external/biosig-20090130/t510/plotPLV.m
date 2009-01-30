function plotPLV(r)
% Plots the PLV over the time (x axis) and frequency (y axis).
%
% This function plots the PLV values over time and frequency as calculated by
% calcPLV.m.
% 
% Usage:
%   plotPLV(r);
% 
% Input parameters:
%   r ... Input structure

% Copyright by Clemens Brunner
% $Revision: 1.1 $ $Date: 2009-01-30 06:04:51 $
% E-Mail: clemens.brunner@tugraz.at

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

plv = r{1}.plv;
sig = r{1}.sig;
fs = r{1}.fs;
f = r{1}.f;

plv = flipud(plv');
sig = flipud(sig');

if (isempty(sig))
    plot_signal = plv;
else
    plot_signal = plv.*sig;
end;

imagesc(plot_signal);
set(gca,'YTick',[f(2)-floor(f(2)/10)*10+1:10:f(2)-f(1)+1]);
set(gca,'YTickLabel',[floor(f(2)/10)*10:-10:ceil(f(1)/10)*10]);
% if (~isempty(sc))
%     caxis(sc);
% end;

if (~isempty(sig))
    colormap([[1 1 1];jet]);  % Draw all nans white
end;
    
line([fs*3,fs*3],[1,f(2)-f(1)+1],'Color','k');  % Mark trigger instant

set(gca,'XTick',[1 fs:fs:7*fs-1]);
set(gca,'XTickLabel',[-3:3]);
