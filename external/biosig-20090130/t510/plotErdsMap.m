function plotErdsMap(r, plotstartend)
% Displays time-frequency (ERDS) maps.
%
% This function plots ERDS maps as calculated by calcErdsMaps.m.
%
% Usage:
%   plotErdsMap(r, plotstartend);
%
% Input parameters:
%   r ... ERDS map structure
%
% Optional input parameters:
%   plotstartend ... Start and end point in the plotted maps (s)

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

if (nargin < 1)
    error('No ERDS map specified.');
end;
if (nargin < 2)
    plotstartend = [];
end;

plot_index = r{1}.plot_index;
n_rows = r{1}.n_rows;
n_cols = r{1}.n_cols;
fs = r{1}.fs;
f = r{1}.f;
name = r{1}.name;
trials = r{1}.trials;
class = r{1}.class;
triallen = r{1}.triallen;
ref = r{1}.ref;
alpha = r{1}.alpha;

% Plot ERDS maps
fig = figure('Color', [1 1 1]);
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperType', 'A4');
set(gcf, 'PaperPosition', [1, 1, 27.7, 19]);
set(gcf, 'DefaultAxesFontSize', 8);

% Invert color map so that ERS is blue and ERD is red
c = flipud(jet(256));
if (~isnan(alpha))  % Plot non-significant values in white
    c(129,:) = [1 1 1];
end;
colormap(c);

for k = 1:length(plot_index)
    subplot(n_rows, n_cols, plot_index(k));
    
    if (~isempty(plotstartend))
        time_axis = r{k}.times(r{k}.times >= plotstartend(1) * 1000 & r{k}.times <= plotstartend(2) * 1000);
        
        start_time_index = find(r{k}.times==time_axis(1));
        stop_time_index = find(r{k}.times==time_axis(end));
        
        imagesc(time_axis, r{k}.freqs(r{k}.dispf), r{k}.PP(r{k}.dispf,start_time_index:stop_time_index), [-1.5 1.5]);
    else
        imagesc(r{k}.times, r{k}.freqs(r{k}.dispf), r{k}.PP(r{k}.dispf,:), [-1.5 1.5]);
    end;

    v = axis;
%     rectangle('Position', [r{1}.baseline(1) v(3) ...
%               abs(r{1}.baseline(1)-r{1}.baseline(2)) abs(v(4)-v(3))], ...
%               'FaceColor', [0.95 0.95 0.95], 'EdgeColor', 'none', ...
%               'EraseMode', 'xor');
%     temp = get(gca, 'Children');
%     set(gca, 'Children', temp([2:end 1]));  % Move rectangle to background
    line([r{1}.baseline(1) r{1}.baseline(1)], [v(3) v(4)], 'LineStyle', ':', 'Color', 'k');
    line([r{1}.baseline(2) r{1}.baseline(2)], [v(3) v(4)], 'LineStyle', ':', 'Color', 'k');

    if ~isempty(r{k}.cue)
        v = axis;
        line([r{k}.cue r{k}.cue], [v(3) v(4)], 'Color', 'k');
    end;

    set(gca, 'ydir', 'normal');
    set(gca, 'XTick', 0:1000:triallen*1000, 'FontSize', 8);
    set(gca, 'XTickLabel', 0:triallen, 'FontSize', 8);
    set(gca, 'YTick', 0:5:round(f(2)), 'FontSize', 8);
    set(gca, 'YTickLabel', 0:5:round(f(2)), 'FontSize', 8);
    set(gca, 'Tag', num2str(k));
    
    axis square;
end;

% Edit plot
if (isempty(ref))
    ref = [0 triallen/fs];
end;

[nrows, ncols] = size(class);
if (ncols == 1)
    class = class';
end;

annotation('textbox', [0 0 0.3 0.1], 'LineStyle', 'none', 'Interpreter', ...
           'none', 'String', {[name ', ' date], ['fs = ' num2str(fs) ...
           ' Hz, alpha = ' num2str(alpha) ', f = ' ...
           num2str(round(r{1}.freqs(r{1}.dispf(1))*10)/10) '-' ...
           num2str(round(r{1}.freqs(r{1}.dispf(end))*10)/10) ' Hz'], ...
           ['trials = ' num2str(trials) ', length = ' ...
           num2str(triallen/fs) ' s, reference = ' num2str(ref(1)) '-' ...
           num2str(ref(2)) ' s'], ['class = ' num2str(class)]});

% Enable interactively exploring ERDS maps:
% Clicking on a small map enlarges it in a new figure
% Clicking in this enlarged map opens a new figure displaying the time course in
% a specific frequency band
hndl = findobj('parent', fig, 'type', 'axes');
for a = 1:length(hndl)
    set(findobj('parent', hndl(a)), 'ButtonDownFcn', 'zoomMap(r)');
end;