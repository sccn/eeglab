function zoomErds(r)
% Zooms out specific ERDS time course.
%
% This function is used to interactively click inside enlarged ERDS maps and
% creates a new figure containing the selected ERDS time course.

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

set(gca, 'Units', 'normalized');
mouse_pos = get(gca, 'CurrentPoint');
temp = round(mouse_pos(3));
curr_freq = find(abs(r{1}.freqs-temp) == min(abs(r{1}.freqs-temp)));

axiscolor = get(gcf, 'Color');  % Background color

nr = str2num(get(gca, 'Tag'));

newfig = figure;
set(newfig, 'Color', axiscolor);
plot(r{nr}.times, r{nr}.P(curr_freq,:));
hold on;
line([0 r{nr}.times(end)], [r{nr}.sig(1) r{nr}.sig(1)], 'LineStyle', ':', 'Color', 'k');
line([0 r{nr}.times(end)], [r{nr}.sig(2) r{nr}.sig(2)], 'LineStyle', ':', 'Color', 'k');

set(gca, 'XTick', 0:1000:8000, 'FontSize', 8);
set(gca, 'XTickLabel', 0:8, 'FontSize', 8);
v = axis;
axis([r{nr}.times(1) r{nr}.times(end) v(3) v(4)]);

set(findobj('parent', newfig, 'type', 'axes'), 'FontSize', 12);

rectangle('Position', [r{1}.baseline(1) v(3) abs(r{1}.baseline(1)-r{1}.baseline(2)) abs(v(4)-v(3))], 'FaceColor', [0.95 0.95 0.95], 'EdgeColor', 'none', 'EraseMode', 'xor');
temp = get(gca, 'Children');
set(gca, 'Children', temp([2:end 1]));  % Move rectangle to background

line([r{1}.baseline(1) r{1}.baseline(1)], [v(3) v(4)], 'LineStyle', ':', 'Color', 'k');
line([r{1}.baseline(2) r{1}.baseline(2)], [v(3) v(4)], 'LineStyle', ':', 'Color', 'k');

if ~isempty(r{1}.cue)
    v = axis;
    line([r{1}.cue r{1}.cue], [v(3) v(4)], 'Color', 'k');
end;

% Turn off interactive clicking
hndl = findobj('parent', newfig, 'type', 'axes');
for a = 1:length(hndl)
    set(findobj('parent', hndl(a)), 'ButtonDownFcn', '');
end;

% title(['Channel ' num2str(nr) ', f = ' num2str(round(r{1}.freqs(curr_freq)*10)/10) ' Hz']);
title(['Channel ' num2str(nr) ', f = ' num2str(round(r{1}.freqs(curr_freq)*10)/10) ' Hz, Pref = ' num2str(round(r{nr}.Pref(curr_freq)*10)/10)]);