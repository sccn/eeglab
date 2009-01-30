function zoomMap(r)
% Zooms out a specific channel within an ERDS map.
%
% This function is used to interactively click inside an ERDS maps and
% creates a new map containing the selected channel.

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

sel = gca;

% Get position for new figure
% set(gcf, 'Units', 'normalized');
% place = get(gcf, 'Position');
axiscolor = get(gcf, 'Color');  % Background color
cmap = colormap;

nr = str2num(get(gca, 'Tag'));

newfig = figure;

% Copy object to new figure
set(newfig, 'Color', axiscolor);
copyobj(sel, newfig);
set(gca, 'Position', [0.130 0.110 0.775 0.815]);
colormap(cmap);

% Increase font size
set(findobj('parent', newfig, 'type', 'axes'), 'FontSize', 12);

% Clicking in the map displays the corresponding ERDS curve
hndl = findobj('parent', newfig, 'type', 'axes');
for a = 1:length(hndl)
    set(findobj('parent', hndl(a)), 'ButtonDownFcn', 'zoomErds(r)');
end;

title(['Channel ' num2str(nr)]);