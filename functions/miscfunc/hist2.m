% hist2() - draw superimposed histograms
%
% Usage:
%   >> hist2(data1, data2);
%
% Inputs:
%   data1   - data to plot first process
%   data2   - data to plot second process
%
% Author: Arnaud Delorme (SCCN, UCSD)

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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


% draw superimposed histograms
% ---------------
function hist2(data1, data2, bins);

hist(data1, bins);
hold on; hist(data2, bins);
%figure; hist( [ measure{:,5} ], 20);
%hold on; hist([ measure{:,2} ], 20);
c = get(gca, 'children');
numfaces = size(get(c(1), 'Vertices'),1);
set(c(1), 'FaceVertexCData', repmat([1 0 0], [numfaces 1]), 'Cdatamapping', 'direct', 'facealpha', 0.5, 'edgecolor', 'none');
numfaces = size(get(c(2), 'Vertices'),1);
set(c(2), 'FaceVertexCData', repmat([0 0 1], [numfaces 1]), 'Cdatamapping', 'direct', 'facealpha', 0.5, 'edgecolor', 'none');
ylabel('Number of values');
