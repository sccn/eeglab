% hist2() - draw superimposed histograms
%
% Usage:
%   >> hist2(data1, data2, bins);
%
% Inputs:
%   data1   - data to plot first process
%   data2   - data to plot second process
%
% Optional inputs:
%   bins    - vector of bin center
%
% Author: Arnaud Delorme (SCCN, UCSD)

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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


% draw superimposed histograms
% ---------------
function hist2(data1, data2, bins);

if nargin < 1
    help hist2;
    return;
end
if nargin < 3
    bins = linspace(min(min(data1), min(data2)), max(max(data1), max(data2)), 100);
elseif length(bins) == 1
    bins = linspace(min(min(data1), min(data2)), max(max(data1), max(data2)), bins);
end

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
xlim([bins(1) bins(end)]);

yl = ylim;
xl = xlim;
line([xl(1) xl(1)]+(xl(2)-xl(1))/2000, yl, 'color', 'k');
line(xl, [yl(1) yl(1)]+(yl(2)-yl(1))/2000, 'color', 'k');


