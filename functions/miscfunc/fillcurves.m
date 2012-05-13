% fillcurves() - fill the space between 2 curves
%
% Usage:
%   h=fillcurves( Y1, Y2);
%   h=fillcurves( X, Y1, Y2);
%
% Example:
%   a = rand(1, 50);
%   b = rand(1, 50)+2; b(10) = NaN;
%   figure; fillcurves([51:100], a, b);
%
% Author: A. Delorme, SCCN, INC, UCSD/CERCO, CNRS

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function h = fillcurves(X, Y1, Y2, color);
    
    if nargin < 2
        help fillcurves;
        return;
    end;
    
    if nargin < 3
        Y2 = Y1;
        Y1 = X;
        X = [1:length(Y1)];
    end;
    if nargin < 4
        color = 'r';
    end;
    X1 = X(:)';
    X2 = X(:)';
    
    % remove NaNs
    tmpnan1 = find(isnan(Y1));
    tmpnan2 = find(isnan(Y2));
    Y1(tmpnan1) = []; X1(tmpnan1) = [];
    Y2(tmpnan2) = []; X2(tmpnan2) = [];
    
    % remove 0s if log scale
    if strcmpi(get(gca, 'yscale'), 'log')
        tmp1 = find(~Y1);
        tmp2 = find(~Y2);
        Y1(tmp1) = []; X1(tmp1) = [];
        Y2(tmp2) = []; X2(tmp2) = [];
    end;
    
    % plot
    allpointsx = [X1 X2(end:-1:1)]';
    allpointsy = [Y1 Y2(end:-1:1)]';
    h = fill(allpointsx, allpointsy, color); 
    set(h, 'edgecolor', 'r');
    xlim([ X(1) X(end) ]);
    
