% windows() - Returns handle to window function or window
%
% Usage:
%   >> h = windows(t);
%   >> h = windows(t, m);
%   >> h = windows(t, m, a);
%
% Inputs:
%   t - char array 'rectangular', 'bartlett', 'hann', 'hamming',
%       'blackman', 'blackmanharris', or 'kaiser'
%
% Optional inputs:
%   m - scalar window length
%   a - scalar or vector with window parameter(s)
%
% Output:
%   h - function handle or column vector window
%
% Author: Andreas Widmann, University of Leipzig, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

function h = windows(t, m, a)

    if nargin < 1
        error('Not enough input arguments.');
    end
    h = str2func(t);
    switch nargin
        case 2
            h = h(m);
        case 3
            h = h(m, a);
    end
end
    
function w = rectangular(m)
    w = ones(m, 1);
end

function w = bartlett(m)
    w = 1 - abs(-1:2 / (m - 1):1)';
end

% von Hann
function w = hann(m);
    w = hamming(m, 0.5);
end

% Hamming
function w = hamming(m, a)
    if nargin < 2 || isempty(a)
        a = 25 / 46;
    end
    m = [0:1 / (m - 1):1]';
    w = a - (1 - a) * cos(2 * pi * m);
end

% Blackman
function w = blackman(m, a)
    if nargin < 2 || isempty(a)
        a = [0.42 0.5 0.08 0];
    end
    m = [0:1 / (m - 1):1]';
    w = a(1) - a(2) * cos (2 * pi * m) + a(3) * cos(4 * pi * m) - a(4) * cos(6 * pi * m);
end

% Blackman-Harris
function w = blackmanharris(m)
    w = blackman(m, [0.35875 0.48829 0.14128 0.01168]);
end

% Kaiser
function w = kaiser(m, a)
    if nargin < 2 || isempty(a)
        a = 0.5;
    end
    m = [-1:2 / (m - 1):1]';
    w = besseli(0, a * sqrt(1 - m.^2)) / besseli(0, a);
end

% Tukey    
function w = tukey(m, a)
    if nargin < 2 || isempty(a)
        a = 0.5;
    end
    if a <= 0
        w = ones(m, 1);
    elseif a >= 1
        w = hann(m);
    else
        a = (m - 1) / 2 * a;
        tapArray = (0:a)' / a;
        w = [0.5 - 0.5 * cos(pi * tapArray); ...
             ones(m - 2 * length(tapArray), 1); ...
             0.5 - 0.5 * cos(pi * tapArray(end:-1:1))];
    end
end
