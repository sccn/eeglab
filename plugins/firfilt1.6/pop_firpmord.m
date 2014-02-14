% pop_firpmord() - Estimate Parks-McClellan filter order and weights
%
% Usage:
%   >> [m, wtpass, wtstop] = pop_firpmord(f, a); % pop-up window mode
%   >> [m, wtpass, wtstop] = pop_firpmord(f, a, dev);
%   >> [m, wtpass, wtstop] = pop_firpmord(f, a, dev, fs);
%
% Inputs:
%   f         - vector frequency band edges
%   a         - vector desired amplitudes on bands defined by f
%   dev       - vector allowable deviations on bands defined by f
%
% Optional inputs:
%   fs        - scalar sampling frequency {default 2}
%
% Output:
%   m         - scalar estimated filter order
%   wtpass    - scalar passband weight
%   wtstop    - scalar stopband weight
%
% Note:
%   Requires the signal processing toolbox. Convert passband ripple from
%   dev to peak-to-peak dB: rp = 20 * log10((1 + dev) / (1 - dev)).
%   Convert stopband attenuation from dev to dB: rs = 20 * log10(dev).
%
% Author: Andreas Widmann, University of Leipzig, 2005
%
% See also:
%   pop_firpm, firpm, firpmord

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

function [m, wtpass, wtstop] = pop_firpmord(f, a, dev, fs)

m = [];
wtpass = [];
wtstop = [];

if exist('firpmord') ~= 2
   error('Requires the signal processing toolbox.');
end

if nargin < 2 || isempty(f) || isempty(a)
    error('Not enough input arguments');
end

% Sampling frequency
if nargin < 4 || isempty(fs)
    fs = 2;
end

% GUI
if nargin < 3 || isempty(dev)
    drawnow;
    uigeom = {[1 1] [1 1]};
    uilist = {{'style' 'text' 'string' 'Peak-to-peak passband ripple (dB):'} ...
              {'style' 'edit'} ...
              {'style' 'text' 'string' 'Stopband attenuation (dB):'} ...
              {'style' 'edit'}};
    result = inputgui(uigeom, uilist, 'pophelp(''pop_firpmord'')', 'Estimate filter order and weights -- pop_firpmord()');
    if length(result) == 0, return, end

    if ~isempty(result{1})
        rp = str2num(result{1});
        rp = (10^(rp / 20) - 1) / (10^(rp / 20) + 1);
        dev(find(a == 1)) = rp;
    else
        error('Not enough input arguments.');
    end
    if ~isempty(result{2})
        rs = str2num(result{2});
        rs = 10^(-abs(rs) / 20);
        dev(find(a == 0)) = rs;
    else
        error('Not enough input arguments.');
    end
end

[m, fo, ao, w] = firpmord(f, a, dev, fs);
wtpass = w(find(a == 1, 1));
wtstop = w(find(a == 0, 1));
