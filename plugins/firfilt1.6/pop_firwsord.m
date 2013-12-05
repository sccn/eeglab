% pop_firwsord() - Estimate windowed sinc filter order depending on
%                  window type and requested transition band width
%
% Usage:
%   >> [m, dev] = pop_firwsord; % pop-up window mode
%   >> m = pop_firwsord(wtype, fs, df);
%   >> m = pop_firwsord('kaiser', fs, df, dev);
%
% Inputs:
%   wtype - char array window type. 'rectangular', 'bartlett', 'hann',
%           'hamming', {'blackman'}, or 'kaiser'
%   fs    - scalar sampling frequency {default 2}
%   df    - scalar requested transition band width
%   dev   - scalar maximum passband deviation/ripple (Kaiser window
%           only)
%
% Output:
%   m     - scalar estimated filter order
%   dev   - scalar maximum passband deviation/ripple
%
% References:
%   [1] Smith, S. W. (1999). The scientist and engineer's guide to
%       digital signal processing (2nd ed.). San Diego, CA: California
%       Technical Publishing.
%   [2] Proakis, J. G., & Manolakis, D. G. (1996). Digital Signal
%       Processing: Principles, Algorithms, and Applications (3rd ed.).
%       Englewood Cliffs, NJ: Prentice-Hall
%   [3] Ifeachor E. C., & Jervis B. W. (1993). Digital Signal
%       Processing: A Practical Approach. Wokingham, UK: Addison-Wesley
%
% Author: Andreas Widmann, University of Leipzig, 2005
%
% See also:
%   pop_firws, firws, pop_kaiserbeta, windows

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

function [m, dev] = pop_firwsord(wtype, fs, df, dev)

    m = [];

    wtypes = {'rectangular' 'bartlett' 'hann' 'hamming' 'blackman' 'kaiser'};

    % Window type
    if nargin < 1 || isempty(wtype)
        wtype = 5;
    elseif ~ischar(wtype) || isempty(strmatch(wtype, wtypes))
        error('Unknown window type');
    else
        wtype = strmatch(wtype, wtypes);
    end

    % Sampling frequency
    if nargin < 2 || isempty(fs)
        fs = 2;
    end

    % Transition band width
    if nargin < 3
        df = [];
    end

    % Maximum passband deviation/ripple
    if nargin < 4 || isempty(dev)
        devs = {0.089 0.056 0.0063 0.0022 0.0002 []};
        dev = devs{wtype};
    end

    % GUI
    if nargin < 3 || isempty(df) || (wtype == 6 && isempty(dev))
        drawnow;
        uigeom = {[1 1] [1 1] [1 1] [1 1]};
        uilist = {{'style' 'text' 'string' 'Sampling frequency:'} ...
                  {'style' 'edit' 'string' fs} ...
                  {'style' 'text' 'string' 'Window type:'} ...
                  {'style' 'popupmenu' 'string' wtypes 'tag' 'wtypepop' 'value' wtype 'callback' {@comwtype, dev}} ...
                  {'style' 'text' 'string' 'Transition bandwidth (Hz):'} ...
                  {'style' 'edit' 'string' df} ...
                  {'style' 'text' 'string' 'Max passband deviation/ripple:' 'tag' 'devtext'} ...
                  {'style' 'edit' 'tag' 'devedit' 'createfcn' {@comwtype, dev}}};
        result = inputgui(uigeom, uilist, 'pophelp(''pop_firwsord'')', 'Estimate filter order -- pop_firwsord()');

        if length(result) == 0, return, end
        if ~isempty(result{1})
            fs = str2num(result{1});
        else
            fs = 2;
        end
        wtype = result{2};
        if ~isempty(result{3})
            df = str2num(result{3});
        else
            error('Not enough input arguments.');
        end
        if ~isempty(result{4})
            dev = str2num(result{4});
        elseif wtype == 6
            error('Not enough input arguments.');
        end
    end

    if length(fs) > 1 || ~isnumeric(fs) || ~isreal(fs) || fs <= 0
        error('Sampling frequency must be a positive real scalar.');
    end
    if length(df) > 1 || ~isnumeric(df) || ~isreal(df) || fs <= 0
        error('Transition bandwidth must be a positive real scalar.');
    end

    df = df / fs; % Normalize transition band width

    if wtype == 6
        if length(dev) > 1 || ~isnumeric(dev) || ~isreal(dev) || dev <= 0
            error('Passband deviation/ripple must be a positive real scalar.');
        end
        devdb = -20 * log10(dev);
        m = 1 + (devdb - 8) / (2.285 * 2 * pi * df);
    else
        dfs = [0.9 2.9 3.1 3.3 5.5];
        m = dfs(wtype) / df;
    end

    m = ceil(m / 2) * 2; % Make filter order even (type 1)

function comwtype(obj, evt, dev)
    enable = {'off' 'off' 'off' 'off' 'off' 'on'};
    devs = {0.089 0.056 0.0063 0.0022 0.0002 dev};
    wtype = get(findobj(gcbf, 'tag', 'wtypepop'), 'value');
    set(findobj(gcbf, 'tag', 'devtext'), 'enable', enable{wtype});
    set(findobj(gcbf, 'tag', 'devedit'), 'enable', enable{wtype}, 'string', devs{wtype});
