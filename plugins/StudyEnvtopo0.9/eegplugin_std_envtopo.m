% eegplugin_std_envtopo() - std_envtopo plugin for updated std_envtopo
%
% Usage:
%   >> eegplugin_std_envtopo(fig, try_strings, catch_strings);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Notes:
%   To create a new plugin, simply create a file beginning with "eegplugin_"
%   and place it in your eeglab folder. It will then be automatically
%   detected by eeglab. See also this source code internal comments.
%   For eeglab to return errors and add the function's results to
%   the eeglab history, menu callback must be nested into "try" and
%   a "catch" strings. For more information on how to create eeglab
%   plugins, see http://www.sccn.ucsd.edu/eeglab/contrib.html
%
% Author: Makoto Miyakoshi & Arnaud Delorme, JSPS/SCCN, INC, UCSD
%
% See also: pop_std_envtopo, std_envtopo, envtopo, pop_envtopo

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2011, Makoto Miyakoshi
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
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA
%
% History
%
% 01/12/2012 ver 1.2 by Makoto. Button renamed with caution.
% 12/23/2011 ver 1.1 by Makoto. Renamed the function.
% 08/01/2011 ver 1.0 by Makoto. Created.

function vers = eegplugin_std_evntopo( fig, trystrs, catchstrs);


vers = 'STUDY envtopo 0.9';
if nargin < 3
    error('eegplugin_std_evntopo requires 3 arguments');
end;

% add folder to path
% ------------------
if exist('std_envtopo', 'file')
    p = which('eegplugin_std_evntopo.m');
    p = p(1:findstr(p,'eegplugin_std_evntopo.m')-1);
    addpath(p);
end;

% add menu to clustedit
% ---------------------
menu = findobj(fig, 'Label', 'Edit/plot clusters');
structure.uilist = { { } ...
    {'style' 'pushbutton' 'string' 'Plot Envtopo plugin (beta)' 'Callback' 'args = { STUDY ALLEEG}; [STUDY LASTCOM] = pop_std_envtopo(args{:}); clear tmpargs;' } { } { } };
structure.geometry = { [1] [1 0.3 1] };
arg = vararg2str( { structure } );
cb_clustedit   = [ trystrs.no_check 'ALLEEGTMP = ALLEEG; [STUDYTMP LASTCOM] = pop_clustedit(STUDY, ALLEEG, [], ' arg ');' catchstrs.update_study];
set(menu, 'callback', cb_clustedit);
