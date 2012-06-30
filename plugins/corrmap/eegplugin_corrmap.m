% eegplugin_corrmap() - EEGLAB plugin to cluster ICs based on a correlation approach.
%
% Usage:
%   >> eegplugin_corrmap(fig, try_strings, catch_strings);
%
% Inputs:
%   fig           - [integer]  EEGLAB figure
%   try_strings   - [struct] "try" strings for menu callbacks.
%   catch_strings - [struct] "catch" strings for menu callbacks.
%
% See also:  pop_corrmap(), corrmap()
%
% Author: F. Campos Viola, MRC-IHR, 19/10/2007

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) F. Campos Viola, MRC-IHR
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

% revised by F Campos-Viola - corrmap1.01 (30/01/2009)

function vers=eegplugin_corrmap( fig, try_strings, catch_strings);

 vers = 'corrmap1.02';
    if nargin < 3
        error('eegplugin_corrmap requires 3 arguments');
    end
    
std = findobj(fig, 'tag', 'study');
uimenu(std, 'label', 'Cluster components by correlation (CORRMAP)', 'callback', ...
    [try_strings.check_ica '[CORRMAP STUDY ALLEEG LASTCOM]= pop_corrmap(STUDY,ALLEEG);' catch_strings.add_to_hist ], 'userdata', 'startup:off;study:on');