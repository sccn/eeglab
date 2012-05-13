% eegplugin_procom() - EEGLAB plugin for importing Procom Infinity data files.
%
% Usage:
%   >> eegplugin_procom(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.


% Copyright (C) 2011 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function vers = eegplugin_procom(fig, trystrs, catchstrs)

    vers = 'procom';
    if nargin < 3
        error('eegplugin_procom requires 3 arguments');
    end;
    
    % add folder to path
    % ------------------
    if exist('pop_importpi', 'file')
        p = which('eegplugin_procom.m');
        p = p(1:findstr(p,'eegplugin_procom.m')-1);
        addpath(p);
    end;
    
    % find import data menu
    % ---------------------
    menu = findobj(fig, 'tag', 'import data');
    
    % menu callbacks
    % --------------
    comcnt = [ trystrs.no_check '[EEG LASTCOM] = pop_importpi;' catchstrs.new_and_hist ];
    
    % create menus
    % ------------
    uimenu( menu, 'label', 'From Procom Infinity Text File', 'callback', comcnt, 'separator', 'on');
