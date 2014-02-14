% eegplugin_firfilt() - EEGLAB plugin for filtering data using linear-
%                       phase FIR filters
%
% Usage:
%   >> eegplugin_firfilt(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Author: Andreas Widmann, University of Leipzig, Germany, 2005

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

function vers = eegplugin_firfilt(fig, trystrs, catchstrs)

    vers = 'firfilt1.6.1';
    if nargin < 3
        error('eegplugin_firfilt requires 3 arguments');
    end

    % add folder to path
    % -----------------------
    if ~exist('pop_firws')
        p = which('eegplugin_firfilt');
        p = p(1:findstr(p,'eegplugin_firfilt.m')-1);
        addpath([p vers]);
    end

    % find import data menu
    % ---------------------
    menu = findobj(fig, 'tag', 'filter');

    % menu callbacks
    % --------------
    comfirfiltnew = [trystrs.no_check '[EEG LASTCOM] = pop_eegfiltnew(EEG);' catchstrs.new_and_hist];
    comfirws = [trystrs.no_check '[EEG LASTCOM] = pop_firws(EEG);' catchstrs.new_and_hist];
    comfirpm = [trystrs.no_check '[EEG LASTCOM] = pop_firpm(EEG);' catchstrs.new_and_hist];
    comfirma = [trystrs.no_check '[EEG LASTCOM] = pop_firma(EEG);' catchstrs.new_and_hist];

    % create menus if necessary
    % -------------------------
    uimenu( menu, 'Label', 'Basic FIR filter (new, default)', 'CallBack', comfirfiltnew, 'Separator', 'on', 'position', 1);
    uimenu( menu, 'Label', 'Windowed sinc FIR filter', 'CallBack', comfirws, 'position', 2);
    uimenu( menu, 'Label', 'Parks-McClellan (equiripple) FIR filter', 'CallBack', comfirpm, 'position', 3);
    uimenu( menu, 'Label', 'Moving average FIR filter', 'CallBack', comfirma, 'position', 4);
  