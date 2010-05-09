% eegplugin_PatID() - EEGLAB plugin for visually editing event markers and identifying bad channels using eegplot.
%
% Usage:
%   >> eegplugin_PatID(fig, try_strings, catch_stringss);
%
% Inputs:
%   fig            - [integer]  EEGLAB figure
%   try_strings    - [struct] "try" strings for menu callbacks.
%   catch_strings  - [struct] "catch" strings for menu callbacks.
%
% Creates Edit menu option "Visually edit events and identify bad channels"
% and calls pop_VisEd(EEG). 
%
%
% Copyright (C) <2008> <James Desjardins> Brock University
%
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



function eegplugin_VisEd(fig,try_strings,catch_strings)


% find EEGLAB tools menu.
% ---------------------
editmenu=findobj(fig,'label','Edit');

% Create "pop_VisEd" callback cmd.
%---------------------------------------
VisEd_cmd='[EEG,LASTCOM] = pop_VisEd(EEG);';

% add "Visual edit" submenu to the "Edit" menu.
%--------------------------------------------------------------------
uimenu(editmenu, 'label', 'Visually edit events and identify bad channels', 'callback', VisEd_cmd, 'separator', 'on');
