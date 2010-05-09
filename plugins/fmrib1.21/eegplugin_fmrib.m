% eegplugin_fmrib() -       EEGLAB plugin from the Centre for 
%                           Functional MRI of the Brain (FMRIB), University of Oxford
%                           with tools to remove artifacts associated with
%                           simultaneous EEG/FMRI acquisition from EEG data.
%
% Usage:
%   >> eegplugin_fmrib(fig, try_strings, catch_strings);
%   
%  This is the main plugin program for the FMRIB EEGLAB plugin.
%  To use this plug in, place the plugin folder (fmribx.x) in the EEGLAB 
%  directory inside the 'plugins' folder.  EEGLAB should detect the plugin
%  on start up.  The FMRIB plugin commands will appear under the 'Tools' menu
%  as 'FMRIB tools'.
%  There will be three options:
%      1.  FASTR: Remove fmri gradient artifacts
%      2.  Detect QRS events
%      3.  Remove pluse artifacts
%
%  See the help for each function for more info on usage.
%
% Inputs:
%   fig            - [integer]  handle to EEGLAB figure
%   try_strings    - [struct] "try" strings for menu callbacks.
%   catch_strings  - [struct] "catch" strings for menu callbacks. 
%
% Notes:
%   See Contents.m for the contents of this plugin.
% 
% Create a plugin:
%   For more information on how to create an EEGLAB plugin see the
%   help message of eegplugin_besa() or visit http://www.sccn.ucsd.edu/eeglab/contrib.html
%
%   Author:  Rami K. Niazy
%   
%   Copyright (c) 2006 University of Oxford

% Copyright (C) 2006 University of Oxford
% Author:   Rami K. Niazy, FMRIB Centre
%           rami@fmrib.ox.ac.uk
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

% APR 26, 2006
% Version update

% JAN 13, 2006
% Version update

% SEP 16, 2005
% Version update

% AUG 16, 2005
% Version update

% DEC 23, 2004
% (c) & Version update


function vers = eegplugin_fmrib(fig,try_strings,catch_strings)

vers='fmrib1.21';

if nargin < 3
    error('eegplugin_fmrib requires 3 arguments');
end;
    
% add folder to path
% ------------------
if ~exist('pop_fmrib_fastr.m','file')
    p = which('eegplugin_fmrib.m');
    p = p(1:findstr(p,'eegplugin_fmrib.m')-1);
    addpath([ p vers ] );
end;

% find tools menu
%------------------
toolsmenu=findobj(fig,'tag','tools');

% construct command
%---------------------
fastrcmd = [ try_strings.no_check '[EEG LASTCOM] = pop_fmrib_fastr(EEG);' catch_strings.new_and_hist ];
pascmd=[ try_strings.no_check '[EEG LASTCOM] = pop_fmrib_pas(EEG);' catch_strings.new_and_hist ];
qrsdetectcmd=[ try_strings.no_check '[EEG LASTCOM] = pop_fmrib_qrsdetect(EEG);' catch_strings.new_and_hist ];

% add menu
%----------
fmribmenu=uimenu(toolsmenu,'label','FMRIB Tools','separator','on','tag','fmrib tools');
fastrmenu=uimenu(fmribmenu,'label','FASTR: Remove FMRI gradient artifacts','tag','fastr menu','callback',fastrcmd);
qrsdetectmenu=uimenu(fmribmenu,'label','Detect QRS events','separator','on','tag','qrsdetect menu','callback',qrsdetectcmd);
pasmenu=uimenu(fmribmenu,'label','Remove pulse artifacts','tag','pas menu','callback',pascmd);
