% eegplugin_ascinstep() - EEGLAB plugin for importing ASC INStep files.
%
% Usage:
%   >> eegplugin_ascinstep(menu);
%   >> eegplugin_ascinstep(menu, trystrs, catchstrs);
%
% Inputs:
%   menu       - [float]  EEGLAB menu handle
%   trystrs    - [struct] "try" strings for menu callbacks. See notes on EEGLab plugins.
%                (http://www.sccn.ucsd.edu/eeglab/contrib.html)
%   catchstrs  - [struct] "catch" strings for menu callbacks. See notes on EEGLab plugins. 
%                (http://www.sccn.ucsd.edu/eeglab/contrib.html)   
%
% 
% Notes:
%   This plugins consist of the following Matlab files:
%   pop_loadeep.m           pop_loadeep_avg.m
%   loadeep.m               loadeep_avg.m
%   read_eep_cnt.m          read_eep_avr.m
%   read_eep_cnt.mexglx     read_eep_avr.mexglx
%   read_eep_cnt.dll        read_eep_avr.dll
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, 2006

% Copyright (C) 2006 Arnaud Delorme
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

function eegplugin_ascinstep(fig, trystrs, catchstrs)

    if nargin < 3
        error('eegplugin_ascinstep requires 3 arguments');
    end;
  
    % add folder to path
    % ------------------
    if ~exist('pop_loadeep')
        p = which('eegplugin_ascinstep.m');
        p = p(1:findstr(p,'eegplugin_ascinstep.m')-1);
        addpath( p );
    end;
    
    % find import data menu
    % ---------------------
    menu = findobj(fig, 'tag', 'import data');

    % menu callbacks
    % --------------
    comcnt = [ trystrs.no_check '[EEG LASTCOM] = pop_loadascinstep;' catchstrs.new_and_hist ];
      
    uimenu( menu, 'label', 'From INStep .ASC file', 'callback', comcnt, 'separator', 'on' );

