% eegplugin_bva_io() - EEGLAB plugin for importing Brainvision 
%                         .vhdr data files.
%
% Usage:
%   >> eegplugin_bva_io(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Author: Andreas Widmann for binary import, 2004
%         Arnaud Delorme for Matlab import and EEGLAB interface
%
% See also: pop_loadbv()

% Copyright (C) 2004 Andreas Widmann & Arnaud Delorme
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

% $Id: eegplugin_bva_io.m 53 2010-05-22 21:57:38Z arnodelorme $

function vers = eegplugin_bva_io(fig, trystrs, catchstrs)

    vers = 'bva_io1.58';
    if nargin < 3
        error('eegplugin_bva_io requires 3 arguments');
    end;
    
    % add folder to path
    % ------------------
    if ~exist('eegplugin_bva_io')
        p = which('eegplugin_bva_io.m');
        p = p(1:findstr(p,'eegplugin_bva_io.m')-1);
        addpath( p );
    end;
    
    % find import data menu
    % ---------------------
    menui = findobj(fig, 'tag', 'import data');
    menuo = findobj(fig, 'tag', 'export');
    
    % menu callbacks
    % --------------
    icadefs;
    versiontype = 1;
    if exist('EEGLAB_VERSION')
        if EEGLAB_VERSION(1) == '4'
            versiontype = 0;
        end;
    end;
    if versiontype == 0
        comcnt1 = [ trystrs.no_check '[EEGTMP LASTCOM] = pop_loadbv;'  catchstrs.new_non_empty ];
        comcnt2 = [ trystrs.no_check '[EEGTMP LASTCOM] = pop_loadbva;' catchstrs.new_non_empty ];
    else
        comcnt1 = [ trystrs.no_check '[EEG LASTCOM] = pop_loadbv;'  catchstrs.new_non_empty ];
        comcnt2 = [ trystrs.no_check '[EEG LASTCOM] = pop_loadbva;' catchstrs.new_non_empty ];
    end;
    comcnt3 = [ trystrs.no_check 'LASTCOM = pop_writebva(EEG);'  catchstrs.add_to_hist ];
                
    % create menus
    % ------------
    uimenu( menui, 'label', 'From Brain Vis. Rec. .vhdr file',  'callback', comcnt1, 'separator', 'on' );
    uimenu( menui, 'label', 'From Brain Vis. Anal. Matlab file', 'callback', comcnt2 );
    uimenu( menuo, 'label', 'Write Brain Vis. exchange format file',  'callback', comcnt3, 'separator', 'on' );
