% eegplugin_loreta() - EEGLAB plugin for exporting/importing component
%                      scalp maps to LORETA.
% Usage:
%   >> eegplugin_loreta(fig);
%   >> eegplugin_loreta(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] EEGLAB figure
%   trystrs    - [struct]  "try" strings for menu callbacks. See notes. 
%   catchstrs  - [struct]  "catch" strings for menu callbacks. See notes. 
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2005
%
% See also: eeglab()

% Copyright (C) 2005 Arnaud Delorme, SCCN, INC, UCSD, 2005 arno@salk.edu
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

function vers = eegplugin_loreta(fig, trystrs, catchstrs)
    
    vers = 'loreta1.0';
    if nargin < 3
        error('eegplugin_loreta requires 3 arguments');
    end;
    
    % add loreta folder to path
    % -----------------------
    if ~exist('loreta_importcomp')
        p = which('eegplugin_loreta');
        p = p(1:findstr(p,'eegplugin_loreta.m')-1);
        addpath([ p vers ] );
    end;
    
    % find tools menu
    % ---------------
    menu = findobj(fig, 'tag', 'tools'); 
    % tag can be 
    % 'import data'  -> File > import data menu
    % 'import epoch' -> File > import epoch menu
    % 'import event' -> File > import event menu
    % 'export'       -> File > export
    % 'tools'        -> tools menu
    % 'plot'         -> plot menu
    
    % command to check that the '.source' is present in the EEG structure 
    % -------------------------------------------------------------------
    check_loreta = trystrs.no_check;
    
    % menu callback commands
    % ----------------------
    catchstrs.store_and_hist = [ ...
                        '[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); ' ...
                        'catch, errordlg2(lasterr, ''EEGLAB error''); LASTCOM= ''''; clear EEGTMP; end;' ...
                        'h(LASTCOM); disp(''Done.''); eeglab(''redraw'');' ];
    comexport = [ trystrs.check_ica 'pop_eeglab2loreta(EEG);' catchstrs.add_to_hist ]; 
    %comimport = [ trystrs.check_ica 'EEG = loretaimport(EEG);' catchstrs.store_and_hist ];
    %complot1  = [ check_loreta 'pop_dipplot(EEG);' catchstrs.add_to_hist ];
    
    % create menus
    % ------------
    submenu = uimenu( menu, 'Label', 'Locate dipoles using LORETA');
    uimenu( submenu, 'Label', 'Export components to LORETA'     , 'CallBack', comexport);
    %uimenu( submenu, 'Label', 'Import dipoles from LORETA'      , 'CallBack', '', 'enable', 'off');
    %uimenu( submenu, 'Label', 'Plot dipoles on LORETA head'     , 'CallBack', '', 'enable', 'off');
 
