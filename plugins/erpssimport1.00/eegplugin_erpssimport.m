% eegplugin_erpssimport() - EEGLAB plugin for importing ERPSS data files.
%
% Usage:
%   >> eegplugin_eepimport(menu);
%   >> eegplugin_eepimport(menu, trystrs, catchstrs);
%
% Inputs:
%   menu       - [float]  EEGLAB menu handle
%   trystrs    - [struct] "try" strings for menu callbacks. See notes on EEGLab plugins.
%                (http://www.sccn.ucsd.edu/eeglab/contrib.html)
%   catchstrs  - [struct] "catch" strings for menu callbacks. See notes on EEGLab plugins. 
%                (http://www.sccn.ucsd.edu/eeglab/contrib.html)   
%
% Authors: see associated functions

% Copyright (C) 2009 Arnaud Delorme, arno@salk.edu
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

function vers = eegplugin_erpssimport(fig, trystrs, catchstrs)

    vers = 'erpssimport1.00';
    
    if nargin < 3
        error('eegplugin_erpssimport requires 3 arguments');
    end;
  
    % add folder to path
    % ------------------
    if ~exist('pop_loadeep')
        p = which('eegplugin_erpssimport.m');
        p = p(1:findstr(p,'eegplugin_erpssimport.m')-1);
        addpath( p );
    end;
    
    % find import data menu
    % ---------------------
    menu = findobj(fig, 'tag', 'import data');

    % menu callbacks
    % --------------
 	cb_read_erpss  = [ trystrs.no_check '[EEG LASTCOM] = pop_read_erpss;' catchstrs.new_non_empty ]; 
  
    uimenu( menu, 'label','From ERPSS .RAW or .RDF file', 'CallBack', cb_read_erpss, 'Separator', 'on'); 

%   future export
%   uimenu( submenu, 'label', 'To EEProbe CNT file format' , 'callback', expcnt );
%   uimenu( submenu, 'label', 'To EEProbe AVR file format' , 'callback', expavr );

