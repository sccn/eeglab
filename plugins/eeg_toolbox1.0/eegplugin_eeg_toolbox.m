% eegplugin_eeg_toolbox() - EEGLAB plugin for plotting ERPs using the
%             EEG toolbox. WARNING: YOU MUST INSTALL THE EEG TOOLBOX 
%             (ALSO KNOWN NOW AS BIOELECTROMAGNETISM TOOLBOX) YOURSELF 
%             TO MAKE THIS PLUGIN WORK. DOWNLOAD IT FROM 
%
%             http://eeg.sourceforge.net (requires a CVS client)
%
%             OR DOWNLOAD AN OLD COPY (Nov 2006) AT 
%
%             http://ftp/.sccn.ucsd.edu/pub/bioelectromagnetism.zip
%
% Usage:
%   >> eegplugin_eeg_toolbox(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Author: Arnaud Delorme, SCC, INC, UCSD

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


function vers = eegplugin_eeg_toolbox(fig, trystrs, catchstrs)

    vers = 'EEG toolbox ERP plotting';
    if nargin < 3
        error('eegplugin_eeg_toolbox requires 3 arguments');
    end;
    
    % add folder to path
    % ------------------
    if ~exist('eegplugin_eeg_toolbox')
        p = which('eegplugin_eeg_toolbox.m');
        p = p(1:findstr(p,'eegplugin_eeg_toolbox.m')-1);
        addpath( p );
    end;
    
    if ~exist('gui_erp_plot')
        vers = 'EEG toolbox not installed, plugin disabled';
        return;
    end;
        
    % find import data menu
    % ---------------------
    menui = findobj(fig, 'tag', 'tools');
    
    % menu callbacks
    % --------------
    comcnt = [ 'gui_erp_plot_eeglab(EEG);' ];
                
    % create menus
    % ------------
    uimenu( menui, 'label', 'Peak detection using EEG toolbox', 'callback', comcnt);
