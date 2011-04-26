% eegplugin_dipfit() - DIPFIT plugin version 2.0 for EEGLAB menu. 
%                      DIPFIT is the dipole fitting Matlab Toolbox of 
%                      Robert Oostenveld (in collaboration with A. Delorme).
%
% Usage:
%   >> eegplugin_dipfit(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Notes:
%   To create a new plugin, simply create a file beginning with "eegplugin_"
%   and place it in your eeglab folder. It will then be automatically 
%   detected by eeglab. See also this source code internal comments.
%   For eeglab to return errors and add the function's results to 
%   the eeglab history, menu callback must be nested into "try" and 
%   a "catch" strings. For more information on how to create eeglab 
%   plugins, see http://www.sccn.ucsd.edu/eeglab/contrib.html
%
% Author: Arnaud Delorme, CNL / Salk Institute, 22 February 2003
%
% See also: eeglab()

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

function vers = eegplugin_dipfit(fig, trystrs, catchstrs)
    
    vers = 'dipfit2.2';
    if nargin < 3
        error('eegplugin_dipfit requires 3 arguments');
    end;
    
    % add dipfit folder to path
    % -----------------------
    if ~exist('ft_dipolefitting')
        p = which('eegplugin_dipfit');
        disp('Fieldtrip functions for dipole localization not found, removing Dipfit2');
        p = p(1:findstr(p,'eegplugin_dipfit.m')-1);
        rmpath(p);
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
    check_dipfit = [trystrs.no_check 'if ~isfield(EEG, ''dipfit''), error(''Run the dipole setting first''); end;'  ...
                    'if isempty(EEG.dipfit), error(''Run the dipole setting first''); end;'  ];
    check_dipfitnocheck = [ trystrs.no_check 'if ~isfield(EEG, ''dipfit''), error(''Run the dipole setting first''); end; ' ];
    check_chans = [ '[EEG tmpres] = eeg_checkset(EEG, ''chanlocs_homogeneous'');' ...
                       'if ~isempty(tmpres), eegh(tmpres), end; clear tmpres;' ];
    
    % menu callback commands
    % ----------------------
    comsetting = [ trystrs.check_ica check_chans '[EEG LASTCOM]=pop_dipfit_settings(EEG);'    catchstrs.store_and_hist ]; 
    combatch   = [ check_dipfit check_chans  '[EEG LASTCOM] = pop_dipfit_gridsearch(EEG);'    catchstrs.store_and_hist ];
    comfit     = [ check_dipfitnocheck check_chans [ 'EEG = pop_dipfit_nonlinear(EEG); ' ...
                        'LASTCOM = ''% === History not supported for manual dipole fitting ==='';' ]  catchstrs.store_and_hist ];
    comauto    = [ check_dipfit check_chans  '[EEG LASTCOM] = pop_multifit(EEG);'        catchstrs.store_and_hist ];
    % preserve the '=" sign in the comment above: it is used by EEGLAB to detect appropriate LASTCOM
    complot    = [ check_dipfit check_chans 'LASTCOM = pop_dipplot(EEG);'                     catchstrs.add_to_hist ];

    
    % create menus
    % ------------
    submenu = uimenu( menu, 'Label', 'Locate dipoles using DIPFIT 2.x', 'separator', 'on');
    uimenu( submenu, 'Label', 'Head model and settings'  , 'CallBack', comsetting);
    uimenu( submenu, 'Label', 'Coarse fit (grid scan)'   , 'CallBack', combatch);
    uimenu( submenu, 'Label', 'Fine fit (iterative)'     , 'CallBack', comfit);
    uimenu( submenu, 'Label', 'Autofit (coarse fit, fine fit & plot)', 'CallBack', comauto);
    uimenu( submenu, 'Label', 'Plot component dipoles'   , 'CallBack', complot, 'separator', 'on');
