% eegplugin_dipfit2_0() - DIPFIT plugin version 2.0 for EEGLAB menu. 
%                      DIPFIT is the dipole fitting Matlab Toolbox of 
%                      Robert Oostenveld (in collaboration with A. Delorme).
%
% Usage:
%   >> eegplugin_dipfit2_0(fig, trystrs, catchstrs);
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

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.12  2005/03/17 19:03:50  arno
% same
%
% Revision 1.11  2005/03/17 19:03:03  arno
% same
%
% Revision 1.10  2005/03/17 19:01:37  arno
% searching for gridsearch
%
% Revision 1.9  2005/03/17 02:37:05  arno
% comment
%
% Revision 1.8  2005/03/17 02:36:21  arno
% debug storing dataset
%
% Revision 1.7  2005/03/17 02:29:24  arno
% fixing history
%
% Revision 1.6  2005/03/17 02:24:49  arno
% same
%
% Revision 1.5  2005/03/17 02:24:04  arno
% show history
%
% Revision 1.4  2005/03/17 02:08:43  arno
% add path to fieldtrip
%
% Revision 1.3  2005/03/05 02:56:27  arno
% adding fieldtrip folder
%
% Revision 1.2  2005/03/04 23:35:29  arno
% menu text
%
% Revision 1.1  2005/03/04 23:34:36  arno
% Initial revision
%
% Revision 1.24  2003/12/05 01:16:48  arno
% still version
%
% Revision 1.23  2003/12/05 01.08:58  arno
% version
%
% Revision 1.22  2003/12/05 01.04:20  arno
% revision number
%
% Revision 1.21  2003/12/04 01:42:25  arno
% adding path
%
% Revision 1.20  2003/11/26 03:37:14  arno
% add dipfit folder to path
%
% Revision 1.19  2003/11/25 19:20:33  arno
% same
%
% Revision 1.18  2003/11/25 19:18:43  arno
% conforming to smart plugin
%
% Revision 1.17  2003/11/17 20:08:10  arno
% menu
%
% Revision 1.16  2003/11.05 16:21:21  arno
% homogenous -> homogeneous
%
% Revision 1.15  2003/10/31 18:09:33  arno
% sparator position
%
% Revision 1.14  2003/10/31 18:08:32  arno
% menus
%
% Revision 1.13  2003/10/29 16:33:31  arno
% manual fitting
%
% Revision 1.12  2003/10/11.01:53:36  arno
% *** empty log message ***
%
% Revision 1.11  2003/10/11.01:52:51  arno
% *** empty log message ***
%
% Revision 1.10  2003/10/11.01:51:48  arno
% *** empty log message ***
%
% Revision 1.9  2003/10/11.01:51:05  arno
% same
%
% Revision 1.8  2003/10/11.01:50:26  arno
% update
% menu
%
% Revision 1.7  2003/10/11.01:22:05  arno
% adding automatic fitting menu
%
% Revision 1.6  2003/10/09 18:20:22  arno
% update munu for grid
%
% Revision 1.5  2003/06/13 16:53:13  arno
% checking electrodes
%
% Revision 1.4  2003/03/12 00:44:04  arno
% adding menu for plotting components
%
% Revision 1.3  2003/02/26 16:23:03  arno
% menu correction
%
% Revision 1.2  2003/02/24 22:32:50  arno
% updating menus
%
% Revision 1.1  2003/02/24 19:53:26  arno
% Initial revision
%

function vers = eegplugin_dipfit2_0(fig, trystrs, catchstrs)
    
    vers = 'dipfit2.0';
    if nargin < 3
        error('eegplugin_dipfit2_0 requires 3 arguments');
    end;
    
    % add dipfit folder to path
    % -----------------------
    if ~exist('dipolefitting')
        p = which('eegplugin_dipfit2_0');
        p = p(1:findstr(p,'eegplugin_dipfit2_0.m')-1);
        dircontent  = dir([ p '..' p(end) '..' p(end) '..'  ]);
        dircontent  = { dircontent.name };
        ind = strmatch('fieldtrip', lower(dircontent));
        if ~isempty(ind)
            addpath([ p '..' p(end) '..' p(end) '..' p(end) dircontent{ind} ] );
        else
            disp('Warning: Add Fieldtrip folder path manualy or dipfit2 will not be functional');
        end;
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
    check_dipfit = ['if ~isfield(EEG, ''dipfit''), error(''Run the dipole setting first''); end;'  ...
                    'if isempty(EEG.dipfit), error(''Run the dipole setting first''); end;'  ...
                    trystrs.no_check ];
    check_dipfitnocheck = ['if ~isfield(EEG, ''dipfit''), error(''Run the dipole setting first''); end; ' ...
                 trystrs.no_check ];
    check_chans = [ '[EEG tmpres] = eeg_checkset(EEG, ''chanlocs_homogeneous'');' ...
                       'if ~isempty(tmpres), h(tmpres), end; clear tmpres;' ];
    
    % menu callback commands
    % ----------------------
    comauto    = [ trystrs.check_ica check_chans  '[EEG LASTCOM] = pop_multifit(EEG);'        catchstrs.store_and_hist ];
    comsetting = [ trystrs.check_ica check_chans '[EEG LASTCOM]=pop_dipfit_settings(EEG);'    catchstrs.store_and_hist ]; 
    combatch   = [ check_dipfit check_chans  '[EEG LASTCOM] = pop_dipfit_gridsearch(EEG);'    catchstrs.store_and_hist ];
    comfit     = [ check_dipfitnocheck check_chans [ 'EEG = pop_dipfit_nonlinear(EEG); ' ...
                        'LASTCOM = ''% === History not supported for manual dipole fitting ==='';' ]  catchstrs.store_and_hist ];
    % preserve the '=" sign in the comment above: it is used by EEGLAB to detect appropriate LASTCOM
    complot    = [ check_dipfit check_chans 'LASTCOM = pop_dipplot(EEG);'                     catchstrs.add_to_hist ];

    
    % create menus
    % ------------
    submenu = uimenu( menu, 'Label', 'Locate dipoles using DIPFIT 2.x');
    uimenu( submenu, 'Label', 'Autofit components'       , 'CallBack', comauto);
    uimenu( submenu, 'Label', 'Head model and settings'  , 'CallBack', comsetting, 'separator', 'on');
    uimenu( submenu, 'Label', 'Coarse fit (grid scan)'   , 'CallBack', combatch);
    uimenu( submenu, 'Label', 'Fine fit (iterative)'     , 'CallBack', comfit);
    uimenu( submenu, 'Label', 'Plot component dipoles'   , 'CallBack', complot, 'separator', 'on');
