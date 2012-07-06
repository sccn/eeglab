% eegplugin_ERPLAB_filters() - EEGLAB plugin for filtering data using linear-
%                       phase FIR filters
%
% Usage:
%   >> eegplugin_ERPLAB_filters(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Author: ERPLAB Team, Universidad de California, Davis, 2007


function vers = eegplugin_ERPLAB_filters(fig, trystrs, catchstrs)

    vers = 'Butter1.0';
    if nargin < 3
        error('eegplugin_firfilt requires 3 arguments');
    end

    % add folder to path
    % -----------------------
    if ~exist('pop_ERPLAB_butter1','file')   % Este es el archivo que maneja la interface de usuario
        p = which('eegplugin_ERPLAB_filters');
        p = p(1:findstr(p,'eegplugin_ERPLAB_filters.m')-1);
        addpath([p vers]);
    end

    % find import data menu
    % ---------------------
    menu = findobj(fig, 'tag', 'filter');

    % menu callbacks
    % --------------
    combutter1 = [trystrs.no_check '[EEG LASTCOM] = pop_ERPLAB_butter1(EEG);' catchstrs.new_and_hist];
    compolydet = [trystrs.no_check '[EEG LASTCOM] = pop_ERPLAB_polydetrend(EEG);' catchstrs.new_and_hist];
%     comfirpm = [trystrs.no_check '[EEG LASTCOM] = pop_firpm(EEG);' catchstrs.new_and_hist];
%     comfirma = [trystrs.no_check '[EEG LASTCOM] = pop_firma(EEG);' catchstrs.new_and_hist];

    % create menus if necessary
    % -------------------------
    uimenu( menu, 'Label', 'ERPLAB Butterworth Filter', 'CallBack', combutter1, 'Separator', 'on');
    uimenu( menu, 'Label', 'ERPLAB Polynomial Detrending', 'CallBack', compolydet);
    %     uimenu( menu, 'Label', 'Parks-McClellan (equiripple) FIR filter', 'CallBack', comfirpm);
%     uimenu( menu, 'Label', 'Moving average FIR filter', 'CallBack', comfirma);
