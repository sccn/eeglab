% eegplugin_4dneuroimaging() - EEGLAB plugin for importing 4D neuroimaging PDF files.
%
% Usage:
%   >> eegplugin_4dneuroimaging(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks. 
%
% Create a plugin:
%   For more information on how to create an EEGLAB plugin see the
%   help message of eegplugin_besa() or visit http://www.sccn.ucsd.edu/eeglab/contrib.html
%
% Author: Christian Wienbruch, University of Konstanz, Clinical Psychology, 2005
%         (plugin interface: Arnaud Delorme, SCCN, INC, UCSD)

% Copyright (C) 2005 Christian Wienbruch, University of Konstanz, Clinical Psychology
%

function vers = eegplugin_4dneuroimaging(fig, trystrs, catchstrs)
    
  vers = '4dneuroimaging1.00';
  if nargin < 3
      error('eegplugin_4dneuroimaging requires 3 arguments');
  end;
    
  % add folder to path
  % ------------------
  if ~exist('pop_read4d')
      p = which('eegplugin_4dneuroimaging.m');
      p = p(1:findstr(p,'eegplugin_4dneuroimaging.m')-1);
      addpath([ p vers ] );
  end;

  % find import data menu
  % ---------------------
  menu = findobj(fig, 'tag', 'import data');
  
  % menu callbacks
  % --------------
  com = [ trystrs.no_check '[EEGTMP LASTCOM]= pop_read4d;' catchstrs.new_non_empty ];
  
  % create menus
  % ------------
  uimenu( menu, 'label',  'From 4D .m4d pdf file', 'callback', com, 'separator', 'on' );
