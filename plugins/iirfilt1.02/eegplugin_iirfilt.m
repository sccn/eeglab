% eegplugin_iirfilt() - EEGLAB plugin for importing data using IIRFILT Matlab toolbox
%
% Usage:
%   >> eegplugin_iirfilt(fig, trystrs, catchstrs);
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
% Authors: Maksym Pozdin (mpozdin.ece04@gtalumni.org, IOL/ONRC,2004), 
%          with Arnaud Delorme and Scott Makeig (SCCN/INC/UCSD, La Jolla CA)

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
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function vers = eegplugin_iirfilt(fig, trystrs, catchstrs)
    
    vers = 'iirfilt1.01';
  if nargin < 3
      error('eegplugin_iirfilt requires 3 arguments');
  end;
    
  % add besa folder to path
  % -----------------------
  if ~exist('pop_iirfilt')
      p = which('eegplugin_iirfilt');
      p = p(1:findstr(p,'eegplugin_iirfilt.m')-1);
      addpath([ p vers ] );
  end;
  
  % find import data menu
  % ---------------------
  menu = findobj(fig, 'tag', 'filter');
  
  % menu callbacks
  % --------------
  combio = [ trystrs.no_check '[EEG LASTCOM] = pop_iirfilt(EEG);' catchstrs.new_and_hist ]; 
  
  % create menus if necessary
  % -------------------------
  uimenu( menu, 'Label', 'Short non-linear IIR filter',  'CallBack', combio, 'Separator', 'on'); 
