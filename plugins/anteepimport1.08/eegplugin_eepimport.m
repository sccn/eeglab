% eegplugin_eepimport() - EEGLAB plugin for importing ANT EEProbe data files.
%                         With this menu it is possible to import and export a continuous CNT file (*.cnt) or
%                         an averaged file (*.avr), linked with an event file (*.trg).
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
%
% Notes:
%   This plugins consist of the following Matlab files:
%   pop_loadeep.m           pop_loadeep_avg.m
%   loadeep.m               loadeep_avg.m
%   read_eep_cnt.m          read_eep_avr.m
%   read_eep_cnt.mexglx     read_eep_avr.mexglx
%   read_eep_cnt.dll        read_eep_avr.dll
%
% Authors:
% Maarten-Jan Hoeve, ANT, The Netherlands / www.ant-software.nl, 3 October 2003
% Maarten van de Velde, ANT, The Netherlands / www.ant-neuro.com, 21 July 2005
%
% See also: eeglab(), pop_loadeep(), loadeep(), read_eep_cnt(), pop_loadeep_avg(), loadeep_avg(), read_eep_avr()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005 ANT Software, The Netherlands, eeprobe@ant-neuro.com / info@ant-neuro.com
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

% $Log: eegplugin_eepimport.m,v $
% Revision 1.6  2006-09-25 14:25:25  mvelde
% updated version number
%
% Revision 1.5  2006-09-25 14:25:25  mvelde
% updated version number
%
% Revision 1.4  2006-09-25 14:04:02  mvelde
% updated for EEGLAB 5.03
%
% fix reading data for the new EEGLAB version - Arno August 2006
%
% Revision 1.3  2005/07/21 10:11:40  mvelde
% updated menu structure
%
% Revision 1.2  2005/06/08 08:16:37  mvelde
% converted files to unix format
%
% Revision 1.1  2004/11/26 13:22:22  jwiskerke
% Added eeglab to cvs.
%
% Revision 1.3  2003/10/24 13:24:54  Maarten-Jan Hoeve
% Added GNU Licence and updated revision history
%
%
% Advanced Neuro Technology (ANT) BV, The Netherlands, www.ant-neuro.com / info@ant-neuro.com
%

function vers = eegplugin_eepimport(fig, trystrs, catchstrs)

    vers = 'eepimport1.06';

    if nargin < 3
        error('eegplugin_eepimport requires 3 arguments');
    end;

    % add folder to path
    % ------------------
    if ~exist('pop_loadeep')
        p = which('eegplugin_eepimport.m');
        p = p(1:findstr(p,'eegplugin_eepimport.m')-1);
        addpath( p );
    end;

    % find import data menu
    % ---------------------
    menu = findobj(fig, 'tag', 'import data');

    % menu callbacks
    % --------------
    comcnt = [ trystrs.no_check '[EEG LASTCOM] = pop_loadeep;' catchstrs.new_non_empty ];
    comavr = [ trystrs.no_check '[EEG LASTCOM] = pop_loadeep_avg;' catchstrs.new_non_empty ];


    uimenu( menu, 'label', 'From ANT EEProbe .CNT file', 'callback', comcnt, 'separator', 'on' );
    uimenu( menu, 'label', 'From ANT EEProbe .AVR file' , 'callback', comavr );

%   future export
%   uimenu( submenu, 'label', 'To EEProbe CNT file format' , 'callback', expcnt );
%   uimenu( submenu, 'label', 'To EEProbe AVR file format' , 'callback', expavr );

