% eegplugin_stderpimage() - EEGLAB plugin for importing stderpimage MATLAB data files.
%
% Usage:
%   >> eegplugin_stderpimage(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Notes:
%   This plugins consist of the following Matlab files:
%
% Create a plugin:
%   For more information on how to create an EEGLAB plugin see the
%   help message of eegplugin_besa() or visit http://www.sccn.ucsd.edu/eeglab/contrib.html
%
% Author: Arnaud Delorme (SCCN, UCSD) and Daren Weber (DNL, UCSF)
%         for the EEGLAB interface. Fredrick Carver (NIH) and Daren Weber for
%         the CTF Matlab reading functions.
%
% See also: pop_ctf_read(), ctf_read(), ctf_readmarkerfile()

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

function vers = eegplugin_stderpimage(fig, trystrs, catchstrs)

    vers = 'stderpimage';
    if nargin < 3
        error('eegplugin_stderpimage requires 3 arguments');
    end;
    
    % add folder to path
    % ------------------
    if exist('readbdf', 'file')
        p = which('eegplugin_stderpimage.m');
        p = p(1:findstr(p,'eegplugin_stderpimage.m')-1);
        addpath(p);
    end;

    % add menu to clustedit
    % ---------------------
    menu = findobj(gcf, 'Label', 'Edit/plot clusters');
    structure.uilist = { { } ...
        {'style' 'pushbutton' 'string' 'Plot ERPimage'    'Callback' 'tmpargs = { ''onecomp''   gcf }; stderpimageplugin_plot(tmpargs{:}); clear tmpargs;' } ...
        {'style' 'pushbutton' 'string' 'Params'           'Callback' 'tmpargs = { ''params''    gcf }; stderpimageplugin_plot(tmpargs{:}); clear tmpargs;' } ... 
        {'style' 'pushbutton' 'string' 'Plot ERPimage(s)' 'Callback' 'tmpargs = { ''oneclust''  gcf }; stderpimageplugin_plot(tmpargs{:}); clear tmpargs;' } };
    structure.geometry = { [1] [1 0.3 1] };
    arg = vararg2str( { structure } );
    cb_clustedit   = [ trystrs.no_check 'ALLEEGTMP = ALLEEG; [STUDYTMP LASTCOM] = pop_clustedit(STUDY, ALLEEG, [], ' arg ');' catchstrs.update_study];
    set(menu, 'callback', cb_clustedit);
    
    % add menu to pop_chanplot
    % ------------------------
    menu = findobj(gcf, 'Label', 'Plot channel measures');
    structure.uilist = { { } ...
        {'style' 'pushbutton' 'string' 'Plot ERPimage'    'Callback' 'tmpargs = { ''changrp''   gcf }; stderpimageplugin_plot(tmpargs{:}); clear tmpargs;' } ...
        {'style' 'pushbutton' 'string' 'Params'           'Callback' 'tmpargs = { ''params''    gcf }; stderpimageplugin_plot(tmpargs{:}); clear tmpargs;' } ... 
        {'style' 'pushbutton' 'string' 'Plot ERPimage(s)' 'Callback' 'tmpargs = { ''onechan''   gcf }; stderpimageplugin_plot(tmpargs{:}); clear tmpargs;' } };
    structure.geometry = { [1] [1 0.3 1] };
    arg = vararg2str( { structure } );
    cb_clustedit   = [ trystrs.no_check 'ALLEEGTMP = ALLEEG; [STUDYTMP LASTCOM] = pop_chanplot(STUDY, ALLEEG, ' arg ');' catchstrs.update_study];
    set(menu, 'callback', cb_clustedit);
    
    
    