% plugin_installstartup() - install popular toolboxes when EEGLAB starts
%
% Usage:
%   >> restartEeglabFlag = plugin_installstartup; % pop up window
%
% Outputs:
%   restartEeglabFlag - [0|1] restart EEGLAB
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, Oct. 29, 2013-

% Copyright (C) 2013 Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
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

function restartEeglabFlag = plugin_installstartup 

%                 { 'style' 'checkbox' 'string' 'Neuroelectromagnetic Forward Head Modeling Toolbox (NFT): Advanced source localization tools - beta (Zeynep Akalin Acar, 100 Mb)' 'value' 1 'enable' 'on' } ...

    uilist = { { 'style' 'text' 'String' 'Download and install some popular third party EEGLAB plug-ins' 'fontweight', 'bold','tag', 'title'} ...
                 { } ...
                 { 'style' 'checkbox' 'string' 'Brain Vision Analyser data import plugin (40 Kb)' 'value' 1 'enable' 'on' } ...
                 { 'style' 'checkbox' 'string' 'ANT data import plugin (800 Kb)' 'value' 1 'enable' 'on' } ...
                 { 'style' 'checkbox' 'string' 'Measure Projection Toolbox (MPT) for ulti-subject ICA analysis (415 Mb)' 'value' 1 'enable' 'on' } ...
                 { 'style' 'checkbox' 'string' 'Source Information Flow Toolbox (SIFT) for causal analysis of EEG sources (600 Mb)' 'value' 1 'enable' 'on' } ...
                 { 'style' 'checkbox' 'string' 'BCILAB for real-time and offline BCI platform (200 Mb)' 'value' 1 'enable' 'on' } ...
                 { 'style' 'checkbox' 'string' 'LIMO for linear analysis of MEEG data using single trials (2.5 Mb)' 'value' 1 'enable' 'on' } ...
                 { } ...
                 { 'style' 'radiobutton' 'string' 'Do not show this query again' 'value' 0 } ...
                 {} ...
                 { 'style' 'text' 'String' 'Note: manage plug-in tools using EEGLAB menu item, "File > Plug-ins > Manage plug-ins"' } ...
                 {} ...
                 { 'Style', 'pushbutton', 'string', 'Do not install now', 'tag' 'cancel' 'callback', 'set(gcbf, ''userdata'', ''cancel'');' } { } ...
                 { 'Style', 'pushbutton', 'string', 'Install plugins now', 'tag', 'ok', 'callback', 'set(gcbf, ''userdata'', ''ok'');' } ...
                 };
             
    geomline = [1];        
    geom     = { [1] [1] geomline geomline geomline geomline geomline geomline [1] geomline [1] [1] [1] [1 1 1]};
    geomvert = [ 1   0.5 1        1        1        1        1        1        0.3 1        0.3 1   1   1]; 

    %result = inputgui( 'geometry', geom, 'uilist', uilist, 'helpcom', 'pophelp(''plugin_installstartup'')', 'title', 'Install popular plugins', 'geomvert', geomvert, 'eval', evalstr);
    fig = figure('visible', 'off');
    [tmp1 tmp2 handles] = supergui( 'geomhoriz', geom, 'uilist', uilist, 'title', 'Install popular plugins', 'geomvert', geomvert, 'fig', fig); %, 'eval', evalstr);
    set(findobj(fig, 'tag', 'title'), 'fontsize', 16);
    waitfor( fig, 'userdata');
    
    % decode inputs
    handles
    results = cellfun(@(x)(get 
    get(handles, 'value')
    