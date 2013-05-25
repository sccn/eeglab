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


    uilist = { { 'style' 'text' 'String' 'You may now automatically download and install some recommended EEGLAB plug-ins:' 'fontweight', 'bold','tag', 'title'} ...
                 { } ...
                 { 'style' 'checkbox' 'string' 'Measure Projection Toolbox (MPT): Multi-subject analysis and visualization of ICA component differences (Nima Bigdely-Shamlo, 415 Mb)' 'value' 1 'enable' 'on' } ...
                 { 'style' 'checkbox' 'string' 'Source Information Flow Toolbox (SIFT): Causal analysis of EEG sources using Granger causality and other methods (Tim Mullen, 600 Mb)' 'value' 1 'enable' 'on' } ...
                 { 'style' 'checkbox' 'string' 'BCILAB: real-time and offline BCI platform including cross-validation tools (Christian Kothe, 200 Mb)' 'value' 1 'enable' 'on' } ...
                 { 'style' 'checkbox' 'string' 'Neuroelectromagnetic Forward Head Modeling Toolbox (NFT): Advanced source localization tools - beta (Zeynep Akalin Acar, 100 Mb)' 'value' 1 'enable' 'on' } ...
                 { 'style' 'checkbox' 'string' 'LIMO: Linear processing of MEEG data using single trials and hierarchical linear models (Cyril Pernet, 2.5 Mb)' 'value' 1 'enable' 'on' } ...
                 { } ...
                 { 'style' 'checkbox' 'string' 'Do not show this query again' 'value' 0 } ...
                 {} ...
                 { 'style' 'text' 'String' 'Note: You may also download or update these and other plug-in tools using the EEGLAB menu item, "File > Plug-ins > Manage plug-ins".' } };
             
    geomline = [1];        
    geom     = { [1] [1] geomline geomline geomline geomline geomline [1] geomline [1] [1] };
    geomvert = [ 1   0.5 1        1        1        1        1        0.3 1        0.3 1 ]; 
    evalstr  = 'set(findobj(gcf, ''tag'', ''title''), ''fontsize'', 16);';

    result = inputgui( 'geometry', geom, 'uilist', uilist, 'helpcom', 'pophelp(''plugin_installstartup'')', 'title', 'Install popular plugins', 'geomvert', geomvert, 'eval', evalstr);
    if length(result) == 0 return; end;
    