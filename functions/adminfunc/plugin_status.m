% plugin_status()- Given the name of a plugin(or part of it), returns 
%                  the status (e.g. [0]not installed, [1]installed, [2]deactivaded)                                                               
%  Inputs:
%         pluginname   - {string} Name (or part of it) of plugin of interest
%
% Optional inputs:
%        exactmatch  - [0,1] Force the function to look for an exact match of
%                       the name provided in the input 'pluginname' 
%                       {defaul: 0 = not exact match enforced}
%  Outputs:
%        
%        status         - [Vector 1xnumber of plugins]. 0 = Not installed,
%                         1 = Installed, 2 = Deactivated.
%       pluginnameout   - Name of the plugins in 'status'.
%       pluginstruct    - If the plugin(s) in 'pluginname' is installed, the function
%                         provide here a cell (one cell per plugin in 'status') array 
%                         of structures with the following fields:
%                         {plugin, version, foldername, funcname, status, versionfunc} 
% 
% Author: Ramon Martinez-Cancino SCCN, 2018
%         
% Copyright (C) 2018  Ramon Martinez-Cancino,INC, SCCN
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

function[status, pluginnameout, pluginstruct] = plugin_status(pluginname, varargin)
status = []; pluginstruct = []; pluginnameout = [];

try
    options = varargin;
    if ~isempty( varargin )
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else
        g = [];
    end
catch
    disp('plugin_status() error: calling convention {''key'', value, ... } error'); return;
end

try g.exactmatch;   catch, g.exactmatch  = 0;  end % Enforces exact match with the name provided

% Look in the plugin list for the plugin name provided
global PLUGINLIST
if ~isempty(PLUGINLIST)
    if g.exactmatch
        hitindx = find(strcmp(pluginname,{PLUGINLIST.plugin}));
        pluginnameout = pluginname;
    else
        hitindx = find(~cellfun(@isempty,strfind(lower({PLUGINLIST.plugin}),lower(pluginname))));
    end
    
    if isempty(hitindx)
        status = 0;
        return;
    else
        pluginstatus = {PLUGINLIST(hitindx).status};
        for i =1:length(pluginstatus)
            if strcmp(pluginstatus{i}, 'ok')
                status(i) = 1;
            else
                status(i) = 0;
            end
        end
        pluginnameout = {PLUGINLIST(hitindx).plugin};
        pluginstruct = PLUGINLIST(hitindx);
    end
end