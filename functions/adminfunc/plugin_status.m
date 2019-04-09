% plugin_status()- Given the name of a plugin(or part of it), returns 
%                  the status (e.g. [0]not installed, [1]installed)                                                               
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
%                         1 = Installed
%       pluginnameout   - Name of the plugins in 'status'.
%       pluginstruct    - If the plugin(s) in 'pluginname' is installed, the function
%                         provide here a cell (one cell per plugin in 'status') array 
%                         of structures with the following fields:
%                         {plugin, version, foldername, funcname, status} 
% 
% Author: Ramon Martinez-Cancino SCCN, 2018
%         
% Copyright (C) 2018  Ramon Martinez-Cancino,INC, SCCN
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

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
