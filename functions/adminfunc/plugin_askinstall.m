% plugin_askinstall() - install EEGLAB plugin from a GUI or command line
%                       call.
% Usage:
%  plugin_askinstall(pluginName, pluginFunc, force);
%
% Inputs:
%  pluginName - [string] name of the plugin
%  pluginFunc - [string] function belonging to the plugin
%  force       - [boolean] force install (even if already installed)
%
% Example:
%   Force install BVA-IO
%   plugin_askinstall('bva-io', [], true);

% Copyright (C) 2012- Arnaud Delorme
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

function installRes = plugin_askinstall(pluginName, pluginFunc, forceInstall)

if nargin < 3, forceInstall = false; end

if nargin < 2 || ~exist(char(pluginFunc))
    
    if ~forceInstall
        db = dbstack;
        if length(db) > 2 && ~strcmpi(db(end).file, 'checkouteeglab.m')
            error([ 'Cannot find ' pluginName ' extension, use EEGLAB Extension Manager to install it' ]);
        end

        installRes = 0;

        % check is deactivated
        try, PLUGINLIST = evalin('base', 'PLUGINLIST'); catch, PLUGINLIST = []; end
        if ~isempty(PLUGINLIST) && isfield(PLUGINLIST, 'plugin')
            indPlugin = strmatch(lower(pluginName), lower({ PLUGINLIST.plugin }), 'exact');
            if ~isempty(indPlugin) && strcmpi(PLUGINLIST(indPlugin(1)).status, 'deactivated')
                res = questdlg2( [ pluginName ' extension is de-activated. Do you want to reactivate it now?' ], [ pluginName ' extension installation' ], 'No', 'Yes', 'Yes' );
                if strcmpi(res, 'no'), return, end
                plugin_reactivate(PLUGINLIST(indPlugin(1)).foldername);
                evalin('base', 'eeglab rebuild');
                installRes = 1;
                return;
            end
        end

        % check for installing
        res = questdlg2( [ pluginName ' extension is not installed. Do you want to download it now?' ], [ pluginName ' extension installation' ], 'No', 'Yes', 'Yes' );
    else
        res = 'yes';
    end
    
    if strcmpi(res, 'no'), return, end
    try
        plugins = plugin_getweb('plugin_install', []);
    catch
        error('Issue with retrieving statistics for extensions, maybe check your connection');
    end
    if isempty(plugins)
        error('Cannot download the extension. Please check your internet connection');
    end
    indPlugin = strmatch(lower(pluginName), lower({ plugins.name }), 'exact');
    if isempty(indPlugin)
        error([ pluginName ' extension not found' ]);
    end
    result = plugin_install(plugins(indPlugin(1)).zip, plugins(indPlugin(1)).name, plugins(indPlugin(1)).version, forceInstall);
    if result == 1, installRes = 1; end
    
    evalin('base', 'eeglab rebuild');
else
    installRes = 1;
end
