% plugin_askinstall() - install EEGLAB plugin from a GUI or comamnd line
%                       call.
% Usage:
%  plugin_install(pluginName, pluginFunc, force);
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

function installRes = plugin_askinstall(pluginName, pluginFunc, forceInstall)

if nargin < 3, forceInstall = false; end

if nargin < 2 || ~exist(char(pluginFunc))
    
    if ~forceInstall
        db = dbstack;
        if length(db) > 2 && ~strcmpi(db(end).name, 'checkouteeglab.m')
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
    plugins = plugin_getweb('import', []);
    indPlugin = strmatch(lower(pluginName), lower({ plugins.name }), 'exact');
    if isempty(indPlugin), 
        plugins = plugin_getweb('process', []);
        indPlugin = strmatch(lower(pluginName), lower({ plugins.name }));
        if isempty(indPlugin), 
            error([ pluginName ' extension not found' ]); 
        end
    end
    result = plugin_install(plugins(indPlugin(1)).zip, plugins(indPlugin(1)).name, plugins(indPlugin(1)).version, forceInstall);
    if result == 1, installRes = 1; end
    
    evalin('base', 'eeglab rebuild');
else
    installRes = 1;
end
