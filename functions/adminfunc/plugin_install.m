% plugin_install() - install EEGLAB plugin. Called by plugin_askinstall().
%
% Usage:
%  plugin_install(zipfilelink, name, version, force);
%
% Inputs:
%  zipfilelink - [string] web link to zip file
%  name        - [string] name of the plugin
%  version     - [string] version of the plugin
%  force       - [boolean] force install (even if already installed)
%
% See also: plugin_askinstall()

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

function result = plugin_install(zipfilelink, name, version, forceInstall)

    if nargin < 4, forceInstall = false; end
    result = 1;

    % get plugin path
    % ---------------
    %version(find(version == '.'))  = '_';
    generalPluginPath = fullfile(fileparts(which('eeglab.m')), 'plugins');
    newPluginPath     = fullfile(generalPluginPath, [ name version ]);

    % download plugin
    % ---------------
    zipfile = 'tmp.zip';
%     [tmp zipfile ext] = fileparts(zipfilelink);
%     zipfile = [ zipfile ext ];
%     equalPos = find(zipfile == '=');
%     if ~isempty(equalPos) zipfile  = zipfile(equalPos(end)+1:end); end
    depth = length(dbstack);
    try
        pluginSize = plugin_urlsize(zipfilelink);
        pluginSizeStr = num2str(round(pluginSize/100000)/10);
        if pluginSize > 500000 && depth > 1 && ~forceInstall
             res = questdlg2( [ 'Extension ' name ' size is ' pluginSizeStr 'MB. Are you sure' 10 ...
                                'you want to download this extension?' ], 'Warning', 'No', 'Yes', 'Yes');
             if strcmpi(res, 'no'), fprintf([ 'Skipping ' name ' extension instalation\n' ]); result = -1; return; end;               
        end
    catch
        msg = [ 'Could not download extension. Host site might be' 10 ...
                'unavailable, too slow or you do not have permission' 10 ...
                'to write in the EEGLAB plugin folder. Try again' 10 ... 
                'just in case or use a faster connection.' 10 10 ...
                'Alternatively install the plugin manually by downloading' 10 ...
                'it at http://sccn.ucsd.edu/wiki/EEGLAB_Extensions_and_plug-ins' 10 ...
                'unziping it in the eeglab/plugin folder and restarting eeglab'];
        if ~forceInstall
            warndlg2(msg);
        else
            disp(msg);
        end
        result = -1;
        return;
    end
    disp([ 'Downloading extension ' name '(' pluginSizeStr 'Mb)...' ]);
    
    try
        plugin_urlread(['http://sccn.ucsd.edu/eeglab/plugin_uploader/plugin_increment.php?plugin=' name '&version=' version ]);
        plugin_urlwrite( zipfilelink, fullfile(generalPluginPath, zipfile));
    catch
        msg = [ 'Could not download extension. Host site might be' 10 ...
                'unavailable, too slow or you do not have permission' 10 ...
                'to write in the EEGLAB plugin folder. Try again' 10 ... 
                'just in case or use a faster connection.' 10 10 ...
                'Alternatively install the plugin manually by downloading' 10 ...
                'it at http://sccn.ucsd.edu/wiki/EEGLAB_Extensions_and_plug-ins' 10 ...
                'unziping it in the eeglab/plugin folder and restarting eeglab'];
        if ~forceInstall
            warndlg2(msg);
        else
            disp(msg);
        end
        result = -1;
        return;
    end

    % unzip plugin
    % ------------
    if ~exist(newPluginPath)
        mkdir(newPluginPath);
    else
        msg = [ 'Extension folder already exist ' newPluginPath 10 'Remove it manually before installing extension' ];
        if ~forceInstall
            warndlg2(msg);
        else
            disp(msg);
        end
        result = -1;
        return;
    end
    disp([ 'Unzipping extension file... ']);
    unzip(fullfile(generalPluginPath, zipfile), newPluginPath);
    
    disp('Cleaning up zip file...');
    delete(fullfile(generalPluginPath, zipfile));
    
    % seeing what is in the plugin and moving files if necessary
    % ----------------------------------------------------------
    pluginContent = dir(newPluginPath);
    if length(pluginContent) > 3
        return;
    end
    for index = 1:length(pluginContent)
        if ~strcmpi(pluginContent(index).name, '.') && ~strcmpi(pluginContent(index).name, '..')
            fullFolder = fullfile(newPluginPath, pluginContent(index).name);
            if exist(fullFolder) == 7 % folder detected
                % move files from folder
                movefile(fullfile(fullFolder, '*'), newPluginPath);
                rmdir(fullFolder, 's');
            end
        end
    end
    fprintf('Extension %s version %s now installed\n', name, version); 
    
