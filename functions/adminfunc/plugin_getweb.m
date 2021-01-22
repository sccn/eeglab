% plugin_getweb - support function to get plugin information from the web

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

function [plugin,eeglabVersionStatus] = plugin_getweb(type, pluginOri, varargin)
plugin = [];
eeglabVersionStatus = [];

if nargin < 1, help plugin_getweb; return; end
if nargin < 2, pluginOri = []; end

% convert plugin list format if necessary
if isfield(pluginOri, 'plugin'), pluginOri = plugin_convert(pluginOri); end

% retreiving statistics
eeglab_options;
try
    disp( [ 'Retreiving download statistics...' ] );
    if exist('OCTAVE_VERSION', 'builtin') == 0
        [plugin, status] = plugin_urlread([ 'http://sccn.ucsd.edu/eeglab/plugin_uploader/plugin_getcountall_nowiki_json.php?type=' type '&upload=' num2str(option_showpendingplugins)]);
    else
        [plugin, status] = urlread([ 'http://sccn.ucsd.edu/eeglab/plugin_uploader/plugin_getcountall_nowiki_json.php?type=' type '&upload=' num2str(option_showpendingplugins)]);
    end
    if isempty(plugin)
        disp('Issue with retrieving statistics for extensions');
        return;
    end
    try
        plugin = jsondecode(plugin);
    catch
        disp('Issue with decoding plugin information, Octave 7.x required to decode JSON strings');
    end
catch
    disp('Cannot connect to the Internet to retrieve statistics for extensions');
    return
end

if status == 0
    disp('Cannot connect to the Internet to retrieve extension list');
    return
end

%% rename fields for backward compabitiligy
renameField = { 'plugin' 'name';
                'count'  'downloads';
                'curversion' 'version';
                'link'   'zip' };

if ~isempty(pluginOri)
    currentNames = lower({ pluginOri.name });
else
    currentNames = {};
end
for iRow = 1:length(plugin)
    
    % rename fields
    for iField = 1:size(renameField, 1)
        plugin(iRow).(renameField{iField, 2}) = plugin(iRow).(renameField{iField, 1});
    end
    
    % decode tags
    plugin(iRow).rawtags = plugin(iRow).tags;
    if ~isempty(plugin(iRow).rawtags)
        tmpTags = textscan(plugin(iRow).rawtags, '%s', 'delimiter', ',');
        plugin(iRow).tags = tmpTags{1}';
    end
    
    plugin(iRow).numrating =  str2double(plugin(iRow).numrating);
    plugin(iRow).rating    =  str2double(plugin(iRow).rating);
    plugin(iRow).critical  =  str2double(plugin(iRow).critical);
    plugin(iRow).downloads =  str2double(plugin(iRow).downloads);
    plugin(iRow).size      =  sscanf(plugin(iRow).size, '%f'); % Only numeric part is taken, possible KB or MB additions are ignored
    plugin(iRow).webrating = [ 'https://sccn.ucsd.edu/eeglab/plugin_uploader/simplestar.php?plugin=' plugin(iRow).name '&version=' plugin(iRow).version ];
    
    % match with existiting plugins
    indMatch = strmatch(lower(plugin(iRow).name), currentNames, 'exact');
    if isempty(indMatch)
        plugin(iRow).currentversion  = '-';
        plugin(iRow).installed       = 0;
        plugin(iRow).installorupdate = 1;
        plugin(iRow).status          = 'notinstalled';
    else
        if length(indMatch) > 1
            disp([ 'Warning: duplicate extension ' plugin(iRow).name ' instaled' ]); 
        end
        plugin(iRow).currentversion = pluginOri(indMatch).currentversion;
        plugin(iRow).foldername     = pluginOri(indMatch).foldername;
        plugin(iRow).status         = pluginOri(indMatch).status;
        plugin(iRow).installed      = 1;
        if strcmpi(plugin(iRow).currentversion, plugin(iRow).version)
            plugin(iRow).installorupdate = 0;
        else
            plugin(iRow).installorupdate = 1;
        end
    end
end
plugin = rmfield(plugin, renameField(:,1)');

% handle the special case of EEGLAB version
indEEGLAB = cellfun(@(x)isequal(x, 'eeglab'), lower( { plugin.name } ));
eeglabVersionStatus = plugin(indEEGLAB);
plugin(indEEGLAB) = [];

% remove plugins with no zip file
% -------------------------------
plugin(cellfun(@isempty, {plugin.zip })) = [];

function plugin = plugin_convert(pluginOri)

for iRow = 1:length(pluginOri)
    plugin(iRow).currentversion = pluginOri(iRow).version;
    plugin(iRow).foldername     = pluginOri(iRow).foldername;
    plugin(iRow).status         = pluginOri(iRow).status;
    plugin(iRow).name           = pluginOri(iRow).plugin;
    plugin(iRow).installed      = 1;
end
