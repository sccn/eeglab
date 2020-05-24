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

function plugin = plugin_getweb(type, pluginOri, varargin)
plugin = [];

if nargin < 1, help plugin_getweb; return; end
if nargin < 2, pluginOri = []; end

% convert plugin list format if necessary
if isfield(pluginOri, 'plugin'), pluginOri = plugin_convert(pluginOri); end

% retreiving statistics
try
    disp( [ 'Retreiving download statistics...' ] );
    [stats, status] = plugin_urlread('http://sccn.ucsd.edu/eeglab/plugin_uploader/plugin_getcountall_nowiki.php');
    stats = textscan(stats, '%s%d%s%s%f%d%s%s%s%s%s%f', 'delimiter', char(9));
    if length(unique(cellfun(@length, stats))) > 1
        disp('Issue with retrieving statistics for extensions');
        return;
    end
catch
    disp('Cannot connect to the Internet to retrieve statistics for extensions');
    return
end

if status == 0
    disp('Cannot connect to the Internet to retrieve extension list');
    return
end

% decode stats into plugins
% -------------------------
if ~isempty(pluginOri)
     currentNames = lower({ pluginOri.name });
else currentNames = {};
end
for iRow = 1:length(stats{1})
    plugin(iRow).name      = stats{1}{iRow};
    plugin(iRow).downloads = stats{2}(iRow);
    plugin(iRow).version   = stats{3}{iRow};
    plugin(iRow).zip       = stats{4}{iRow};
    plugin(iRow).rating    = stats{5}(iRow);
    plugin(iRow).numrating = stats{6}(iRow);
    plugin(iRow).description  = stats{7}{iRow};
    plugin(iRow).rawtags      = stats{8}{iRow};
    if ~isempty(plugin(iRow).rawtags)
        tmpTags = textscan(plugin(iRow).rawtags, '%s', 'delimiter', ',');
        plugin(iRow).tags         = tmpTags{1}';
    end
    plugin(iRow).contactname  = stats{9}{iRow};
    plugin(iRow).contactemail = stats{10}{iRow};
    plugin(iRow).webdoc    = stats{11}{iRow};
    plugin(iRow).size      = stats{12}(iRow);
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
