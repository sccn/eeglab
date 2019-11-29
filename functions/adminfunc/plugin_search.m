% plugin_search() - support function for plugin_menu

% Copyright (C) 2019 Arnaud Delorme
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

function plugin_search(fig)

tmpobj = get(fig, 'userdata');
allPlugins = tmpobj.allplugins;

uistyle      = get(findobj(fig, 'tag', 'search'), 'style');
if strcmpi(uistyle, 'pushbutton')
    set(findobj(fig, 'tag', 'search'), 'style', 'edit', 'string', '');
    return;
end

searchString = lower(get(findobj(fig, 'tag', 'search'), 'string'));
if isempty(searchString), return; end
searchString = textscan(searchString, '%s');
searchString = searchString{1};

% reset filters
set(findobj(fig, 'tag', 'filter1'), 'value', 1);
set(findobj(fig, 'tag', 'filter2'), 'value', 1);

% search
if isempty(searchString)
    selectedPlugins = allPlugins;
else
    selectedPluginsIndex = [];
    for iSearch = 1:length(searchString)
        res = strfind( { allPlugins.strsearch }, searchString{iSearch});
        if isempty(selectedPluginsIndex)
            selectedPluginsIndex = find(~cellfun(@isempty, res));
        else
            selectedPluginsIndex = intersect(selectedPluginsIndex, find(~cellfun(@isempty, res)));
        end
    end
    selectedPlugins = allPlugins(selectedPluginsIndex);
end

% update GUI
set(findobj(fig, 'tag', 'pluginlist'), 'string', { selectedPlugins.text }, 'value', []);
tmpobj.selection = [];
tmpobj.selectedplugins = selectedPlugins;

set(fig, 'userdata', tmpobj);
