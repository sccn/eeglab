% pluguin_uiupdate() - support function for plugin_menu

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

function plugin_uiupdate(fig)

tmpobj = get(fig, 'userdata');
listboxVal = get(findobj(fig, 'tag', 'pluginlist'), 'value');

allPlugins = tmpobj.allplugins;
selectedPlugins = tmpobj.selectedplugins;
selection       = tmpobj.selection;
if length(listboxVal) == 1
    lastSelect = listboxVal;
else
    lastSelect = setdiff(listboxVal, selection);
    
    % check if selection is consistent
    if ~all([selectedPlugins([ listboxVal ]).installed]) && ~all(~[selectedPlugins([ listboxVal ]).installed] )
        listboxVal = lastSelect(1);
        set(findobj(fig, 'tag', 'pluginlist'), 'value', listboxVal);
    end
end
if isempty(lastSelect), lastSelect = listboxVal(1); end
if length(lastSelect) > 1, lastSelect = lastSelect(1); end

tmpstr = selectedPlugins(lastSelect).description; 

% set text
hDescription = findobj(fig, 'tag', 'description');
hTags        = findobj(fig, 'tag', 'tags');
hStatus      = findobj(fig, 'tag', 'status');
hSize        = findobj(fig, 'tag', 'size');
tmpstr = textwrap(hDescription, {tmpstr}); 
set(hDescription, 'string', tmpstr);
set(hTags       , 'string', selectedPlugins(lastSelect).rawtags );
set(hSize       , 'string', [num2str(ceil(selectedPlugins(lastSelect).size/100)/10) ' MB'] );

% installed status
if selectedPlugins(lastSelect).installed
    if selectedPlugins(lastSelect).installorupdate
        set(hStatus, 'string', ['version ' selectedPlugins(lastSelect).currentversion ' installed - an update is available'] );
    else
        set(hStatus , 'string', 'installed' );
    end
else
    set(hStatus, 'string', 'not installed');
end

% set buttons
set(findobj(fig, 'tag', 'install'), 'enable', 'on');
if selectedPlugins(lastSelect).installed
    set(findobj(fig, 'tag', 'rmbut'), 'enable', 'on');
else
    set(findobj(fig, 'tag', 'rmbut'), 'enable', 'off');
end
if ~isempty(selectedPlugins(lastSelect).webdoc)
    set(findobj(fig, 'tag', 'documentation'), 'enable', 'on', 'callback', ['web(''' selectedPlugins(lastSelect).webdoc ''', ''-browser'');' ]);
else
    set(findobj(fig, 'tag', 'documentation'), 'enable', 'off')
end
set(findobj(fig, 'tag', 'rating'),  'callback', ['web(''' selectedPlugins(lastSelect).webrating ''', ''-browser'');' ]);

tmpobj.selection = listboxVal;
set(fig, 'userdata', tmpobj);
