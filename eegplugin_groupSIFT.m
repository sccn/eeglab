% eegplugin_groupSIFT(): A plugin for EEGLAB to integrate SIFT results in
%                        the group-level. It's one of distribution packages
%                        of Nima Bigdely-Shamlo's Network Projection.
%
% Author: Makoto Miyakoshi, SCCN,INC,UCSD
%
% History
% 04/04/2023 Makoto. startup:on;study:on added to force EEGLAB to show it up in the GUI menu.

% Copyright (C) 2016, Makoto Miyakoshi (mmiyakoshi@ucsd.edu) , SCCN,INC,UCSD
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

function vers = eegplugin_groupSIFT(fig, try_strings, catch_strings)

% 04/04/2023.
vers = '0.31';

if nargin < 3
    error('eegplugin_groupSIFT requires 3 arguments');
end;

% create a highLevelMenu
highLevelMenu = findobj(fig, 'tag', 'tools');
set(highLevelMenu, 'UserData', 'startup:on;study:on'); % This unlocks 'Tools' menu without loading .set data haha.
submenu       = uimenu(highLevelMenu, 'label', 'groupSIFT','separator','on');

% add submenu
uimenu( submenu, 'label', '1.Run SIFT batch',                     'callback', 'pop_groupSIFT_runSiftBatch');
uimenu( submenu, 'label', '2.Validate AR models',                 'callback', 'pop_groupSIFT_validateArModels');
uimenu( submenu, 'label', '3.Convert to group anatomical ROIs',   'callback', 'pop_groupSIFT_convertToGroupAnatomicalRois');
uimenu( submenu, 'label', '4.Compute t-stats & p-values',         'callback', 'pop_groupSIFT_computeTstatsAndPvalues');
uimenu( submenu, 'label', '5 Show pre-selected ROIs',             'callback', 'pop_groupSIFT_showPreselectedRois');
uimenu( submenu, 'label', '6.View results & Export for movie',    'callback', 'pop_groupSIFT_viewResultsAndExportForMovie');