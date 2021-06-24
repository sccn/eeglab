% eeg_hedremoveunicode() - remove unicode characters in HED strings.
%
% Usage:
%   >> EEGOUT = eeg_hedremoveunicode(EEGIN);
%
% Input:
%   EEGOUT - EEG dataset structure or array of structures
%
% Output:
%   EEGOUT - EEG dataset structure or array of structures
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, 2021

% Copyright (C) 2021 Arnaud Delorme
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

function EEG = eeg_hedremoveunicode(EEG)

unicodeList = {};
if isfield(EEG.event, 'usertags')
    for iEvent = 1:length(EEG.event)
        [EEG.event(iEvent).usertags,unicodeList{end+1}] = removeunicode(EEG.event(iEvent).usertags);
    end
end
if isfield(EEG.urevent, 'usertags')
    for iEvent = 1:length(EEG.urevent)
        [EEG.urevent(iEvent).usertags,unicodeList{end+1}] = removeunicode(EEG.urevent(iEvent).usertags);
    end
end

% check ETC structure
if isfield(EEG.etc, 'tags')
    if isfield(EEG.etc.tags, 'map')
        if isfield(EEG.etc.tags.map, 'values')
            for iMap = 1:length(EEG.etc.tags.map)
                for iTag1 = 1:length(EEG.etc.tags.map(iMap).values)
                    [EEG.etc.tags.map(iMap).values(iTag1).tags,unicodeList{end+1}] = removeunicode(EEG.etc.tags.map(iMap).values(iTag1).tags);
                end
            end
        end
    end
end

% show character list
unicodeList = unicodeList(~cellfun(@isempty, unicodeList));
listUni = unique(unicodeList);
listUni = [ listUni{:} ];
for iUni = 1:length(listUni)
    fprintf('Found unicode char: %s (%d)\n', listUni(iUni), listUni(iUni)+0);
end

% remove data and attempt saving v6
TMP = EEG;
TMP.data = [];
TMP.times = [];
TMP.icaact = [];
try
    save('tmptmpfile.mat', '-mat', '-v6', 'TMP');
    delete('tmptmpfile.mat');
catch
    disp('Could not remove all unicode');
end

function [str, unicodechars] = removeunicode(str)
unicodechars = [];
if iscell(str)
    for iStr = 1:length(str)
        [str{iStr}, unicodechars] = removeunicode(str{iStr});
    end
else
    posUni = find(str > 255);
    unicodechars = str(posUni);
    str(posUni) = ' ';
end


