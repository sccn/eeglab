% std_chaninds() - look up channel indices in a STUDY
%
% Usage:
%         >> inds = std_chaninds(STUDY,  channames);
%         >> inds = std_chaninds(EEG, channames);
%         >> inds = std_chaninds(chanlocs, channames);
% Inputs:
%         STUDY    - studyset structure containing a changrp substructure.
%         EEG      - EEG structure containing channel location structure
%         chanlocs - channel location structure
%     channames - [cell] channel names
%
% Outputs:
%       inds - [integer array] channel indices
%
% Author: Arnaud Delorme, CERCO, 2006-

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function finalinds = std_chaninds(instruct, channames);

    finalinds   = [];
    if isfield(instruct, 'chanlocs')
        EEG = instruct;
        tmpchanlocs = EEG.chanlocs;
        tmpallchans = lower({ tmpchanlocs.labels });
    elseif isfield(instruct, 'filename')
        tmpallchans = lower({ instruct.changrp.name });
    else
        tmpallchans = instruct;
    end
    if ~iscell(channames), channames = { channames }; end
    
    if isempty(channames), finalinds = [1:length(tmpallchans)]; return; end
    for c = 1:length(channames)
        if isnumeric(channames)
            chanind = channames(c);
        else
            chanind = strmatch( lower(channames{c}), tmpallchans, 'exact');
            if isempty(chanind), warning(sprintf([ 'Channel "%s" and maybe others was not' 10 'found in pre-computed data file' ], channames{c})); end
        end
        if length(chanind) > 1
            error(sprintf('Duplicate channel label %s - fix the issue and try again', channames{c}));
        end
        finalinds   = [ finalinds chanind ];
    end
