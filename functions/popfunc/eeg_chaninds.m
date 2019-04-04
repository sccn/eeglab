% std_chaninds() - look up channel indices in a EEG structure
%
% Usage:
%         >> inds = std_chaninds(EEG, channames);
% Inputs:
%         EEG   - EEG structure containing a chanlocs substructure.
%                 the chanlocs structure may also be used as input.
%     channames - [cell] channel names. May also be a string containing
%                 one or several channel names.
%
% Outputs:
%       inds - [integer array] channel indices
%
% Author: Arnaud Delorme, CERCO, 2009-

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

function finalinds = eeg_chaninds(EEG, channames, errorifnotfound);

    if nargin < 2
        help eeg_chaninds;
        return;
    end
    if nargin < 3
        errorifnotfound = 1;
    end
    
    if isfield(EEG, 'chanlocs')
         chanlocs = EEG.chanlocs;
    else chanlocs = EEG;
    end
    
    % decode string if necessary
    % --------------------------
    if ischar(channames)
        channames = parsetxt( channames );
    end

    finalinds   = [];
    if isempty(chanlocs)
         tmpallchans = [];
    else tmpallchans = lower({ chanlocs.labels });
    end
    if isempty(channames), finalinds = [1:length(chanlocs)]; return; end
    for c = 1:length(channames)
        chanind = strmatch( lower(channames{c}), tmpallchans, 'exact');
        if isempty(chanind), 
            chanind = str2double(channames{c});
            if isnan(chanind), chanind = []; end
            if errorifnotfound && isempty(chanind)
                error(sprintf('Channel %s not found', channames{c})); 
            end
        end
        finalinds   = [ finalinds chanind ];
    end
    finalinds = sort(finalinds);
