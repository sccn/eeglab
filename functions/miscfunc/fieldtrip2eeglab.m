% fieldtrip2eeglab - convert Fieldtrip structures to EEGLAB dataset
%
% EEG = fieldtrip2eeglab(header, rawdata, evt);
% EEG = fieldtrip2eeglab(data);
%
% Inputs:
%    header   - Fieldtrip data header 
%    rawdata  - Fieldtrip raw data
%    evt      - Fieldtrip event structure (optional)
%    data     - Fieldtrip data out of ft_preprocessing
%
% Output:
%    EEG     - EEGLAB structure
%
% Author: Arnaud Delorme, UCSD

% Copyright (C) Arnaud Delorme, UCSD 2018
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

function EEG = fieldtrip2eeglab(header,data,evt)

if nargin < 3
    evt = [];
end

if nargin >= 2
    EEG = pop_fileio(header, data, evt);
else
    EEG = pop_fileio(header.hdr, header.trial);
    if iscell(EEG.data) && length(EEG.data) == 1
        EEG.data = EEG.data{1};
    elseif iscell(EEG.data)
        len = cellfun(@(x)size(x,2), EEG.data);
        if length(unique(len)) > 1
            error('Cannot convert epochs of different length');
        end
        EEG.data = [ EEG.data{:} ];
        EEG.pnts = len(1);
    else
        error('Unknown fieldtrip data format');
    end
end   