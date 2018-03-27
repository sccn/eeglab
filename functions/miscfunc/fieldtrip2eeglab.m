% fieldtrip2eeglab - convert Fieldtrip structures to EEGLAB dataset
%
% EEG = fieldtrip2eeglab(header, data, evt);
%
% Inputs:
%    header   - Fieldtrip data header 
%    data     - Fieldtrip raw data
%    evt      - Fieldtrip event structure (optional)
%
% Output:
%    EEG     - EEGLAB structure
%
% Author: Arnaud Delorme, UCSD

function EEG = fieldtrip2eeglab(header,data,evt)

if nargin < 3
    evt = [];
end

EEG = pop_fileio(header, data, evt);
