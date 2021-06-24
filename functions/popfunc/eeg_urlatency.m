% eeg_urlatency() - find the original (ur) latency of a time point in 
%                   the original continuous data.
%
% Usage:
%   >> lat_out = eeg_urlatency( event, lat_in );
%
% Inputs:
%   event      - event structure. If this structure contain boundary
%                events, the length of these events is added to restore
%                the original latency from the relative latency in 'lat_in'
%   lat_in     - relative latency in sample point
%
% Outputs:
%   lat_out     - output latency
%
% Note: the function that finds the latency in the current dataset using (ur)
%       original latencies as input is eeg_latencyur()
% 
% Author: Arnaud Delorme, SCCN, INC, UCSD, April, 15, 2004
%
% See also: eeg_latencyur()

% Copyright (C) 2004 Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
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

function latout = eeg_urlatency( events, latin )
    
    if nargin < 2
        help eeg_urlatency;
        return;
    end
    
    boundevents = { events.type };
    latout      = latin;
    if ~isempty(boundevents) && ischar(boundevents{1})
        try
            indbound = strmatch('boundary', boundevents);
        catch
            % crash in the unlikely case of numerical events
            return
        end
        
        if isfield(events, 'duration') && ~isempty(indbound)
            for index = indbound'
                tmpInds = find(events(index).latency < latin); % the find handles several input latencies
                latout(tmpInds) = latout(tmpInds) + events(index).duration;
            end
        elseif ~isempty(indbound) % boundaries but with no associated duration
            latout = NaN;
        end
    end
    
