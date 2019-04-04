% eeg_latencyur() - transform latency of sample point in the continuous 
%                data into latencies in the transformed dataset.
%
% Usage:
%   >> lat_out = eeg_latencyur( events, lat_in);
%
% Inputs:
%   events     - event structure. If this structure contain boundary
%                events, the length of these events is added to restore
%                the original latency from the relative latency in 'lat_in'
%   lat_in     - sample latency (in point) in the original urEEG.
%
% Outputs:
%   lat_out     - output latency
%
% Note: the function that finds the original (ur) latency in the original
%       dataset using latencies in the current dataset is called 
%       eeg_urlatency()
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2011-
% 
% See also: eeg_urlatency()

% Copyright (C) 2011 Arnaud Delorme, SCCN, INC, UCSD, arno@ucsd.edu
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

function latout = eeg_latencyur( events, latin );
    
    if nargin < 2
        help eeg_latencyur;
        return;
    end

    boundevents = { events.type };
    latout      = latin;
    if ~isempty(boundevents) && ischar(boundevents{1})
        indbound = strmatch('boundary', boundevents);
        
        if isfield(events, 'duration') && ~isempty(indbound)
            for index  = indbound'
                lowerVals = find(latout > events(index).latency);
                latout(lowerVals) = latout(lowerVals)-events(index).duration;
            end
        end
    end
    return;
    
    % build array of 0 and 1 (0 no data)
    boundevents = { events.type };
    latout      = latin;
    if ~isempty(boundevents) && ischar(boundevents{1})
        indbound = strmatch('boundary', boundevents);
        
        if isfield(events, 'duration') && ~isempty(indbound)
            currentadd = 0;
            points     = ones(1, events(end).latency+sum([events(indbound').duration])); % arrays of 1
            for index  = indbound'
                currentlat = events(index).latency+currentadd;
                currentdur = events(index).duration;
                points(round(currentlat):round(currentlat+currentdur)) = 0;
                currentadd = currentadd + currentdur;
            end
        end
    end
    8;
    
    
    
    
