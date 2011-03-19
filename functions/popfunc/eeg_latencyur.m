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
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function latout = eeg_latencyur( events, latin );
    
    if nargin < 2
        help eeg_latencyur;
        return;
    end;

    boundevents = { events.type };
    latout      = latin;
    if ~isempty(boundevents) & isstr(boundevents{1})
        indbound = strmatch('boundary', boundevents);
        
        if isfield(events, 'duration') & ~isempty(indbound)
            for index  = indbound'
                lowerVals = find(latout > events(index).latency);
                latout(lowerVals) = latout(lowerVals)-events(index).duration;
            end;
        end;
    end;
    return;
    
    % build array of 0 and 1 (0 no data)
    boundevents = { events.type };
    latout      = latin;
    if ~isempty(boundevents) & isstr(boundevents{1})
        indbound = strmatch('boundary', boundevents);
        
        if isfield(events, 'duration') & ~isempty(indbound)
            currentadd = 0;
            points     = ones(1, events(end).latency+sum([events(indbound').duration])); % arrays of 1
            for index  = indbound'
                currentlat = events(index).latency+currentadd;
                currentdur = events(index).duration;
                points(round(currentlat):round(currentlat+currentdur)) = 0;
                currentadd = currentadd + currentdur;
            end;
        end;
    end;
    8;
    
    
    
    
