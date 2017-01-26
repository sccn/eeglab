% eeg_insertbound() - insert boundary event in an EEG event structure.
%
% Usage:
%       >> [eventout indnew] = eeg_insertbound( eventin, pnts, ...
%                                                abslatency, duration);
% Required Inputs:
%   eventin    - EEGLAB event structure (EEG.event)
%   pnts       - data points in EEG dataset (EEG.pnts * EEG.trials)
%   abslatency - absolute latency of regions in original dataset. Can
%                also be an array of [beg end] latencies with one row
%                per region removed. Then 'lengths' argument is ignored.
% Optional Inputs:
%   lengths    - lengths of removed regions
%
% Outputs:
%   eventout   - EEGLAB event output structure with added boundaries
%   indnew     - Indices of the new events
%
% Notes:
%   This function performs the following: 
%   1) add boundary events to the 'event' structure; 
%        remove nested boundaries; 
%        recompute the latencies of all events.
%   2) all latencies are given in (float) data points. 
%        e.g., a latency of 2000.3 means 0.3 samples (at EEG.srate)
%              after the 2001st data frame (since first frame has latency 0).
% 
% Author: Arnaud Delorme and Hilit Serby, SCCN, INC, UCSD, April, 19, 2004
%
% See also: eeg_eegrej(), pop_mergeset()

% Copyright (C) 2004 Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
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

function [eventin, newind] = eeg_insertbound( eventin, pnts, regions, lengths)
    
    if nargin < 3
        help eeg_insertbound;
        return;
    end;
    regions = round(regions);
    regions(regions < 1) = 1;
    regions(regions > pnts) = pnts;
    for i=2:size(regions,1)
        if regions(i-1,2) >= regions(i,1)
            regions(i,1) = regions(i-1,2)+1;
        end;
    end;

    if ~isempty(regions)
        fprintf('eeg_insertbound(): %d boundary (break) events added.\n', size(regions, 1));
    else 
        return;
    end;

    % recompute latencies of boundevents (in new dataset)
    % ---------------------------------------------------
    [tmp, tmpsort] = sort(regions(:,1));
    regions        = regions(tmpsort,:);
    lengths = regions(:,2)-regions(:,1)+1;
    
    if ~isempty(eventin)
         eventLatencies = [ eventin.latency ]; 
    else eventLatencies = [];
    end;
    newEventLatencies = eventLatencies;
    oriLen            = length(eventin);
    rmEvent           = [];
	for iRegion = 1:size(regions,1) % sorted in decreasing order

        % find event succeding boundary to insert event
        % at the correct location in the event structure
        % ----------------------------------------------
        tmpind    = find( eventLatencies - regions(iRegion,1) > 0 );
        newEventLatencies(tmpind) = newEventLatencies(tmpind)-lengths(iRegion);
        
        % insert event
        % ------------
        [tmpnest, addlength ]  = findnested(eventin, eventLatencies, regions(iRegion,:));
        rmEvent = [ rmEvent tmpnest ];
        if regions(iRegion,1)>1
            eventin(end+1).type   = 'boundary';
            eventin(end).latency  = regions(iRegion,1)-sum(lengths(1:iRegion-1))-0.5;
            eventin(end).duration = lengths(iRegion,1)+addlength;
        end;
    end

    % copy latencies
    % --------------
    for iEvent = 1:oriLen
        eventin(iEvent).latency = newEventLatencies(iEvent);
    end;
    eventin(rmEvent) = [];
    
    % resort events
    % -------------
    if ~isempty(eventin) && isfield(eventin, 'latency')
        eventin([ eventin.latency ] < 1) = [];
        alllatencies = [ eventin.latency ];
        [tmp, sortind] = sort(alllatencies);
        eventin = eventin(sortind);
        newind = sortind(oriLen+1:end);
    end;
    
    if ~isempty(rmEvent)
        fprintf('eeg_insertbound(): event latencies recomputed and %d events removed.\n', length(rmEvent));
    end;

    
% look for nested events
% retrun indices of nested events and
% their total length
% -----------------------------------
function [ indEvents, addlen ] = findnested(event, eventlat, region)
    indEvents = find( eventlat > region(1) & eventlat < region(2));

    if ~isempty(event) && isstr(event(1).type) && isfield(event, 'duration')
        boundaryInd = strmatch('boundary', { event(indEvents).type });
        addlen      = sum( [ event(indEvents(boundaryInd)).duration ] );
    else
        addlen = 0;
    end;
    