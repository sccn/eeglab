% eeg_insertbound() - insert boundary event in an EEG event structure.
%
% Usage:
%       >> [eventout newind] = eeg_insertbound( eventin, pnts, ...
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
%   newind     - Indices of the new events
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

function [eventin, newind] = eeg_insertbound( eventin, pnts, regions, lengths)
    
    if nargin < 3
        help eeg_insertbound;
        return;
    end
    regions = round(regions);
    regions(regions < 1) = 1;
    regions(regions > pnts) = pnts;
    for i=2:size(regions,1)
        if regions(i-1,2) >= regions(i,1)
            regions(i,1) = regions(i-1,2)+1;
        end
    end

    if ~isempty(regions)
        fprintf('eeg_insertbound(): %d boundary (break) events added.\n', size(regions, 1));
    else 
        return;
    end

    % recompute latencies of boundevents (in new dataset)
    % ---------------------------------------------------
    [tmp, tmpsort] = sort(regions(:,1));
    regions        = regions(tmpsort,:);
    lengths = regions(:,2)-regions(:,1)+1;
    
    if ~isempty(eventin)
         eventLatencies = [ eventin.latency ]; 
    else eventLatencies = [];
    end
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
        %if regions(iRegion,1) % do not remove first event
            eventin(end+1).type   = 'boundary';
            eventin(end).latency  = regions(iRegion,1)-sum(lengths(1:iRegion-1))-0.5;
            eventin(end).duration = lengths(iRegion,1)+addlength;
        %end
    end

    % copy latencies
    % --------------
    for iEvent = 1:oriLen
        eventin(iEvent).latency = newEventLatencies(iEvent);
    end
    eventin(rmEvent) = [];
    
    % resort events
    % -------------
    if ~isempty(eventin) && isfield(eventin, 'latency')
%       eventin([ eventin.latency ] < 1) = [];
        eventin([ eventin.latency ] < 0) = [];
        alllatencies = [ eventin.latency ];
        [tmp, sortind] = sort(alllatencies);
        eventin = eventin(sortind);
        newind = sortind(oriLen+1:end);
    end
    
    if ~isempty(rmEvent)
        fprintf('eeg_insertbound(): boundary events inserted and %d events removed.\n', length(rmEvent));
    end

    
% look for nested events
% retrun indices of nested events and
% their total length
% -----------------------------------
function [ indEvents, addlen ] = findnested(event, eventlat, region)
    indEvents = find( eventlat > region(1) & eventlat < region(2));

    if ~isempty(event) && isfield(event,'type') && ischar(event(1).type) && isfield(event, 'duration')
        boundaryInd = strmatch('boundary', { event(indEvents).type });
        addlen      = sum( [ event(indEvents(boundaryInd)).duration ] );
    else
        addlen = 0;
    end
    
