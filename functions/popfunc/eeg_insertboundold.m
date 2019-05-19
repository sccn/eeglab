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
%   indnew     - array of indices returning new event index for any old 
%                (input eventin) event index
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

function [eventout,indnew] = eeg_insertbound( eventin, pnts, regions, lengths);
    
    if nargin < 3
        help eeg_insertbound;
        return;
    end
    if size(regions,2) ~= 1 && exist('lengths') ~= 1
        lengths = regions(:,2)-regions(:,1)+1;
        regions = regions(:,1);
    end
    if exist('lengths') ~= 1
        lengths = zeros(size(regions));
    end
    
    % resort events
    if ~isempty(eventin)
        [~,sortedEvents] = sort([eventin.latency]);
        eventin = eventin(sortedEvents);
    end
    
    if length(regions)
        fprintf('eeg_insertbound(): %d boundary (break) events added.\n', size(regions, 1));
    else 
        return;
    end

    % recompute latencies of boundevents (in new dataset)
    % ---------------------------------------------------
    [regions tmpsort] = sort(regions);
    lengths           = lengths(tmpsort);
    boundevents       = regions(:,1)-0.5;
    
    % sort boundevents by decreasing order (otherwise bug in new event index)
    % ------------------------------------
    boundevents = boundevents(end:-1:1);
    lengths     = lengths    (end:-1:1);
    eventout    = eventin;
    indnew      = 1:length(eventin);
    allnest     = [];
    countrm     = 0;
	for tmpindex = 1:length(boundevents) % sorted in decreasing order
        if boundevents(tmpindex) >= 0.5 && boundevents(tmpindex) <= pnts
                            
            % find event succeding boundary to insert event 
            % at the correct location in the event structure
            % ----------------------------------------------
            if ~isempty(eventout) && isfield(eventout, 'latency')
                alllats   = [ eventout.latency ] - boundevents(tmpindex);
                tmpind    = find( alllats >= 0 );
                [tmp tmpind2 ] = min(alllats(tmpind));
                tmpind2        = tmpind(tmpind2);
            else
                tmpind2 = [];
            end
            
            % insert event at tmpind2
            % -----------------------
            if ~isempty(tmpind2)
                eventout(end+1).type      = 'boundary';
                tmp = eventout(end);
                eventout(tmpind2+1:end) = eventout(tmpind2:end-1);
                eventout(tmpind2) = tmp;
                indnew(tmpind2:end) = indnew(tmpind2:end)+1;
            else
                tmpind2 = length(eventout)+1;
                eventout(tmpind2).type     = 'boundary';
            end
            eventout(tmpind2).latency  = boundevents(tmpindex);
            eventout(tmpind2).duration = lengths(tmpindex); % just to create field
            
            [ tmpnest addlength ] = findnested(eventout, tmpind2);
            
            % recompute latencies and remove events in the rejected region
            % ------------------------------------------------------------
            eventout(tmpnest) = [];
            countrm           = countrm+length(tmpnest);
            for latind = tmpind2+1:length(eventout)
                eventout(latind).latency = eventout(latind).latency-lengths(tmpindex);
            end
            
            % add lengths of previous events (must be done after above)
            % ---------------------------------------------------------
            eventout(tmpind2).duration = lengths(tmpindex)+addlength;                
            if eventout(tmpind2).duration == 0, eventout(tmpind2).duration=NaN; end
        
        end; 
	end

    if countrm > 0
        fprintf('eeg_insertbound(): event latencies recomputed and %d events removed.\n', countrm);
    end

    
% look for nested events
% retrun indices of nested events and
% their total length
% -----------------------------------
function [ indnested, addlen ] = findnested(event, ind);
    indnested = [];
    addlen = 0;
    tmpind = ind+1;

    while tmpind <= length(event) & ...
        event(tmpind).latency < event(ind).latency+event(ind).duration
        if strcmpi(event(tmpind).type, 'boundary')
            if ~isempty( event(tmpind).duration )
                addlen    = addlen + event(tmpind).duration;
                % otherwise old event duration or merge data discontinuity
            end
        end
        indnested = [ indnested tmpind ];
        tmpind = tmpind+1;
    end
    
% remove urevent and recompute indices
% THIS FUNCTION IS DEPRECATED
% ------------------------------------
function [event, urevent] = removenested(event, urevent, nestind);
    
    if length(nestind) > 0
        fprintf('eeg_insertbound() debug msg: removing %d nested urevents\n', length(nestind));
        nestind = sort(nestind);
        urind = [ event.urevent ]; % this must not be done in the loop
                                             % since the indices are dyanmically updated
    end
    
    for ind = 1:length(nestind)
        % find event urindices higher than the urevent to suppress
        % --------------------------------------------------------
        tmpind = find( urind > nestind(ind) );
        for indevent = tmpind
            event(indevent).urevent = event(indevent).urevent-1;
        end;    
    end
    
    urevent(nestind) = [];
    
