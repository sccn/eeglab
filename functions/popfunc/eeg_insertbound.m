% eeg_insertbound() - insert boundary event in event structure.
%
% Usage:
%   >> EEGOUT = eeg_insertbound( EEG, latency, abslatency, lengths, checkevent);
%
% Inputs:
%   EEG        - EEG dataset
%   latency    - relative latency (in the event structure) of boundary 
%                event(s)
%   abslatency - absolute latency of regions in original dataset. Can
%                also be an array of [beg end] latency with one row
%                per region removed. Then 'lengths' is ignored.
%   lengths    - length of removed regions
%   checkevent - [0|1] check event consistency (1). Default is 1.
%
% Outputs:
%   EEG        - output EEG dataset with 'event' and 'urevent' fields
%                updated
%
% Notes:
%   1) This function performs the following: add boundary events to the 
%   'event' and 'urevent' structures. Update urevent indices in the
%   'event' structure. Add a length field to the 'urevent' boundaries.
%   2) all latencies are given in data point unit.
% 
% Author: Arnaud Delorme, SCCN, INC, UCSD, April, 19, 2004
%
% See also: eeg_eegrej(), pop_mergeset()

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.10  2004/05/06 21:52:43  arno
% debuging for mergeset
%
% Revision 1.9  2004/05/06 17:42:32  arno
% add checkevent input
% .,
%
% Revision 1.8  2004/05/05 01:54:40  arno
% removing debug message
%
% Revision 1.7  2004/05/05 01:50:17  arno
% same
%
% Revision 1.6  2004/05/05 01:49:25  arno
% debug regions input
%
% Revision 1.5  2004/05/05 01:43:11  arno
% don't know
%
% Revision 1.4  2004/05/04 23:19:52  arno
% typo
%
% Revision 1.3  2004/05/04 23:05:02  arno
% debug length
%
% Revision 1.2  2004/05/04 19:01:08  arno
% removing nested urevent boundarie
%
% Revision 1.1  2004/04/20 02:09:31  arno
% Initial revision
%
% Revision 1.1  2004/04/20 01:11:39  arno
% Initial revision
%

function EEG = eeg_insertbound( EEG, boundevents, regions, lengths, checkevent );
    
    if nargin < 2
        help eeg_insertbound;
        return;
    end;
    if nargin < 5
        checkevent = 1;
    end;
    if size(regions,2) ~= 1 & exist('lengths') ~= 1
        lengths = regions(:,2)-regions(:,1)+1;
        regions = regions(:,1);
    end;
    if exist('lengths') ~= 1
        lengths = zeros(size(regions));
    end;
    
    if isfield(EEG, 'urevent')
        newur = EEG.urevent;
    end;
    rmnested = [];
	for tmpindex = 1:length(boundevents)
        if boundevents(tmpindex) >= 1 & boundevents(tmpindex) <= EEG.pnts
            
            % insert event at the correct location in the urevent structure
            % -------------------------------------------------------------
            if isfield(EEG, 'urevent') & isfield(EEG.urevent, 'latency')
                % find event succeding boundary
                % ------------------------------
                urlatency = eeg_urlatency(EEG.urevent, regions(tmpindex)-0.5);
                alllats   = cell2mat( { newur.latency } ) - urlatency;
                tmpind    = find( alllats > 0 );
                [tmp tmpind2 ] = min(alllats(tmpind));
                tmpind2        = tmpind(tmpind2);
                
                % insert event at tmpind2
                % -----------------------
                if ~isempty(tmpind2)
                    newur(tmpind2+1:end+1) = newur(tmpind2:end);
                    newur(tmpind2).type  = 'boundary';
                else
                    newur(end+1).type  = 'boundary';
                    tmpind2 = length(newur);
                end;
                newur(tmpind2).latency = urlatency;
                newur(tmpind2).length  = lengths(tmpindex);
                rmnested = [ rmnested findnestedur(newur, tmpind2) ];
                
                % update indices in event structure
                % ---------------------------------
                for index = 1:length(EEG.event)
                    if EEG.event(index).urevent >= tmpind2
                        EEG.event(index).urevent = EEG.event(index).urevent+1;
                    end;
                end;
            end; 
            
            % update event structure
            % ----------------------
            EEG.event(end+1).type  = 'boundary';
            EEG.event(end).latency = boundevents(tmpindex);
            if isfield(EEG, 'urevent') & isfield(EEG.urevent, 'latency')
                EEG.event(end).urevent = tmpind2;              
            end;
        end; 
	end;
    
    % sort urevents by latency
    % ------------------------
    if isfield(EEG, 'urevent') & isfield(EEG.urevent, 'latency')
        [EEG.event EEG.urevent] = removenested(EEG.event, newur, rmnested);
    end;
    
	EEG = pop_editeventvals( EEG, 'sort', { 'latency', [0] } );
    if checkevent
        EEG = eeg_checkset(EEG, 'eventconsistency');
    end;

% look for nested ur events
% retrun indices of nested events and
% their total length
% -----------------------------------
function [ indnested, addlen ] = findnestedur(newur, ind);
    indnested = [];
    addlen = 0;
    tmpind = ind+1;
    
    while tmpind < length(newur) & ...
        newur(tmpind).latency < newur(ind).latency+newur(ind).length
        if strcmpi(newur(tmpind).type, 'boundary')
            indnested = [ indnested tmpind ];
            addlen    = addlen + newur(tmpind).length;
        end;
        tmpind = tmpind+1;
    end;
    
% remove urevent and recompute indices
% ------------------------------------
function [event, urevent] = removenested(event, urevent, nestind);
    
    if length(nestind) > 0
        fprintf('Debug msg: removing %d nested urevents\n', length(nestind));
        nestind = sort(nestind);
        urind = cell2mat({ event.urevent }); % this must not be done in the loop
                                             % since the indices are dyanmically updated
    end;
    
    for ind = 1:length(nestind)
        % find event urindices higher than the urevent to suppress
        % --------------------------------------------------------
        tmpind = find( urind > nestind(ind) );
        for indevent = tmpind
            event(indevent).urevent = event(indevent).urevent-1;
        end;    
    end;
    
    urevent(nestind) = [];
    