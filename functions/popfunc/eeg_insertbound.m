% eeg_insertbound() - insert boundary event in event structure.
%
% Usage:
%   >> EEGOUT = eeg_insertbound( EEG, latency, regions );
%
% Inputs:
%   EEG        - EEG dataset
%   latency    - relative latency (in the event structure) of boundary 
%                event(s)
%   regions    - removed regions in original dataset
%
% Outputs:
%   EEG        - output EEG dataset with 'event' and 'urevent' fields
%                updated
%
% Notes:
%   This function performs the following: add boundary events to the 
%   'event' and 'urevent' structures. Update urevent indices in the
%   'event' structure. Add a length field to the 'urevent' boundaries.
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
% Revision 1.1  2004/04/20 01:11:39  arno
% Initial revision
%

function EEG = eeg_insertbound( EEG, boundevents, regions );
    
    if nargin < 2
        help eeg_insertbound;
        return;
    end;
    
    if isfield(EEG, 'urevent')
        newur = EEG.urevent;
    end;
	for tmpindex = 1:length(boundevents)
        if boundevents(tmpindex) >= 1 & boundevents(tmpindex) <= EEG.pnts
            
            % insert event at the correct location in the urevent structure
            % -------------------------------------------------------------
            if isfield(EEG, 'urevent') & isfield(EEG.urevent, 'latency')
                % find event succeding boundary
                % ------------------------------
                urlatency = eeg_urlatency(EEG.urevent, regions(tmpindex,1)-0.5);
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
                newur(tmpind2).length  = regions(tmpindex,2)-regions(tmpindex,1)+1;
                
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
        EEG.urevent = newur;
    end;
    
	EEG = pop_editeventvals( EEG, 'sort', { 'latency', [0] } );
	EEG = eeg_checkset(EEG, 'eventconsistency');

