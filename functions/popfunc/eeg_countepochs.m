% eeg_countepochs() Count how many epochs there are of each type
%
% Usage:
%  >> eeg_countepochs(EEG);
%
% Inputs:
%   EEG            - input dataset 
%   epochmarker    - ['type'|'eventtype'] indicates which part of the
%   EEG.epoch structure the different trial types are stored in. Depending
%   on what system the data are from and how they were preprocessed, this
%   may be either EEG.epoch.type or EEG.epoch.eventtype. Defaults to
%   'type'.

%  Outputs:
%    sweeps          - Scalar structure showing, for each epoch type
%    (sweeps.types) the number of sweeps in the dataset with that type
%    (sweeps.counts)
%
% Example:
% eeg_countepochs( EEG, 'type' )
% eeg_countepochs( EEG, 'eventtype' )
%
%  Author: Stephen Politzer-Ahles, University of Kansas, 2013

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

function [sweeps] = eeg_countepochs(EEG, epochmarker);

if nargin < 1
	help eeg_countepochs;
	return;
end;

% Initialize an array which will keep counts
clearvars types counts
types{1} = '';
counts = [0];

% Iterate through all trials
for trial=1:length(EEG.epoch)
	
    % Default 'epochmarker' to 'type' if no input was provided.
    if nargin < 2
        epochmarker = 'type';
    end;
    
    % Look for epochs that have >1 event and find which event is epoched around
    eventindex = 1; % if only 1 event, eventindex=1
    if length(EEG.epoch(trial).eventlatency) > 1
        for eventnum = 1:length(EEG.epoch(trial).eventlatency)
            if EEG.epoch(trial).eventlatency{eventnum}==0
                eventindex = eventnum;
            end % end if
        end % end for
    end % end if
    
	% Get the trigger number for this trial
    if strcmp(epochmarker,'eventtype')                      
        
        % change event type from cell/number event to string
        if iscell(EEG.epoch(trial).eventtype(eventindex))
            type = EEG.epoch(trial).eventtype{eventindex};
        elseif isnumeric(EEG.epoch(trial).eventtype(eventindex))
            type = num2str( EEG.epoch(trial).eventtype(eventindex) );
        else % is char
            type = EEG.epoch(trial).eventtype;
        end
    elseif strcmp(epochmarker,'type')
                      
        % change cell/number event to string
        if iscell(EEG.epoch(trial).type(eventindex))
            type = EEG.epoch(trial).type{eventindex};
        elseif isnumeric(EEG.epoch(trial).type(eventindex))
            type = num2str( EEG.epoch(trial).type(eventindex) );
        else % is char
            type = EEG.epoch(trial).type;
        end
    end; % end if

	% Find the row of the array that corresponds to this trigger (or if this trigger has not been counted yet)
	index_of_type = find( cellfun(@(x) strcmp(x,type), types) );

	% If this trial has a trigger that has already been counted before,
	% 'row_of_type' will be a number. If not, it will be an empty
	% array.
	if index_of_type
		
		% If we already have a count for this trigger, increment that count by 1
		counts(index_of_type) = counts(index_of_type) + 1;
	else

		% If not, create a new count for this type
		types = [types type];
        counts(end+1) = 1;
	end; % end if

end; % end for

% Remove the unnecessary first item
types = types(2:end);
counts = counts(2:end);

% Sorts
[Y,I] = sort( types );

% Add types and sweeps to a scalar structure
sweeps.types = types(I);
sweeps.counts = counts(I);

end
