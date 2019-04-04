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
%
%  Outputs:
%    sweeps          - Scalar structure showing, for each epoch type
%    (sweeps.types) the number of sweeps in the dataset with that type
%    (sweeps.counts)
%
% Example:
% eeg_countepochs( EEG, 'type' )
% eeg_countepochs( EEG, 'eventtype' )
%
% Author: Stephen Politzer-Ahles, University of Kansas, 2013

% Copyright: Stephen Politzer-Ahles, University of Kansas, 2013
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

function [sweeps] = eeg_countepochs(EEG, epochmarker);

if nargin < 1
	help eeg_countepochs;
	return;
end

% Initialize an array which will keep counts
clearvars types counts
types{1} = '';
counts = [0];

% Iterate through all trials
for trial=1:length(EEG.epoch)
	
    % Default 'epochmarker' to 'type' if no input was provided.
    if nargin < 2
        epochmarker = 'type';
    end
    
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
