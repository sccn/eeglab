% eeg_addnewevents() Add new events to EEG structure. Both EEG.event and
% EEG.urevent are updated.
%
% Usage:
%  >> EEG = eeg_addnewevents(EEG, latencies, types, fieldNames, fieldValues);
%
% Inputs:
%   EEG            - input dataset
%   latencies      - cell containing numerical arrays for latencies of new
%                    events, each array corresponds to a different event type.
%   type           - cell array containing name of event types.
%
% Optional Inputs
%   fieldNames     - cell array containing names of fields to be added to event structure.
%   fieldValues    - A cell containing arrays for field values corresponding to fieldNames.
%                    Number of values for each field should be equal to the total number of 
%                    latencies (new events) added to dataset.
%  Outputs:
%    EEG          - EEG dataset with updated event and urevent fields
%
% Example:
% EEG = eeg_addnewevents(EEG, {[100 200] [300 400 500]}, {'type1' 'type2'}, {'field1' 'field2'}, {[1 2 3 4 5] [6 7 8 9]});
%
%  Author: Nima Bigdely Shamlo, SCCN/INC/UCSD, 2008
function EEG = eeg_addnewevents(EEG, eventLatencyArrays, types, fieldNames, fieldValues);



if ~isfield(EEG, 'event')
    EEG.event = [];
    EEG.urevent = [];

    EEG.event(1).type = 'dummy';
    EEG.event(1).latency = 1;
    EEG.event(1).duration = 0;
    EEG.event(1).urevent = 1;

    EEG.urevent(1).type = 'dummy';
    EEG.urevent(1).latency = 1;
    EEG.urevent(1).duration = 0;
end;

% add duration field if it does not exist
if length(EEG.event)>0 && ~isfield(EEG.event(1),'duration')
    EEG.event(1).duration = 0;
    EEG.urevent(1).duration = 0;
end;

if nargin<4
    fieldNames = [];
    fieldValues = [];
end;

newEventLatency = [];
for i=1:length(eventLatencyArrays)
    newEventLatency = [newEventLatency eventLatencyArrays{i}];
end;


newEventType = [];

for i=1:length(eventLatencyArrays{1})
    newEventType{i}  = types{1};
end;


for j=2:length(eventLatencyArrays)
    startIndex = length(newEventType);
    for i=1:length(eventLatencyArrays{j})
        newEventType{startIndex+i}  = types{j};
    end;
end;

% mix new and old events, sort them by latency and put them back in EEG
originalEventLatency = [];
originalEventType = [];
originalFieldNames = [];
for i=1:length(EEG.event)
    originalEventLatency(i) = EEG.event(i).latency;
    originalEventType{i} = EEG.event(i).type;
    originalEventFields(i) = EEG.event(i);
end;


% make sure that originalEventFields has all the new field names
if ~isempty(EEG.event)
    originalFieldNames = fields(originalEventFields);
    for f= 1:length(fieldNames)
        if ~isfield(originalEventFields, fieldNames{f})
            originalEventFields(length(originalEventFields)).(fieldNames{f}) = NaN;
        end;
    end;
end;

% make sure that newEventFields has all the original field names
for i=1:length(originalFieldNames)
    newEventFields(length(newEventLatency)).(originalFieldNames{i}) = NaN;
end;

for i=1:length(newEventLatency)
    newEventFields(i).latency = newEventLatency(i);
    newEventFields(i).type = newEventType{i};
    newEventFields(i).duration = 0;
    for f= 1:length(fieldNames)
        newEventFields(i).(fieldNames{f}) = fieldValues{f}(i);
    end;
end;

if ~isempty(EEG.event)
    %newEventFields = struct('latency', num2cell(newEventLatency), 'type', newEventType);
    combinedFields = [originalEventFields newEventFields];
    combinedLatencies = [originalEventLatency newEventLatency];
    combinedType = [originalEventType newEventType];
else
    combinedFields = newEventFields;
    combinedLatencies = newEventLatency;
    combinedType = newEventType;
end

[sortedEventLatency order] = sort(combinedLatencies,'ascend');
sortedEventType = combinedType(order);
combinedFields = combinedFields(order);
% put events in eeg

%EEG.urevent = [];
%EEG.event = [];
EEG = rmfield(EEG,'event');

for i=1:length(sortedEventLatency)
    %    EEG.urevent(i).latency = sortedEventLatency(i);
    %    EEG.urevent(i).type = sortedEventType{i};

    %    combinedFields(order(i)).urevent = i;

    EEG.event(i) = combinedFields(i);
    %    EEG.event(i).urevent = i;
end;

%% adding new urevents

originalUreventNumber = 1:length(EEG.urevent);
originalUreventLatency = zeros(1, length(EEG.urevent));
originalUreventFields= cell(1, length(EEG.urevent));

for i=1:length(EEG.urevent)
    originalUreventLatency(i) = EEG.urevent(i).latency;
    originalUreventFields{i} =  EEG.urevent(i);
end;

newUreventLatency = [];
newUreventType = [];

for i=1:length(EEG.event)
    if ~isfield(EEG.event,'urevent') || length(EEG.event(i).urevent) == 0 || isnan(EEG.event(i).urevent)
        %   newUreventLatency = [newUreventLatency newEventUrEventLatency(EEG, combinedFields, i)];

        % use eeg_urlatency to calculate the original latency based on
        % EEG.event duartions
        newUreventLatency = [newUreventLatency eeg_urlatency(EEG.event, EEG.event(i).latency)];
    else
        newUreventLatency = [newUreventLatency EEG.urevent(EEG.event(i).urevent).latency];
    end;
  
    newUreventFields{i} = EEG.event(i);
    
    newUreventEventNumber(i) = i;
end;

combinedEventNumber = newUreventEventNumber;%[NaN(1,length(EEG.urevent)) newUreventEventNumber];
combinedUrEventLatencies = newUreventLatency;%[originalUreventLatency newUreventLatency];
[sortedUrEventLatency order] = sort(combinedUrEventLatencies,'ascend');

% make urvent stucture ready
EEG.urevent = [];
EEG.urevent= newUreventFields{order(1)};

for i=1:length(order)
    %if ~isnan(newUreventEventNumber(i))
    
    EEG.urevent(i) = newUreventFields{order(i)};
    
    EEG.urevent(i).latency = combinedUrEventLatencies(order(i));
    EEG.event(newUreventEventNumber(i)).urevent = i;
    %end;
end;

if isfield(EEG.urevent,'urevent')
    EEG.urevent = rmfield(EEG.urevent,'urevent'); % remove urevent field
end;

% turn empty event durations into 0
for i=1:length(EEG.event)
    if isempty(EEG.event(i).duration)
        EEG.event(i).duration = 0;
    end;
end;

for i=1:length(EEG.urevent)
    if isempty(EEG.urevent(i).duration)
        EEG.urevent(i).duration = 0;
    end;
end;

%
% function latency = newEventUrEventLatency(EEG, combinedFields, i)
%
% %% looks for an event with urvent before the new event
% urlatencyBefore = [];
% currentEventNumber = i;
%
% while isempty(urlatencyBefore) && currentEventNumber > 1
%     currentEventNumber = currentEventNumber - 1;
%     if ~(~isfield(combinedFields(currentEventNumber),'urevent') || isempty(combinedFields(currentEventNumber).urevent) || isnan(combinedFields(currentEventNumber).urevent))
%         urlatencyBefore = EEG.urevent(combinedFields(currentEventNumber).urevent).latency;
%     end;
% end
%
% %% if no event with urevent is found before, look for an event with urvent after the new event
% if isempty(urlatencyBefore)
%     urlatencyAfter = [];
%     currentEventNumber = i;
%
%     while isempty(urlatencyAfter) && currentEventNumber < length(combinedFields)
%         currentEventNumber = currentEventNumber + 1;
%         if ~(~isfield(combinedFields(currentEventNumber),'urevent') ||  isempty(combinedFields(currentEventNumber).urevent) || isnan(combinedFields(currentEventNumber).urevent))
%             urlatencyAfter = EEG.urevent(combinedFields(currentEventNumber).urevent).latency;
%         end;
%     end
% end;
% %%
% if ~isempty(urlatencyBefore)
%     latency = urlatencyBefore + combinedFields(i).latency - combinedFields(currentEventNumber).latency;
% elseif ~isempty(urlatencyAfter)
%     latency = urlatencyAfter + combinedFields(currentEventNumber).latency - combinedFields(i).latency;
% else
%     latency = [];
% end;
