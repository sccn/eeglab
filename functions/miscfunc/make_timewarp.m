% make_timewarp() - Select a subset of epochs containing a given event sequence, and return 
%                   a matrix of latencies for time warping the selected epochs to a common 
%                   timebase in newtimef(). Events in the given sequence may be further 
%                   restricted to those with specified event field values.
% Usage:
%       >> timeWarp = make_timewarp(EEG, eventSequence, 'key1',value1,
%       'key2',value2,...);
%
% Inputs:
%   EEG                 - dataset structure
%   eventSequence       - cell array containing a sequence of event type strings. 
%                         For example, to select epochs containing a 'movement Onset' 
%                         event followed by a 'movement peak', use
%                            {'movement Onset' 'movement peak'} 
%                         The output timeWarp matrix will contain the epoch latencies
%                         of the two events for each selected epoch. 
%
% Optional inputs (in 'key', value format):
%   'baselineLatency'     - (ms) the minimum acceptable epoch latency for the first
%                         event in the sequence {default: 0}
%   'eventConditions'     - cell array specifying conditions on event fields.
%                         For example, for a sequence consisting of two
%                         events with (velocity) fields vx and vy, use
%                            {vx>0' 'vx>0 && vy<1'} 
%                         To accept events unconditionally, use empty strings,
%                         for example {'vx>0','',''} {default: no conditions}
%  'maxSTDForAbsolute'    - (positive number of std. devs.) Remove epochs containing events 
%                         whose latencies are more than this number of standard deviations
%                         from the mean in the selected epochs. {default: Inf -> no removal}
%  'maxSTDForRelative'    - (positive number of std. devs.) Remove epochs containing inter-event
%                         latency differences larger or smaller than this number of standard deviations
%                         from the mean in the selected epochs. {default: Inf -> no removal}
% Outputs:
%   timeWarp            - a structure with latencies (time-warp matrix with fields
%                          timeWarp.latencies - an (N, M) timewarp matrix for use in newtimef() 
%                           where    N = number of selected epochs containing the specified sequence,
%                                    M = number of events in specified event sequence.
%                          timeWarp.epochs - a (1, M) vector giving the index of epochs with the 
%                              specified sequence. Only these epochs should be passed to newtimef().
%                          timeWarp.eventSequence - same as the 'eventSequence' input variable.
% Example:
%   % Create a timewarp matrix for a sequence of events, first an event of type 'movement Onset' followed by 
%   % a 'movement peak' event, with event fields vx and vy:
%
%   >> timeWarp = make_timewarp(EEG, {'movement Onset' 'movement%   peak'},'baselineLatency', 0 ...
%                  ,'eventConditions', {'vx>0' 'vx>0 && vy<1'});
%
%   % To remove events with latencies more than 3 standard deviations from the mean OR 2 std. devs. 
%   % from the mean inter-event difference:
%
%   >> timeWarp = make_timewarp(EEG, {'movement Onset' 'movement peak'},'baselineLatency', 0, ... 
%                  'eventConditions', {'vx>0' 'vx>0 && vy<1'},'maxSTDForAbsolute', 3,'maxSTDForRelative', 2);
%
% Author: Nima Bigdely Shamlo, SCCN/INC/UCSD, 2008
% See also: show_events(), newtimef()

function timeWarpStructure = make_timewarp(EEG, eventSequence, varargin)

inputKeyValues = finputcheck(varargin, ...
   {'baselineLatency'     'real'    []          0; ...
    'eventConditions'     'cell'    {}        {}  ; ...
    'maxSTDForAbsolute'        'real'    [0 Inf]          Inf; ...
    'maxSTDForRelative'        'real'    [0 Inf]         Inf...
    });

baselineLatency = [];
maxSTDForAbsolute = 0;
maxSTDForRelative = 0;
eventConditions = [];

% place key values into function workspace variables
inputKeyValuesFields = fieldnames(inputKeyValues);
for i=1:length(inputKeyValuesFields)
    eval([inputKeyValuesFields{i} '= inputKeyValues.' inputKeyValuesFields{i} ';']);
end;


if length(eventConditions) < length(eventSequence)
    for i = (length(eventConditions)+1):length(eventSequence)
        eventConditions{i} = '';
    end;
end;

epochsIsAcceptable = ones(1, length(EEG.epoch));
for epochNumber = 1:length(EEG.epoch)
    eventNameID = 1;
    minimumLatency = baselineLatency;
    timeWarp{epochNumber} = [];
    while eventNameID <= length(eventSequence) % go tthrought event names and find the event that comes after a certain event with correct type and higher latency

        firstLatency = eventsOfCertainTypeAfterCertainLatencyInEpoch(EEG.epoch(epochNumber), eventSequence{eventNameID}, minimumLatency, eventConditions{eventNameID});
        if isempty(firstLatency) % means there were no events, so the epoch is not acceptable
            break;
        else
            timeWarp{epochNumber} = [timeWarp{epochNumber}; firstLatency];
            minimumLatency = firstLatency;
            eventNameID = eventNameID + 1;
        end;
    end;

    if length(timeWarp{epochNumber}) < length(eventSequence)
        epochsIsAcceptable(epochNumber) = false;
    end;
end;


acceptableEpochs = find(epochsIsAcceptable);

if isempty(acceptableEpochs)
    timeWarp = {}; % no epoch meet the criteria
else
    timeWarp = cell2mat(timeWarp(acceptableEpochs));
   
    rejectedEpochesBasedOnLateny = union_bc(rejectEventsBasedOnAbsoluteLatency(timeWarp), rejectEventsBasedOnRelativeLatency(timeWarp));
    timeWarp(:,rejectedEpochesBasedOnLateny) = [];
    acceptableEpochs(rejectedEpochesBasedOnLateny) = [];
    
end;

% since latencies and accepted epochs always have to be together, we put them in one structure
timeWarpStructure.latencies = timeWarp'; % make it suitable for newtimef()

if isempty(timeWarpStructure.latencies) % when empty, it becomes a {} instead of [], so we change it to []
    timeWarpStructure.latencies = [];
end;

timeWarpStructure.epochs = acceptableEpochs;
timeWarpStructure.eventSequence = eventSequence;

    function rejectedBasedOnLateny =  rejectEventsBasedOnRelativeLatency(timeWarp)
        % remeve epochs in which the time warped event is further than n
        % standard deviations to mean of latency distance between events
        timeWarpDiff =  diff(timeWarp);
        rejectedBasedOnLateny = [];
        for eventNumber = 1:size(timeWarpDiff, 1)
            rejectedBasedOnLateny = [rejectedBasedOnLateny find(abs(timeWarpDiff(eventNumber, :)-  mean(timeWarpDiff(eventNumber, :))) > maxSTDForRelative * std(timeWarpDiff(eventNumber, :)))];
        end;
        rejectedEpochesBasedOnLateny = unique_bc(rejectedBasedOnLateny);
    end

    function rejectedBasedOnLateny =  rejectEventsBasedOnAbsoluteLatency(timeWarp)
        % remeve instances in which the time warped event is further than n
        % standard deviations to mean
        rejectedBasedOnLateny = [];
        for eventNumber = 1:size(timeWarp, 1)
            rejectedBasedOnLateny = [rejectedBasedOnLateny find(abs(timeWarp(eventNumber, :)-  mean(timeWarp(eventNumber, :))) >  maxSTDForAbsolute* std(timeWarp(eventNumber, :)))];
        end;
        rejectedEpochesBasedOnLateny = unique_bc(rejectedBasedOnLateny);
    end

    function [firstLatency resultEventNumbers] = eventsOfCertainTypeAfterCertainLatencyInEpoch(epoch, certainEventType, certainLatency, certainCondition)
        resultEventNumbers = [];
        firstLatency = []; % first event latency that meets the critria
        for eventNumber = 1:length(epoch.eventtype)
%            if strcmp(certainEventType, epoch.eventtype(eventNumber)) && epoch.eventlatency{eventNumber} >= certainLatency && eventMeetsCondition(epoch, eventNumber, certainCondition)
             if  eventIsOfType(epoch.eventtype(eventNumber), certainEventType)  && epoch.eventlatency{eventNumber} >= certainLatency && eventMeetsCondition(epoch, eventNumber, certainCondition)
                resultEventNumbers = [resultEventNumbers eventNumber];

                if isempty(firstLatency)
                    firstLatency = epoch.eventlatency{eventNumber};
                end;
            end;
        end;
    end

    function result = eventMeetsCondition(epoch, eventNumber, condition)
        if strcmp(condition,'') || strcmp(condition,'true')
            result = true;
        else

            % get the name thatis before themfor field in the epoch, then remove 'event' name
            epochField = fieldnames(epoch);
            for i=1:length(epochField)
                epochField{i} = strrep(epochField{i},'event',''); % remove event from the beginning of field names
                condition = strrep(condition, epochField{i}, ['cell2mat(epoch.event' epochField{i} '(' num2str(eventNumber) '))' ]);
            end;

            result = eval(condition);
        end;

    end

    function result = eventIsOfType(eventStr, types)
        % if events are numbers, turn them into strings before comparison
        if ischar(types)
            if iscell(eventStr) && isnumeric(cell2mat(eventStr))
                eventStr = num2str(cell2mat(eventStr));
            end;
            
            result = strcmp(eventStr, types);
        else % it must be a cell of strs
            result = false;
            for i=1:length(types)
                result = result ||  strcmp(eventStr, types{i});                
            end;
        end;
    end
        
end
