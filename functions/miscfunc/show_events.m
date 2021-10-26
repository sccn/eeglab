% show_events() - Display events in epochs. Events selected by 
%                 make_timewarp() function can be optionally highlighted.
%                 Each epoch is visualized as a row in the output image with
%                 events marked by colored rectangles.
%
% Usage:
%       >> im = show_events(EEG, 'key1',value1,'key2',value2,...);
%
% Inputs:
%
%   EEG                 - dataset structure.
%
%
% Optional inputs (in 'key', value format):
%
%   'eventThicknessCoef'  - (positive number) adjust the thickness of event 
%                           markers in the image. The default value is 1. 
%                           Lower values reduce thickness and higher values 
%                           increase it. For example, to decrease event
%                           rectangle width by half, set this parameter to
%                           0.5. {default: 1}
%   'eventNames'          - a cell array containing names of events to be
%                           displayed. Order and repetition are ignored. 
%                           {default: all events if timeWarp is not provided)
%
%   'timeWarp'            - a structure with latencies (time-warp matrix
%                           for newtimef function) and epochs with correct
%                           sequence created by make_timewarp() function.
%                           the subsets {default: false|0)
% Outputs:
%
%   im                    - color (RGB) image showing events of interest in
%                           epochs.
%
% Example:
%
%   % To display all events in data structure EEG
%   >> show_events(EEG);
%
%   % To highlight events selected by make_timewarp() function with thin
%   %  event markers
%  >> show_events(EEG, 'eventThicknessCoef', 0.5, 'eventNames', timeWarp.eventSequence,
%   'timeWarp', timeWarp);
%
% Version 1.1
% Author: Nima Bigdely Shamlo, SCCN/INC/UCSD, 2008
% See also: make_timewarp(), newtimef()

% Copyright (C) Nima Bigdely Shamlo, SCCN/INC/UCSD, 2008
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

function  im = show_events(EEG, varargin)


%% check inputs
EEG = change_events_to_string(EEG);

inputKeyValues = finputcheck(varargin, ...
   {'eventThicknessCoef'     'real'    [0 Inf]          1; ...
    'eventNames'     'cell'    {}        {}  ; ...
    'timeWarp'        'struct'    []          struct([]); ...
    });

eventThicknessCoef = 0;
eventNames = {};
timeWarp=[];

% place key values into function workspace variables
inputKeyValuesFields = fieldnames(inputKeyValues);
for i=1:length(inputKeyValuesFields)
    eval([inputKeyValuesFields{i} '= inputKeyValues.' inputKeyValuesFields{i} ';']);
end

if length(eventNames) == 0
    if length(timeWarp) == 0
        eventNames = [];
        for i=1:length(EEG.epoch)
            eventNames =  [eventNames EEG.epoch(i).eventtype];
        end
    else
        for i=1:length(timeWarp.eventSequence)
            if ischar(timeWarp.eventSequence{i})
                eventNames = [eventNames timeWarp.eventSequence{i}];
            else % in case it is a cell of strings
                for j=1:length(timeWarp.eventSequence{i})
                    eventNames = [eventNames timeWarp.eventSequence{i}{j}];
                end
            end
        end
        %        eventNames = timeWarp.eventSequence; % if event names is not provided only show timeWarp events
    end
end

if length(timeWarp) == 0
    timeWarp=[];
    timeWarp.latencies = [];
    timeWarp.epochs = 1:length(EEG.epoch);
end
%% set image parameters

imWidth = 300*4;
imHeight = 240*4;
im = double(zeros(imHeight, imWidth, 3));

%% find latencies for events of interest in each epoch

uniqueEventNames = unique_cell_string(eventNames);
for epochNumber = 1:length(EEG.epoch)
    for eventNumber = 1:length(uniqueEventNames)
       ids = find(strcmp(EEG.epoch(epochNumber).eventtype,uniqueEventNames{eventNumber})); % index for events with this type in the current epoch
       latency{epochNumber, eventNumber} = cell2mat(EEG.epoch(epochNumber).eventlatency(ids));
    end
end

%% find a good default value for event marker width based on quantiles
% of inter-event latencies.

intervals = [];  
for i=1:size(latency,1) 
    intervals = [intervals diff(sort(cell2mat(latency(i,:)),'ascend'))]; %#ok<AGROW>
end
 
if isempty(intervals)
    eventLineWidth = 0.05*eventThicknessCoef;
else
    eventLineWidth = round(eventThicknessCoef * 0.5 * imWidth * -quantile(-intervals,0.8) /  ((EEG.xmax - EEG.xmin)*1000));
end


%% place colored rectangle block o image to display events
eventColors = lines(length(uniqueEventNames));
for epochNumber = 1:length(EEG.epoch)
    for eventNumber = 1:length(uniqueEventNames)
        for eventInstanceNumber = 1:length(latency{epochNumber, eventNumber})
            startHeight =  max(1,round((epochNumber - 1) * imHeight/ length(EEG.epoch)));
            endHeight = min(imHeight, round((epochNumber) * imHeight/ length(EEG.epoch)));
            startWidth = max(1,round(-eventLineWidth/2 + imWidth * (latency{epochNumber, eventNumber}(eventInstanceNumber) - (EEG.xmin*1000))/ ((EEG.xmax - EEG.xmin)*1000)));
            endWidth =  min(imWidth,startWidth + eventLineWidth);
            colorForEvent =  eventColorBasedOnTypeAndAcceptance(eventColors, eventNumber, epochNumber, latency{epochNumber, eventNumber}(eventInstanceNumber), timeWarp);
            for c = 1:3 % to make sure the highlighted events are not overwritten by other events in the image, we use 'max' and give the brightest color priority               
                im(startHeight:endHeight,startWidth:endWidth,c) = max(im(startHeight:endHeight,startWidth:endWidth,c), colorForEvent(c));
            end
        end
    end
end

%% plot image with legend and axis labels

figure;
hold on;
for i=1:length(uniqueEventNames)
    plot([0 0],rand * [0 0],'Color',eventColors(i,:),'linewidth',10);
    uniqueEventNames{i} = strrep(uniqueEventNames{i},'_','-'); % replace _ with - so it is displayed correctly in the legend
end

legend(uniqueEventNames, 'Location', 'NorthWest');


image(round([EEG.xmin*1000 EEG.xmax*1000 ]), [1 length(EEG.epoch)],im);
xlim(round([EEG.xmin*1000 EEG.xmax*1000 ]))
ylim([1 length(EEG.epoch)])
axis normal;
xlabel('Latency', 'fontsize',16);
xlabel('Latency', 'fontsize',16);
ylabel('Epochs', 'fontsize',16);

if nargout == 0  % suppress output if not requested
    clear im ;
end
    


function color = eventColorBasedOnTypeAndAcceptance(eventColors, eventNumber, epochNumber, eventLatency,timeWarp)
accepted = ismember_bc(epochNumber, timeWarp.epochs);
if ~isempty(timeWarp.latencies)
    matchedEpoch = find(timeWarp.epochs ==  epochNumber);
    accepted = accepted && ~isempty(find(timeWarp.latencies(matchedEpoch,:) == eventLatency, 1));
end

color = eventColors(eventNumber,:);

if ~accepted
    color = color*0.3; % dim the color of unaccepted events.
end

function EEG = change_events_to_string(EEG)
needChange = false; 
for i=1:length(EEG.event)
    if ~ischar(EEG.event(i).type)
        EEG.event(i).type = num2str( EEG.event(i).type );
        needChange = true;
    end
end

if needChange
    for e=1:length(EEG.epoch)
        for i=1:length(EEG.epoch(e).eventtype)
            EEG.epoch(e).eventtype(i) = {num2str(cell2mat(EEG.epoch(e).eventtype(i)))};
        end
    end
end
