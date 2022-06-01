% eeg_eegrej() - reject porition of continuous data in an EEGLAB 
%                dataset
%
% Usage:
%   >> EEGOUT = eeg_eegrej( EEGIN, regions );
%
% Inputs:
%   INEEG      - input dataset
%   regions    - array of regions to suppress. number x [beg end]  of 
%                regions. 'beg' and 'end' are expressed in term of points
%                in the input dataset. Size of the array is
%                number x 2 of regions.
%
% Outputs:
%   INEEG      - output dataset with updated data, events latencies and 
%                additional boundary events.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 8 August 2002
%
% See also: eeglab(), eegplot(), pop_rejepoch()

% Copyright (C) 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [EEG, com] = eeg_eegrej( EEG, regions);

com = '';
if nargin < 2
	help eeg_eegrej;
	return;
end
if nargin<3
    probadded = [];
end
if isempty(regions)
	return;
end

% regions = sortrows(regions,3); % Arno and Ramon on 5/13/2014 for bug 1605

% Ramon on 5/29/2014 for bug 1619
if size(regions,2) > 2
    regions = sortrows(regions,3);
else
    regions = sortrows(regions,1);
end

try
    % For AMICA probabilities...Temporarily add model probabilities as channels
    %-----------------------------------------------------
    if isfield(EEG.etc, 'amica') && ~isempty(EEG.etc.amica) && isfield(EEG.etc.amica, 'v_smooth') && ~isempty(EEG.etc.amica.v_smooth) && ~isfield(EEG.etc.amica,'prob_added')
        if isfield(EEG.etc.amica, 'num_models') && ~isempty(EEG.etc.amica.num_models)
            if size(EEG.data,2) == size(EEG.etc.amica.v_smooth,2) && size(EEG.data,3) == size(EEG.etc.amica.v_smooth,3) && size(EEG.etc.amica.v_smooth,1) == EEG.etc.amica.num_models

                EEG = eeg_formatamica(EEG);
                %-------------------------------------------

                [EEG com] = eeg_eegrej(EEG,regions);

                %-------------------------------------------

                EEG = eeg_reformatamica(EEG);
                EEG = eeg_checkamica(EEG);
                return;
            else
                disp('AMICA probabilities not compatible with size of data, probabilities cannot be epoched')
                disp('Load AMICA components before extracting epochs')
                disp('Resuming rejection...')
            end
        end

    end
    % ------------------------------------------------------
catch
    warnmsg = strcat('your dataset contains amica information, but the amica plugin is not installed.  Continuing and ignoring amica information.');
    warning(warnmsg)
end

% handle regions from eegplot
% ---------------------------
if size(regions,2) > 2, regions = regions(:, 3:4); end
regions = combineregions(regions);

% remove events within regions
% ----------------------------
if ~isempty(EEG.event)
    allEventLatencies = [ EEG.event.latency];
    allEventFlag      = zeros(1,length(allEventLatencies));
    for iRegion = 1:size(regions,1)
        allEventFlag = allEventFlag | ( allEventLatencies >= regions(iRegion,1) & allEventLatencies <= regions(iRegion,2));
    end
    EEG.event(allEventFlag) = [];
end

% reject data
% -----------
[EEG.data, EEG.xmax, event2, boundevents] = eegrej( EEG.data, regions, EEG.xmax-EEG.xmin, EEG.event);
oldEEGpnts = EEG.pnts;
oldEEGevents = EEG.event;
EEG.pnts   = size(EEG.data,2);
EEG.xmax   = EEG.xmax+EEG.xmin;
if length(event2) > 1 && event2(1).latency == 0, event2(1) = []; end
if length(event2) > 1 && event2(end).latency == EEG.pnts, event2(end) = []; end
if length(event2) > 2 && event2(end).latency == event2(end-1).latency
    if isfield(event2, 'type') && isequal(event2(end).type, event2(end-1).type)
        event2(end) = []; 
    end
end


% add boundary events
% -------------------
[ EEG.event ] = eeg_insertbound(EEG.event, oldEEGpnts, regions);
EEG = eeg_checkset(EEG, 'eventconsistency');
if ~isempty(EEG.event) && EEG.trials == 1 && EEG.event(end).latency > EEG.pnts
    EEG.event(end) = []; % remove last event if necessary
end

% double check event latencies
% the function that insert boundary events and recompute latency is
% delicate so we do it twice using different methods and check
% the results. It is longer, but accuracy is paramount.
eeglab_options;
if isfield(EEG.event, 'latency') && length(EEG.event) < 3000
    % assess difference between old and new event latencies
    [ eventtmp ] = eeg_insertboundold(oldEEGevents, oldEEGpnts, regions);
    if ~isempty(eventtmp)
        [~,indEvent] = sort([ eventtmp.latency ]);
        eventtmp = eventtmp(indEvent);
    end
    if ~isempty(eventtmp) && length(eventtmp) > length(EEG.event) && isfield(eventtmp, 'type') && eeg_isboundary(eventtmp(1))
        eventtmp(1) = [];
    end
    if isfield(eventtmp, 'duration')
        for iEvent=1:length(eventtmp)
            if isempty(eventtmp(iEvent).duration)
                eventtmp(iEvent).duration = 0;
            end
        end
    end
    if ~isempty(eventtmp) && eventtmp(end).latency > EEG.pnts
        eventtmp(end) = [];
    end
    % add initial event to eventtmp when missing
    if ~isempty(eventtmp) && ~isempty(EEG.event) && eeg_isboundary(eventtmp(1)) && ...
            EEG.event(1).latency == 0.5 && eventtmp(1).latency ~= 0.5
        if size(eventtmp,2) > 1
            eventtmp = [ eventtmp(1) eventtmp(1:end) ];
        else
            eventtmp = [ eventtmp(1) eventtmp(1:end)' ];
        end
        eventtmp(1).type = eeg_boundarytype(eventtmp);
        eventtmp(1).latency = 0.5;
        eventtmp(1).duration = EEG.event(1).duration;
    end
        
    differs = 0;
    for iEvent=1:min(length(EEG.event), length(eventtmp)-1)
        if ~issameevent(EEG.event(iEvent), eventtmp(iEvent)) && ~issameevent(EEG.event(iEvent), eventtmp(iEvent+1)) 
            differs = differs+1;
        end
    end
    
    if 100*differs/length(EEG.event) > 50
        db = dbstack;
        if length(db) > 1 
            fprintf(['bug 1971 warning for scritps older than 2017: see https://sccn.ucsd.edu/wiki/EEGLAB_bug1971\n' ]);
        end
    end
    
    alllats = [ EEG.event.latency ];
    if ~isempty(event2)
        otherlatencies = [event2.latency];
        if ~isequal(alllats, otherlatencies)
            warning([ 'Minor discrepency when recomputing event latency.' 10 '(this will not affect the accuracy of analyses).' 10 'Try to reproduce the problem and send us your dataset' ]);
        end
    end
end

% double check boundary event latencies
if ~isempty(EEG.event) && length(EEG.event) < 3000 && ischar(EEG.event(1).type) && isfield(EEG.event, 'duration') && isfield(event2, 'duration')
    try
        if ~isempty(EEG.event) && ischar(EEG.event(1).type)
            indBound1 = find(cellfun(@(x)strcmpi(num2str(x), 'boundary'), { EEG.event(:).type }));
            indBound2 = find(cellfun(@(x)strcmpi(num2str(x), 'boundary'), { event2(:).type }));
        else
            indBound1 = find([ EEG.event(:).type ] == -99);
            indBound2 = find([ event2(:).type ] == -99);
        end
        duration1 = [EEG.event(indBound1).duration]; duration1(isnan(duration1)) = [];
        duration2 = [event2(indBound2).duration]; duration2(isnan(duration2)) = [];
        if ~isequal(duration1, duration2)
            duration1(duration1 == 0) = [];
            if ~isequal(duration1, duration2)
                warning(['Inconsistency in boundary event duration.' 10 '(this will not affect the accuracy of analyses).' 10 'Try to reproduce the problem and send us your dataset' ]); 
            end
        end
    catch, warning('Unknown error when checking event latency - please send us your dataset');
    end
end

% debugging code below
% regions, n1 = 1525; n2 = 1545; n = n2-n1+1;
% a = zeros(1,n); a(:) = 1; a(strmatch('boundary', { event2(n1:n2).type })') = 8; 
% [[n1:n2]' alllats(n1:n2)' [event2(n1:n2).latency]' alllats(n1:n2)'-[event2(n1:n2).latency]' otherorilatencies(n1:n2)' a']
% figure; ev = 17; range = [-1000:1000]; plot(EEG.data(1,EEG.event(ev).latency+range)); hold on; plot(tmpdata(1,tmpevent(EEG.event(ev).urevent).latency+range+696), 'r'); grid on;

com = sprintf('EEG = eeg_eegrej( EEG, %s);', vararg2str({ regions })); 

% combine regions if necessary
% it should not be necessary but a 
% bug in eegplot makes that it sometimes is
% ----------------------------
% function newregions = combineregions(regions)
% newregions = regions;
% for index = size(regions,1):-1:2
%     if regions(index-1,2) >= regions(index,1)
%         disp('Warning: overlapping regions detected and fixed in eeg_eegrej');
%         newregions(index-1,:) = [regions(index-1,1) regions(index,2) ];
%         newregions(index,:)   = [];
%     end
% end

function res = issameevent(evt1, evt2)

res = true;
if isequal(evt1,evt2)
    return;
else
    if isfield(evt1, 'type') && isnumeric(evt2.type) && ~isnumeric(evt1.type) 
        evt2.type = num2str(evt2.type);
        if isequal(evt1,evt2)
            return;
        end
    end
    if isfield(evt1, 'duration') && isfield(evt2, 'duration')
        if isnan(evt1.duration) && isnan(evt2.duration)
            evt1.duration = 1;
            evt2.duration = 1;
        end
        if abs( evt1.duration - evt2.duration) < 1e-10
            evt1.duration = 1;
            evt2.duration = 1;
        end
        if isequal(evt1,evt2)
            return;
        end
    end
    if isfield(evt1, 'latency') && isfield(evt2, 'latency')
        if abs( evt1.latency - evt2.latency) < 1e-10
            evt1.latency = 1;
            evt2.latency = 1;
        end
        if isequal(evt1,evt2)
            return;
        end
    end
end
res = false;
return;

function newregions = combineregions(regions)
% 9/1/2014 RMC
regions = sortrows(sort(regions,2)); % Sorting regions
allreg = [ regions(:,1)' regions(:,2)'; ones(1,numel(regions(:,1))) -ones(1,numel(regions(:,2)')) ].';
allreg = sortrows(allreg,1); % Sort all start and stop points (column 1),

mboundary = cumsum(allreg(:,2)); % Rationale: regions will start always with 1 and close with 0, since starts=1 end=-1
indx = 0; count = 1;

while indx ~= length(allreg) 
    newregions(count,1) = allreg(indx+1,1);
    [tmp,I]= min(abs(mboundary(indx+1:end)));
    newregions(count,2) = allreg(I + indx,1);
    indx = indx + I ;
    count = count+1;
end

% Verbose
if size(regions,1) ~= size(newregions,1)
    disp('Warning: overlapping regions detected and fixed in eeg_eegrej');
end
