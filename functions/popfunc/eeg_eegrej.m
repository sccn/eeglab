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

function [EEG, com] = eeg_eegrej( EEG, regions);

com = '';
if nargin < 2
	help eeg_eegrej;
	return;
end;
if nargin<3
    probadded = [];
end
if isempty(regions)
	return;
end;

% regions = sortrows(regions,3); % Arno and Ramon on 5/13/2014 for bug 1605

% Ramon on 5/29/2014 for bug 1619
if size(regions,2) > 2
    regions = sortrows(regions,3);
else
    regions = sortrows(regions,1);
end;

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


if isfield(EEG.event, 'latency'),
   	 tmpevent = EEG.event;
     
     tmpdata = EEG.data; % REMOVE THIS, THIS IS FOR DEBUGGING %
     
     tmpalllatencies = [ tmpevent.latency ];

else tmpalllatencies = []; 
end;

% handle regions from eegplot
% ---------------------------
if size(regions,2) > 2, regions = regions(:, 3:4); end;
regions = combineregions(regions);

[EEG.data, EEG.xmax, event2, boundevents] = eegrej( EEG.data, regions, EEG.xmax-EEG.xmin, EEG.event);
oldEEGpnts = EEG.pnts;
oldEEGevents = EEG.event;
EEG.pnts   = size(EEG.data,2);
EEG.xmax   = EEG.xmax+EEG.xmin;

% add boundary events
% -------------------
[ EEG.event ] = eeg_insertbound(EEG.event, oldEEGpnts, regions);
EEG = eeg_checkset(EEG, 'eventconsistency');
if ~isempty(EEG.event) && EEG.trials == 1 && EEG.event(end).latency > EEG.pnts
    EEG.event(end) = []; % remove last event if necessary
end;

% double check event latencies
% the function that insert boundary events and recompute latency is
% delicate so we do it twice using different methods and check
% the results. It is longer, but accuracy is paramount.
if isfield(EEG.event, 'latency') && length(EEG.event) < 3000
    % assess difference between old and new event latencies
    [ eventtmp ] = eeg_insertboundold(oldEEGevents, oldEEGpnts, regions);
    if ~isempty(eventtmp) && length(eventtmp) > length(EEG.event) && isfield(eventtmp, 'type') && isequal(eventtmp(1).type, 'boundary')
        eventtmp(1) = [];
    end
    if isfield(eventtmp, 'duration')
        for iEvent=1:length(eventtmp)
            if isempty(eventtmp(iEvent).duration)
                eventtmp(iEvent).duration = 0;
            end
        end
    end
    differs = 0;
    for iEvent=1:min(length(EEG.event), length(eventtmp)-1)
        if ~issameevent(EEG.event(iEvent), eventtmp(iEvent)) && ~issameevent(EEG.event(iEvent), eventtmp(iEvent+1)) 
            differs = differs+1;
        end
    end
    if 100*differs/length(EEG.event) > 50
        fprintf(['BUG 1971 WARNING: IF YOU ARE USING A SCRIPT WITTEN FOR A PREVIOUS VERSION OF\n' ...
                'EEGLAB TO CALL THIS FUNCTION, BECAUSE YOU ARE REJECTING THE ONSET OF THE DATA,\n' ...
                'EVENTS WERE CORRUPTED. EVENT LATENCIES ARE NOW CORRECT (SEE https://sccn.ucsd.edu/wiki/EEGLAB_bug1971);\n' ]);
    end
    
    alllats = [ EEG.event.latency ];
    if ~isempty(event2)
        otherlatencies = [event2.latency];
        if ~isequal(alllats, otherlatencies)
            warning([ 'Discrepency when recomputing event latency.' 10 'Try to reproduce the problem and send us your dataset' ]);
        end
    end
end

% double check boundary event latencies
if ~isempty(EEG.event) && length(EEG.event) < 3000 && ischar(EEG.event(1).type) && isfield(EEG.event, 'duration') && isfield(event2, 'duration')
    try
        indBound1 = find(cellfun(@(x)strcmpi(num2str(x), 'boundary'), { EEG.event(:).type }));
        indBound2 = find(cellfun(@(x)strcmpi(num2str(x), 'boundary'), { event2(:).type }));
        duration1 = [EEG.event(indBound1).duration]; duration1(isnan(duration1)) = [];
        duration2 = [event2(indBound2).duration]; duration2(isnan(duration2)) = [];
        if ~isequal(duration1, duration2)
            warning(['Inconsistency in boundary event duration.' 10 'Try to reproduce the problem and send us your dataset' ]); 
        end;
    catch, warning('Unknown error when checking event latency - please send us your dataset');
    end;
end;

% debuging code below
% regions, n1 = 1525; n2 = 1545; n = n2-n1+1;
% a = zeros(1,n); a(:) = 1; a(strmatch('boundary', { event2(n1:n2).type })') = 8; 
% [[n1:n2]' alllats(n1:n2)' [event2(n1:n2).latency]' alllats(n1:n2)'-[event2(n1:n2).latency]' otherorilatencies(n1:n2)' a']
% figure; ev = 17; range = [-1000:1000]; plot(EEG.data(1,EEG.event(ev).latency+range)); hold on; plot(tmpdata(1,tmpevent(EEG.event(ev).urevent).latency+range+696), 'r'); grid on;

com = sprintf('%s = eeg_eegrej( %s, %s);', inputname(1), inputname(1), vararg2str({ regions })); 

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
%     end;
% end;

function res = issameevent(evt1, evt2)

res = true;
if isequal(evt1,evt2)
    return;
elseif isfield(evt1, 'duration') && isnan(evt1.duration) && isfield(evt2, 'duration') && isnan(evt2.duration)
    evt1.duration = 1;
    evt2.duration = 1;
    if isequal(evt1,evt2)
        return;
    end;
end;
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