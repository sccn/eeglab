% eegrej() - reject/excise arbitrary periods from continuous EEG data 
%            (e.g., EEG.data).
%
% Usage:
%   >> [outdata newt newevents boundevents] = ...
%            eegrej( indata, regions, timelength, eventlatencies);
%
% Inputs:
%   indata     - input data (channels, frames). If indata is a string, 
%                the function use the disk to perform the rejection
%   regions    - array of regions to suppress. [beg end] x number of 
%                regions. 'beg' and 'end' are expressed in term of points
%                in the input dataset. The size() of the array should be
%                (2, number of regions).
%   timelength - length in time (s) of the input data. Only used to compute 
%                new total data length after rejections (newt).
%   eventlatencies - vector of event latencies in data points. 
%                    Default []=none.
%
% Outputs:
%   outdata    - output dataset
%   newt       - new total data length 
%   newevents  - new event latencies. If the event was in a removed
%                region, NaN is returned.
%   boundevents - boundary events latencies 
%
% Exemple: 
%   [outdat t] = eegrej( 'EEG.data', [1 100; 200 300]', [0 10]);
%   this command pick up two regions in EEG.data (from point 1 to
%   point 100, and point 200 to point 300) and put the result in 
%   the EEG.data variable in the global workspace (thus allowing
%   to perform the operation even if EEG.data takes all the memory)
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [indata, times, events, boundevents] = eegrej( indata, regions, times, events )

if nargin < 2
   help eegrej;
	return;
end;

if nargin < 4
    events = [];
end;

if isstr(indata)
  datlen = evalin('base', [ 'size(' indata ',2)' ]);
else
  datlen = size(indata, 2);
end;

reject = zeros(1,datlen);

% Checking for extreme values in regions (correcting round) RMC
regions = round(regions);
regions(regions > size(indata, 2)) = size(indata, 2);
regions(regions < 1) = 1;
regions = sortrows(sort(regions,2));        % Sorting regions %regions = sort(regions,1); RMC
for i=2:size(regions,1)
    if regions(i-1,2) >= regions(i,1)
        regions(i,1) = regions(i-1,2)+1;
    end;
end;

for i=1:size(regions,1)
    reject(regions(i,1):regions(i,2)) = 1;
end;

% recompute event times
% ---------------------
rmEvent = [];
rejectedEvents = {};
oriEvents = events;
if ~isempty(events)
    eventLatencies = [ events.latency ];
    oriEventLatencies = eventLatencies;
    for i=1:size(regions,1)
        
        % indices of removed events
        rejectedEvents{i} = find( oriEventLatencies > regions(i,1) & oriEventLatencies < regions(i,2) );
        rmEvent = [ rmEvent rejectedEvents{i} ];
        
        % remove events within the selected regions
        eventLatencies( oriEventLatencies > regions(i,1) ) = eventLatencies( oriEventLatencies > regions(i,1) )-(regions(i,2)-regions(i,1)+1);
        %if i == 24, dsafds; end;
        
    end;
    for iEvent = 1:length(events)
        events(iEvent).latency = eventLatencies(iEvent);
    end;
        
    events(rmEvent) = [];
end;

% generate boundaries latencies
% -----------------------------
boundevents = regions(:,1)-1;
for iRegion1=1:size(regions,1)
    duration(iRegion1)    = regions(iRegion1,2)-regions(iRegion1,1)+1;
    
    % add nested boundary events
    if ~isempty(events) && isstr(events(1).type) && isfield(events, 'duration')
        selectedEvent = oriEvents(rejectedEvents{iRegion1});
        indBound      = strmatch('boundary', { selectedEvent.type });
        duration(iRegion1) = duration(iRegion1) + sum([selectedEvent(indBound).duration]);
    end;
    
	for iRegion2=iRegion1+1:size(regions)
		boundevents(iRegion2) = boundevents(iRegion2) - (regions(iRegion1,2)-regions(iRegion1,1)+1);
    end;
end;
boundevents = boundevents+0.5;
boundevents(boundevents < 0) = [];

% reject data
% -----------
if isstr(indata)
  disp('Using disk to reject data');
  increment = 10000;
  global elecIndices;
  evalin('base', 'global elecIndices');
  elecIndices = find(reject == 0);
  evalin('base', 'fid = fopen(''tmpeegrej.fdt'', ''w'')');
  evalin('base', ['numberrow = size(' indata ',1)']);
  evalin('base', ['for indextmp=1:10000:length(elecIndices),', ...
		  '   endtmp = min(indextmp+9999, length(elecIndices));' ...
		  '   fwrite(fid,' indata '(:,elecIndices(indextmp:endtmp)), ''float'');'...
		  'end']);
  evalin('base', 'fclose(fid)');
  evalin('base', 'clear global elecIndices');  
  evalin('base', [ indata '=[]; clear ' indata '']);  
  evalin('base', 'fid = fopen(''tmpeegrej.fdt'', ''r'')');
  evalin('base', [ indata '= fread(fid, [numberrow ' int2str(datlen) '], ''float'')']);
  evalin('base', 'fclose(fid)');
  evalin('base', 'clear numberrow indextmp endtmp fid');  
  evalin('base', 'delete(''tmpeegrej.fdt'')');  
else
  indata(:,reject == 1) = [];
end;
times = times * size(indata,2) / length(reject);

% merge boundary events
% ---------------------
for iBound = length(boundevents):-1:2
    if boundevents(iBound) == boundevents(iBound-1)
        duration(iBound-1) = duration(iBound-1)+duration(iBound);
        boundevents(iBound) = [];
        duration(iBound) = [];
    end;
end;
if ~isempty(boundevents) && boundevents(end) > size(indata,2)
    boundevents(end) = [];
end;

% insert boundary events
% ----------------------
for iRegion1=1:length(boundevents)
    if boundevents(iRegion1) > 1 && boundevents(iRegion1) < size(indata,2)
        events(end+1).type = 'boundary';
        events(end).latency  = boundevents(iRegion1);
        events(end).duration = duration(iRegion1);
    end;
end;
if ~isempty(events) && isfield(events, 'latency')
    events([ events.latency ] < 1) = [];
    alllatencies = [ events.latency ];
    [tmp, sortind] = sort(alllatencies);
    events = events(sortind);
end;

return;
