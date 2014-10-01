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

function [indata, times, newevents, boundevents] = eegrej( indata, regions, times, eventtimes );

if nargin < 2
   help eegrej;
	return;
end;

if isstr(indata)
  datlen = evalin('base', [ 'size(' indata ',2)' ]);
else
  datlen = size(indata, 2);
end;

reject = zeros(1,datlen);
regions = round(regions);
% Checking for extreme values in regions (correcting round) RMC
if max(regions(:)) > size(indata, 2)
    IndxOut = find(regions(:) > size(indata, 2));
    regions(IndxOut) = size(indata, 2);
end

regions = sortrows(sort(regions,2));        % Sorting regions %regions = sort(regions,1); RMC
Izero = find(regions == 0);                 % Find regions index == 0 to adjust them
if ~isempty(Izero), regions(Izero) = 1;end; % Fractional point below 1 adjusted to 1
for i=1:size(regions,1)
   try
      reject(regions(i,1):regions(i,2)) = 1;
   catch
      error(['Region ' int2str(i) ' out of bound']);
   end;
end;

% recompute event times
% ---------------------
newevents = [];
if exist('eventtimes') == 1 
    if ~isempty(eventtimes)
        try
            rmevent = find( reject(min(length(reject),max(1,round(eventtimes)))) == 1); % cko: sometimes, events may have latencies < 0.5 or >= length(reject)+0.5 (e.g. after resampling)
	    catch, error('Latency event out of bound'); end;
        eventtimes(rmevent) = NaN;
	    newevents = eventtimes;
	    for i=1:size(regions,1)
	        for index=1:length( eventtimes )
	            if ~isnan(eventtimes(index))
	                if regions(i,2) < eventtimes(index)
	                    newevents(index) = newevents(index) - (regions(i,2)-regions(i,1)+1);
	                end;
	            end;
	        end;
	    end;
    end;
end;                    

% generate boundaries latencies
% -----------------------------
boundevents = regions(:,1)-1;
for i=1:size(regions,1)
	for index=i+1:size(regions)
		boundevents(index) = boundevents(index) - (regions(i,2)-regions(i,1)+1);
    end;
end;
boundevents = boundevents+0.5;

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
  timeIndices = find(reject == 1);
  indata(:,timeIndices) = [];
end;
times = times * length(find(reject ==0)) / length(reject);

return;
