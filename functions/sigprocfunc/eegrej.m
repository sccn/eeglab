% eegrej() - reject eeg regions in continuous eeg data.
%
% Usage:
%   >> [outdata newt newevents boundevents] = ...
%            eegrej( indata, regions, timelength, eventlatencies);
%
% Inputs:
%   indata     - input data nbchannel x nbpoints. If indata is a
%                string of character, the function use the disk to
%                perform the rejection
%   regions    - array of regions to suppress. [beg end] x number of 
%                regions. 'beg' and 'end' are expressed in term of points
%                in the input dataset. Size of the array is
%                2xnumber of regions.
%   timelength - time length. Only used to compute new total time length.
%   eventlatencies - vector of event latencies in data points. 
%                    Default []=none.
%
% Outputs:
%   outdata    - output dataset
%   newt       - new timelength
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

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.3  2002/07/31 16:11:15  arno
% debugging
%
% Revision 1.2  2002/05/21 20:51:29  scott
% removed ; from evalin() calls -sm
%
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%
% 01-25-02 reformated help & license -ad 
% 03-27-02 added event latency recalculation -ad 

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
	    rmevent = find( reject(round(eventtimes)) == 1);
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
  evalin('base', [ indata '= fread(fid, [numberrow ' int2str(nbpoint) '], ''float'')']);
  evalin('base', 'fclose(fid)');
  evalin('base', 'clear numberrow indextmp endtmp fid');  
  evalin('base', 'delete(''tmpeegrej.fdt'')');  
else
  elecIndices = find(reject == 0);
  indata = indata(:,elecIndices);
end;
times = times * length(find(reject ==0)) / length(reject);

return;
