% readegi() - read EGI Simple Binary datafile (versions 2,3,4,5,6,7).
%	      Return header info, EEG data, and any event data.
% Usage:
%   >> [head, TrialData, EventData, CatIndex] = readegi(filename, dataChunks, forceversion)
%
% Required Input:
%   filename = EGI data filename
%
% Optional Input:
%   dataChunks = vector containing desired frame numbers(for unsegmented
%                datafiles) or segments (for segmented files). If this
%                input is empty or is not provided then all data will be
%                returned.
%   forceversion = integer from 2 to 7 which overrides the EGI version read from the
%                file header.  This has been occasionally necessary in cases where
%                the file format was incorrectly indicated in the header.
% Outputs:
%   head = struct containing header info (see readegihdr() )
%   TrialData = EEG channel data
%   EventData = event codes
%   CatIndex  = segment category indices
%
% Author: Cooper Roddey, SCCN, 13 Nov 2002
%
% Note: code derived from C source code written by
%       Tom Renner at EGI, Inc.
%
% See also: readegihdr()

% Copyright (C) 2002 , Salk Institute, arno@salk.edu
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

function  [head, TrialData, EventData, SegmentCatIndex] = readegi(filename, dataChunks,forceversion)

if nargin <1 | nargin > 3,
    help readegi;
    return;
end

if nargout < 2 | nargout > 4,
	error('2 to 4 output args required');
end

if ~exist('dataChunks','var')
     dataChunks = [];
end

if ~isvector(dataChunks)
     error('dataChunks must either be empty or a vector');
end

[fid,message] = fopen(filename,'rb','b');
if (fid < 0),
   error(message);
end

% get our header structure
fprintf('Importing binary EGI data file ...\n');
if exist('forceversion')
   head = readegihdr(fid,forceversion);
else
   head = readegihdr(fid);
end

% do we have segmented data?
segmented = 0;
switch(head.version),
    case {2,4,6}
       segmented = 0;
    case {3,5,7}
       segmented = 1;
end

% each type of event has a dedicated "channel"
FrameVals = head.nchan + head.eventtypes;

if (segmented)

    desiredSegments = dataChunks;
    if isempty(desiredSegments),
       desiredSegments = [1:head.segments];
    end

    nsegs = length(desiredSegments);

    TrialData = zeros(FrameVals,head.segsamps*nsegs);
    readexpected = FrameVals*head.segsamps*nsegs; 

else

    desiredFrames = dataChunks;
    if isempty(desiredFrames),
       desiredFrames = [1:head.samples];
    end
    
    nframes = length(desiredFrames);

    TrialData = zeros(FrameVals,nframes);
    readexpected = FrameVals*nframes;

end

% get datatype from version number
switch(head.version)
   case {2,3}
      datatype = 'integer*2';
   case {4,5} 
      datatype = 'float32';
   case {6,7}
      datatype = 'float64';
   otherwise,
      error('Unknown data format');
end

% read in epoch data
readtotal = 0;
j=0;
SegmentCatIndex = [];
if (segmented),

    for i=1:head.segments,
	segcatind = fread(fid,1,'integer*2');
	segtime = fread(fid,1,'integer*4');

	[tdata, count] = ...
	   fread(fid,[FrameVals,head.segsamps],datatype);

       % check if this segment is one of our desired ones
        if ismember(i,desiredSegments) 

           j=j+1;
           SegmentCatIndex(j) = segcatind;
           SegmentStartTime(j) = segtime;
 
           TrialData(:,[1+(j-1)*head.segsamps:j*head.segsamps]) = tdata;
   
           readtotal = readtotal + count;
        end

	if (j >= nsegs), break; end;

    end;

else
    % read unsegmented data

    % if dataChunks is empty, read all frames
    if isempty(dataChunks)

        [TrialData, readtotal] = fread(fid, [FrameVals,head.samples],datatype);     
  
    else   % grab only the desiredFrames
           % This could take a while...

        for i=1:head.samples,

          [tdata, count] = fread(fid, [FrameVals,1],datatype);
  
          % check if this segment is a keeper
          if ismember(i,desiredFrames),
              j=j+1;
              TrialData(:,j) = tdata;
              readtotal = readtotal + count;
          end

          if (j >= nframes), break; end;

        end
    end
end

fclose(fid);

if ~isequal(readtotal, readexpected)
     error('Number of values read not equal to the number expected.');
end 

EventData = [];
if (head.eventtypes > 0),
	EventData = TrialData(head.nchan+1:end,:);
	TrialData = TrialData(1:head.nchan,:);
end 

% convert from A/D units to microvolts
if ( head.bits ~= 0 & head.range ~= 0 )
       TrialData = (head.range/(2^head.bits))*TrialData;
end

% convert event codes to char
% ---------------------------
head.eventcode = char(head.eventcode);


%--------------------------- isvector() ---------------
function retval = isvector(V)

   s = size(V);
   retval = ( length(s) < 3 ) & ( min(s) <= 1 );

%------------------------------------------------------
