% readegi() - read EGI Simple Binary datafile (versions 2,3,4,5,6,7).
%	      Return header info, EEG data, and any event data.
%
% Usage:
%   >> [head, TrialData, EventData, CatIndex] = readegi(filename, dataChunks)
%
% Required Input:
%   filename = EGI data filename
%
% Optional Input:
%   dataChunks = vector containing desired frame numbers(for unsegmented
%                datafiles) or segments (for segmented files). If this
%                input is empty or is not provided then all data will be
%                returned.
% 
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

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.15  2005/02/02 20:06:37  arno
% error msg
%
% Revision 1.14  2005/02/02 20:05:15  arno
% returning segment category indices
%
% Revision 1.13  2003/06/11 07:39:35  cooper
% Added 'dataChunks' input argument.
%
% Revision 1.12  2002/12/07 00:32:13  arno
% converting event codes to char
%
% Revision 1.11  2002/12/04 00:33:18  cooper
% Added conversion of A/D units to microvolts.
%
% Revision 1.10  2002/12/03 21:50:16  cooper
% readegi now compatible with EGI Simple Binary
% datafile versions 2,3,4,5,6, & 7.
%
% Revision 1.9  2002/11/26 20:16:00  cooper
% Added support for version 5 RAW datafiles.
% Added a check that the total # of data points
% read is equal to the # specified in the header.
%
% Revision 1.8  2002/11/14 19:00:13  scott
% help msg
%
% Revision 1.7  2002/11/14 18:01:20  arno
% adding messave
%
% Revision 1.6  2002/11/14 18:00:23  arno
% readEGIhr -> readegihdr
%
% Revision 1.5  2002/11/14 17:54:09  arno
% header typo
%
% Revision 1.4  2002/11/14 17:49:51  arno
% new cooper function
%
% Revision 1.3  2002/11/13 02:33:56  arno
% help
%
% Revision 1.2  2002/11/13 02:23:09  arno
% header ...
%

function  [head, TrialData, EventData, SegmentCatIndex] = readegi(filename, dataChunks)

if nargin <1 | nargin >2,
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
if (fid == 0),
   error(message);
end

% get our header structure
fprintf('Importing binary EGI data file ...\n');
head = readegihdr(fid);

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
	  if ismember(i,desiredSegments),
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
