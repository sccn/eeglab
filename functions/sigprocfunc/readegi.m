% readegi() - read EGI Simple Binary datafile (versions 2,3,4,5,6,7).
%	      Return header info, EEG data, and any event data.
%
% Usage:
%   >> [head, TrialData, EventData] = readegi(filename)
%
% Inputs:
%   filename = EGI data filename
%
% Outputs:
%   head = struct containing header info (see readegihdr() )
%   TrialData = EEG channel data
%   EventData = event codes
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

function  [head, TrialData, EventData] = readegi(filename)

if nargin < 1,
    help readegi;
    return;
end

if nargout < 2 | nargout > 3,
	error('2 or 3 output args required');
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
    TrialData = zeros(FrameVals,head.segsamps*head.segments);
    readexpected = FrameVals*head.segsamps*head.segments; 
else
    TrialData = zeros(FrameVals,head.samples);
    readexpected = FrameVals*head.samples;
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
if (segmented),

    for i=1:head.segments,
	SegmentCatIndex(i) = fread(fid,1,'integer*2');
	SegmentStartTime(i) = fread(fid,1,'integer*4');

	[TrialData(:,[1+(i-1)*head.segsamps:i*head.segsamps]), count] = ...
	   fread(fid,[FrameVals,head.segsamps],datatype);

        readtotal = readtotal + count;
    end;

else
    % read unsegmented data
    [TrialData, readtotal] = fread(fid, [FrameVals,head.samples],datatype); 
end

fclose(fid);

if ~isequal(readtotal, readexpected)
     error('Number of values read not equal to number given in header.');
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