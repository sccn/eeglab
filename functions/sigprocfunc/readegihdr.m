% readegihdr() - read header information from EGI (versions 2,3,4,5,6,7) data file.
%
% Usage:
%   >> [head] = readegihdr(fid)
%
% Input:
%   fid - file identifier of EGI datafile
%
% Output:
%   head - structure containing header information.
%          Structure fields are:
%           version     = 2,3,4,5,6,or 7
%           samp_rate   = sampling rate in samples/s
%           nchan       = number of EEG channels
%           gain        = gain of amplifier
%           bits        = number of bits/sample
%           range       = abs(max. value)
%           segments    = number of epochs
%           categories  = number of categories
%           catname     = cell array of category names
%           segsamps    = number of samples/segment
%           eventtypes  = number of event types
%           eventcode   = string array of event codes
%
% Author: Cooper Roddey, SCCN, 13 Nov 2002
%
% Note: this code derived from C source code written by 
%       Tom Renner at EGI, Inc. (www.egi.com)
%
% See also: readegi()

function head = readegihdr(fid)

if nargin < 1
    help readegihdr;
    return;
end;
    
head.version = fread(fid,1,'integer*4');

if ~( head.version >= 2 & head.version <= 7 ),
        error('EGI Simple Binary Versions 2-7 supported only.');
end;

year = fread(fid,1,'integer*2');
month = fread(fid,1,'integer*2');
day = fread(fid,1,'integer*2');
hour = fread(fid,1,'integer*2');
minute = fread(fid,1,'integer*2');
second = fread(fid,1,'integer*2');
millisecond = fread(fid,1,'integer*4');

head.samp_rate = fread(fid,1,'integer*2');
head.nchan = fread(fid,1,'integer*2');
head.gain = fread(fid,1,'integer*2');
head.bits = fread(fid,1,'integer*2');
head.range = fread(fid,1,'integer*2');

head.samples = 0;
head.segments = 0;
head.segsamps = 0;
head.eventtypes = 0;
head.categories = 0;
head.catname = {};
head.eventcode = '';

switch (head.version)
   case {2,4,6}
        head.samples = fread(fid, 1 ,'integer*4');
   case {3,5,7} 
        head.categories = fread(fid,1,'integer*2');
        if (head.categories),
                for i=1:head.categories,
                        catname_len(i) = fread(fid,1,'uchar');
                        head.catname{i} = char(fread(fid,catname_len(i),'uchar'))';
                end
        end
        head.segments = fread(fid,1,'integer*2');
        head.segsamps = fread(fid,1,'integer*4');
   otherwise
        error('Invalid EGI version');
end

head.eventtypes = fread(fid,1,'integer*2');

if isequal(head.eventtypes,0),
	head.eventcode(1,1:4) = 'none';
else
        for i = 1:head.eventtypes,
                head.eventcode(i,1:4) = fread(fid,[1,4],'uchar');
        end
end
