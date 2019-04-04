% readegihdr() - read header information from EGI (versions 2,3,4,5,6,7) data file.
%
% Usage:
%   >> [head] = readegihdr(fid,forceversion)
%
% Input:
%   fid - file identifier of EGI datafile
%   forceversion - optional integer input to override automatic reading of version number.
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

% Copyright (C) 2002 , Salk Institute, arno@salk.edu
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


function head = readegihdr(fid,forceversion)

if nargin < 1
    help readegihdr;
    return;
end
    
head.version = fread(fid,1,'integer*4');
if exist('forceversion')
    head.version = forceversion;
end

if ~( head.version >= 2 && head.version <= 7 ),
        error('EGI Simple Binary Versions 2-7 supported only.');
end

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
