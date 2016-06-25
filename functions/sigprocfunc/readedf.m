% readedf() - read eeg data in EDF format.
%
% Usage: 
%    >> [data,header] = readedf(filename);
%
% Input:
%    filename - file name of the eeg data
% 
% Output:
%    data   - eeg data in (channel, timepoint)
%    header - structured information about the read eeg data
%      header.length - length of header to jump to the first entry of eeg data
%      header.records - how many frames in the eeg data file
%      header.duration - duration (measured in second) of one frame
%      header.channels - channel number in eeg data file
%      header.channelname - channel name
%      header.transducer - type of eeg electrods used to acquire
%      header.physdime - details
%      header.physmin - details
%      header.physmax - details
%      header.digimin - details
%      header.digimax - details
%      header.prefilt - pre-filterization spec
%      header.samplerate - sampling rate
%
% Author: Jeng-Ren Duann, CNL/Salk Inst., 2001-12-21

% Copyright (C) Jeng-Ren Duann, CNL/Salk Inst., 2001-12-21
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

% 03-21-02 editing header, add help -ad 

function [data,header] = readedf(filename);

if nargin < 1
    help readedf;
    return;
end;
    
fp = fopen(filename,'r','ieee-le');
if fp == -1,
  error('File not found ...!');
  return;
end

hdr = setstr(fread(fp,256,'uchar')');
header.length = str2num(hdr(185:192));
header.records = str2num(hdr(237:244));
header.duration = str2num(hdr(245:252));
header.channels = str2num(hdr(253:256));
header.channelname = setstr(fread(fp,[16,header.channels],'char')');
header.transducer = setstr(fread(fp,[80,header.channels],'char')');
header.physdime = setstr(fread(fp,[8,header.channels],'char')');
header.physmin = str2num(setstr(fread(fp,[8,header.channels],'char')'));
header.physmax = str2num(setstr(fread(fp,[8,header.channels],'char')'));
header.digimin = str2num(setstr(fread(fp,[8,header.channels],'char')'));
header.digimax = str2num(setstr(fread(fp,[8,header.channels],'char')'));
header.prefilt = setstr(fread(fp,[80,header.channels],'char')');
header.samplerate = str2num(setstr(fread(fp,[8,header.channels],'char')'))./header.duration;

fseek(fp,header.length,-1);
data = fread(fp,'int16');
fclose(fp);

data = reshape(data,header.duration*header.samplerate(1),header.channels,header.records);
temp = [];
for i=1:header.records,
  temp = [temp data(:,:,i)'];
end
data = temp;
