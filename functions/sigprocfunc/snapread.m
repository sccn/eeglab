% snapread() - Read data in Snap-Master Standard Binary Data File Format
%              Reads Snap-Master header and data matrix (nchans,nframes).
%              Ref: Users Guide, Snap-Master for Windows (1997) p. 4-19
% Usage:
%           >> data = snapread(filename);  % read .SMA file data
%           >> [data,params,events,head] = snapread(filename,seekframes); 
%                                          % save parameters, event list
% Inputs:
%        filename   = string containing whole filename of .SMA file to read
%        seekframes = skip this many initial time points {default 0}
% Output:
%        data   = data matrix, size(nchans,nframes)
%        params = [nchannels, nframes, srate]
%        events = vector of event frames (lo->hi transitions on event channel);
%                 See source to set EVENT_CHANNEL and EVENT_THRESH. 
%        head   = complete char file header 
%
% Authors: Scott Makeig & Tzyy-Ping Jung, SCCN/INC/UCSD, La Jolla, January 31, 2000 

% Copyright (C) January 31, 2000 from plotdata() Scott Makeig & Tzyy-Ping Jung, SCCN/INC/UCSD,
% scott@sccn.ucsd.edu
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

% Scott Makeig & Tzyy-Ping Jung, CNL / Salk Institute / January 31, 2000
% 7-13-00 fixed fprintf count print bug -sm
% 7-13-00 added arg seekframes -sm & ss
% 7-13-00 added test for file length -sm & ss
% 2-15-02 change close('all') to fclose('all') -ad
% 2-15-02 reformated help & license -ad 

function [data,params,events,head] = snapread(file,seekframes) 

EVENT_CHANNEL = 1;   % This channel assumed to store event pulses only!
EVENT_THRESH  = 2.3; % Default event threshold (may need to adjust!!!!)

if nargin < 2
 seekframes = 0;
end 
data   = [];
events = [];
params = [];
head   = [];
%
%%%%%%%%%%%%%%%%%%%%%%%%%% Open file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
fid = fopen(file,'r', 'ieee-le');
if fid<1
  fprintf('\nsnapread(): Could not open file %s.\n\n',file)
  return
else
  fprintf('\n  Opened file %s for reading...\n',file);
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%% Read header %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read header and extract info

numbegin=0;
head = [];
while ~numbegin,
  line =fgets(fid);
  head = [head line];
  if (length(line)>=8 && all(line(1:8)=='"NCHAN%"'))
      nchans=str2num(line(findstr(line,'=')+1:end-1));
  end
  if (length(line)>= 12 && all(line(1:12)=='"NUM.POINTS"'))
      nframes=str2num(line(findstr(line,'=')+1:end-1));
  end
  if (length(line)>= 10 && all(line(1:10)=='"ACT.FREQ"'))
      srate=str2num(line(findstr(line,'=')+1:end-1));
  end
  if (length(line)>= 4 && all(line(1:4)=='"TR"'))
     head = head(1:length(head)-length(line));
     line =fgets(fid); % get the time and date stamp line
     numbegin=1;
  end
end
params = [nchans, nframes, srate];
fprintf('  Number of channels:    %d\n',nchans);
fprintf('  Number of data points: %d\n',nframes);
fprintf('  Sampling rate:         %3.1f Hz\n',srate);

fseek(fid,1,0); % skip final hex-AA char

%
%%%%%%%%%%%%%%%%%%% Test for correct file length %%%%%%%%%%%%%%%%%%%%
%
datstart = ftell(fid); % save current position in file
fseek(fid,0,'eof');    % go to file end
datend = ftell(fid);   % report position
discrep = (datend-datstart) - nchans*nframes*4;
if discrep ~= 0
   fprintf(' ***> Note: File length does not match header information!\n')
   fprintf('            Difference: %g frames\n',discrep/(4*nchans));
end
fseek(fid,datstart,'bof'); % seek back to data start
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read data %%%%%%%%%%%%%%%%%%%%%%%%%%
%

if seekframes > 0
  fprintf('  Omitting first %d frames of data.\n',seekframes)
  fprintf('moving %d bytes\n',seekframes*nchans*4);
  fsq = fseek(fid,seekframes*nchans*4,'cof')
  if fsq<0
    fsqerr = ferror(fid);
    fprintf('fseek() error: %s\n',fsqerr);
    return
  end
end

[data count] = fread(fid,[nchans,nframes],'float'); % read data frame
if count<nchans*(nframes-seekframes)
  fprintf('\nsnapread(): Could not read %d data frames.\n             Only %g read.\n\n',...
       nframes,count/nchans);
  return
else
  fprintf('\n  Read %d frames, each one sync channel(%d) and %d data channels.\n',...
       nframes,EVENT_CHANNEL,nchans-1);
end

if nargout>2
  fprintf('  Finding event transitions (>%g) on channel %d ...\n',...
                             EVENT_THRESH,EVENT_CHANNEL);
  events = zeros(EVENT_CHANNEL,nframes);
  for n=2:nframes
     if abs(data(EVENT_CHANNEL,n-1)) < EVENT_THRESH ...
        && abs(data(EVENT_CHANNEL,n)) > EVENT_THRESH
             events(n) = 1;
     end
  end
  fprintf('  Total of %d events found.\n',sum(events));
  data = data(2:nchans,:);  % eliminate channel 1
  fprintf('  Event channel %d removed from output.\n',EVENT_CHANNEL);
  params(1) = nchans-1;
end
fprintf('\n');
fclose('all');
