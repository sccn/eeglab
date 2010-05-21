% rmart() - Remove eye artifacts from EEG data using regression with 
%           multiple time lags. Each channel is first made mean-zero. 
%           After JL Kenemans et al., Psychophysiology 28:114-21, 1991.
%
%   Usage: >> rmart('datafile','outfile',nchans,chanlist,eogchan,[threshold])
% Example: >> rmart('noisy.floats','clean.floats',31,[2:31],7)
%
%   Input:   datafile - input float data file, multiplexed by channel
%            outfile  - name of output float data file
%            nchans   - number of channels in datafile
%            chanlist - indices of EEG channel(s) to process (1,...,nchans)
%            eogchan  - regressing channel indices(s) (1,...,nchans)
%            threshold- abs threshold value to trigger regression {def|0 -> 80}
%
% Output: Writes [length(chanlist),size(data,2)] floats to 'outfile'
%
% Note: Regression epoch length and number of lags are set in the script.  
%       Some machines may require a new byte_order value in the script. 
%       note that runica() -> icaproj() should give better results! See
%       Jung et al., Psychophysiology 111:1745-58, 2000.
%
% Author: Tzyy-Ping Jung, SCCN/INC/UCSD, La Jolla, 1997 

% Copyright (C) 1997 Tzyy-Ping Jung, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 2-22-97  Tzyy-Ping Jung  CNL/Salk Institute, La Jolla, CA
% 2-24-97  Formatted for ICA package release -Scott Makeig
% 12-10-97 Changed name from rmartifact to rmart for toolbox inclusion -sm
% 12-11-97 Adapted to read/write a float matrix -sm & sw
% 09-14-00 Added comments and help -sm
% 01-25-02 reformated help & license -ad 
  
function rmart(datafile,outfile,nchans,chanlist,eogchan,threshold)

if nargin < 5
   help rmart
   return
end

%
% The following parameters may be fine-tuned for a data set
%
DEF_THRESHOLD = 80; % trigger regression on blocks exceeding this (default)
epoch         = 80; % remove artifacts in successive blocks of this length
nlags         = 40; % perform multiple regression filtering of this length

byte_order    = 'b';% (machine-dependent) byte order code for fopen();
MAKE_MEAN_ZERO = 1  % 1/0 flag removing mean offset from each channel

fprintf('Performing artifact regression on data in %s.\n',datafile);

if nargin<6 
   threshold = 0;
end
if threshold == 0,
   threshold = DEF_THRESHOLD;
end
fprintf('Regression threshold %g.\n',threshold);

%
% Read the input data
%
[fid,msg]=fopen(datafile,'r',byte_order); % open datafile
if fid < 3, 
  fprintf('rmart() - could not open data file: %s\n',msg);
  exit 1
end
data=(fread(fid,'float'))';
status=fclose('all');
if rem(length(data),nchans) == 0 % check length
   fprintf('rmart() - data length not divisible by %d chans.\n',nchans);
    return
end

data = reshape(data,nchans,length(data)/nchans);
[chans,frames] = size(data);
fprintf('Data of size [%d,%d] read.\n',chans,frames);
eog = data(eogchan,:);
data = data(chanlist,:);
procchans = length(chanlist);

fprintf('Regression epoch length %d frames.\n',epoch);
fprintf('Using %d regression lags.\n',nlags);
if length(eogchan)> 1
  fprintf('Processing %d of %d channels using %d EOG channels.\n',...
                        procchans,chans,length(eogchan));
else
  fprintf('Processing %d of %d channels using EOG channel %d.\n',...
                        procchans,chans,eogchan);
end

%
% Process the data
%
for i=1:procchans
  chan = chanlist(i);
  idx=[];
  frame=1+epoch/2+nlags/2;
  if MAKE_MEAN_ZERO
    data(chan,:) = data(chan,:) - mean(data(chan,:)); % make mean-zero
  end

  % Search the EOG & EEG records for values above threshold, 
  % Selected frame numbers are registered in the variable "idx". 
  % The entries in "idx" are at least epoch apart to avoid double 
  % compensation (regression) on the same portion of the EEG data.

  while frame <= length(eog)-epoch/2-nlags/2,  % foreach epoch in channel
    stop = min(frame+epoch-1,eogframes);
    tmp= ...
        find( abs(eog(frame:stop)) >= threshold ...
             | abs(data(chan,frame:stop)) >= threshold);
                                           % find beyond-threshold values
    if length(tmp) ~= 0
      mark = tmp(1)+frame-1;
      if  length(idx) ~= 0
        if mark-idx(length(idx))  < epoch,
           idx=[idx idx(length(idx))+epoch]; % To guarantee idx(i) & idx(i-1)
                                             % are at least EPOCH points apart
           frame = idx(length(idx))+epoch/2;
        else
           idx=[idx mark];
           frame = mark + epoch/2;
        end
      else
        idx=[idx mark];
        frame = mark + epoch/2;
      end 
    else
      frame=frame+epoch;
    end
  end % while

  % For each registered frame, take "epoch" points
  % surrounding it from the EEG, and "epoch + lag" points
  % from the EOG channel. Then perform multivariate 
  % linear regression on EEG channel.

  for j=1:length(idx);
     art=ones(1,epoch);
     eogtmp=eog(idx(j)-epoch/2-nlags/2:idx(j)+epoch/2-1+nlags/2);

     % Collect EOG data from lag/2 points before to lag/2 points 
     % after the regression window.
     for J=nlags:-1:1,
       art=[art ; eogtmp(J:J+epoch-1)];
     end
     eegtmp=data(chan,idx(j)-epoch/2:idx(j)+epoch/2-1);

     eegeog=eegtmp*art';        % perform the regression here
     eogeog=art*art';
     b=eegeog/eogeog;
     eegtmp=eegtmp-b*art;
     data(chan,idx(j)-epoch/2:idx(j)+epoch/2-1)=eegtmp;
   end % j
end % i

%
% Write output file
%
[fid,msg]=fopen(outfile,'w',byte_order);
if fid < 3
  fprintf('rmart() - could not open output file: %s\n',msg);
  return
end
count = fwrite(fid,data,'float');
if count == procchans*frames,
  fprintf('Output file "%s" written, size = [%d,%d] \n\n',...
             outfile,procchans,frames);
else
  fprintf('rmart(): Output file "%s" written, SIZE ONLY [%d,%g]\n',...
             outfile,procchans,count/procchans);
end
fclose('all');














