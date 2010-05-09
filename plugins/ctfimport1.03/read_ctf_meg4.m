function [meg] = read_ctf_meg4(fname, hdr, begsample, endsample, chanindx)

% READ_CTF_MEG4 reads specified samples from a CTF continous datafile
% It neglects all trial boundaries as if the data was acquired in
% non-continous mode.
%
% [meg] = read_ctf_meg4(filename, hdr, begsample, endsample)
%
% Where
%   filename	name of the datafile, including the .meg4 extension
%   header      with all data information (from read_ctf_meg4)
%   begsample   index of the first sample to read
%   endsample   index of the last sample to read
%   chanindx	index of channels to read (optional, default is all)
%
% Due to a non-dislosure agreement between CTF and the F.C. Donders Centre, 
% the allowed use of this function is limited. Do not distribute this function.

% This program is provided to users of CTF MEG systems as a courtesy only. Please
% do not redistribute it without permission from CTF Systems Inc.
% This program has no warranty whatsoever.

% Copyright (C) 2003, Robert Oostenveld
%
% use global flag for feedback
global fb
if isempty(fb)
  fb = 0;
end

if begsample<1
  error('cannot read before the start of the data');
end

if endsample>hdr.nSamples*hdr.nChans*hdr.nTrials
  error('cannot read beyond the end of the data');
end

if begsample>endsample
  error('cannot read a negative number of samples');
end

if nargin<5
  % select all channels
  chanindx = 1:hdr.nChans;
end

if isempty(chanindx)
  error('no channels were specified for reading CTF data')
end

fid = fopen(fname,'r','ieee-be');

if fid == -1
  error('could not open datafile');
end

CTFformat = setstr(fread(fid,8,'char'))';
if (strcmp(CTFformat(1,1:7),'MEG41CP')==0),
    error('datafile is not in CTF MEG4 format')
end 

% the data is not channel multiplexed, but stored in trials
% in each trial first all samples for channel 1 are given, then all samples of channel 2 ...
% this means that we have to read each channel for each trial

begtrial = ceil(begsample/hdr.nSamples);
endtrial = ceil(endsample/hdr.nSamples);

% this counts the files and offset if the 2GB file size boundary is encountered
multifilecount  = 0;
multifileoffset = 0;
fseek(fid, 0, 'eof');
multifilelength = ftell(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read only the selected data, channel-wise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if endtrial==begtrial
  rawbegsample = begsample - (begtrial-1)*hdr.nSamples;
  rawendsample = endsample - (begtrial-1)*hdr.nSamples;
  for chan=1:length(chanindx)
    % jump to the begin of this channel in this trial
    channeloffset = 8*(multifilecount+1)+(begtrial-1)*hdr.nSamples*4*hdr.nChans+(chanindx(chan)-1)*hdr.nSamples*4-multifileoffset;
    if channeloffset>=multifilelength
      % data goes beyond 2GB file boundary, jump to the next file
      channeloffset   = channeloffset - multifilelength+8;	% change the current offset, keep the header 
      multifileoffset = multifileoffset + multifilelength;	% remember for the next trial
      multifilecount  = multifilecount + 1;			% increase the file counter
      nextname = sprintf('%s.%d_meg4', fname(1:(end-5)), multifilecount);
      fprintf('data goes beyond 2GB file boundary, continuing with %s\n', nextname);
      fclose(fid);
      fid = fopen(nextname,'r','ieee-be');
      fseek(fid, 0, 'eof');
      multifilelength = ftell(fid);				% determine the length of the current file
    end
    fseek(fid, channeloffset, 'bof');
    % jump to the first sample of interest
    fseek(fid, 4*(rawbegsample-1), 'cof');
    % read the data from this channel
    [tmp, count] = fread(fid,[(endsample-begsample+1),1],'int32');
    raw(:,chan) = tmp;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read and concatenate the raw data of all trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  raw = zeros((endtrial-begtrial+1)*hdr.nSamples, length(chanindx));
  for trial=begtrial:endtrial

    if length(chanindx)==hdr.nChans & all(chanindx(:)'==1:hdr.nChans)
      % read the data from all channels
      rawbegsample = (trial-begtrial)*hdr.nSamples + 1;
      rawendsample = (trial-begtrial)*hdr.nSamples + hdr.nSamples;
      % jump to the begin of this trial
      trialoffset = 8*(multifilecount+1)+(trial-1)*hdr.nSamples*4*hdr.nChans-multifileoffset;
      if trialoffset>=multifilelength
        % data goes beyond 2GB file boundary, jump to the next file
        trialoffset   = trialoffset - multifilelength+8;	% change the current offset, keep the header 
        multifileoffset = multifileoffset + multifilelength;	% remember for the next trial
        multifilecount  = multifilecount + 1;			% increase the file counter
        nextname = sprintf('%s.%d_meg4', fname(1:(end-5)), multifilecount);
        fprintf('data goes beyond 2GB file boundary, continuing with %s\n', nextname);
        fclose(fid);
        fid = fopen(nextname,'r','ieee-be');
        fseek(fid, 0, 'eof');
        multifilelength = ftell(fid);				% determine the length of the current file
      end
      fseek(fid, trialoffset, 'bof');
      [tmp, count] = fread(fid,[hdr.nSamples,hdr.nChans],'int32');
      raw(rawbegsample:rawendsample,:) = tmp;

    else
      % read the data from the selected channels
      rawbegsample = (trial-begtrial)*hdr.nSamples + 1;
      rawendsample = (trial-begtrial)*hdr.nSamples + hdr.nSamples;
      for chan=1:length(chanindx)
        % jump to the begin of this channel in this trial
        channeloffset = 8*(multifilecount+1)+(trial-1)*hdr.nSamples*4*hdr.nChans+(chanindx(chan)-1)*hdr.nSamples*4-multifileoffset;
        if channeloffset>=multifilelength
          % data goes beyond 2GB file boundary, jump to the next file
          channeloffset   = channeloffset - multifilelength+8;	% change the current offset, keep the header 
          multifileoffset = multifileoffset + multifilelength;	% remember for the next trial
          multifilecount  = multifilecount + 1;			% increase the file counter
          nextname = sprintf('%s.%d_meg4', fname(1:(end-5)), multifilecount);
          fprintf('data goes beyond 2GB file boundary, continuing with %s\n', nextname);
          fclose(fid);
          fid = fopen(nextname,'r','ieee-be');
          fseek(fid, 0, 'eof');
          multifilelength = ftell(fid);				% determine the length of the current file
        end
        fseek(fid, channeloffset, 'bof');
        [tmp, count] = fread(fid,[hdr.nSamples,1],'int32');
        raw(rawbegsample:rawendsample,chan) = tmp;
      end %chan

    end %read all channels
  end %trial

  % select the raw data corresponding to the samples of interest
  rawoffset    = (begtrial-1)*hdr.nSamples;
  rawbegsample = begsample - rawoffset;
  rawendsample = endsample - rawoffset;
  raw = raw(rawbegsample:rawendsample, :);
end
fclose(fid);

% multiply the dimensionless values with the calibration value
gain = hdr.gainV(chanindx);	% only for selected channels
meg = raw';			% transpose the raw data
for i=1:size(meg,1)
  meg(i,:) = gain(i)*meg(i,:);
end

