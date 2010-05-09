function [meg, hdr, hc] = read_ctf_ds(fname, trial)

% READ_CTF_DS reads specified trials from a CTF dataset
%
% [meg, hdr, hc] = read_ctf_ds(filename, trial)
%
% Where
%   filename	name of the dataset directory
%   trial	number of trials to read (can be multiple)
% and
%   meg		raw channel data
%   hdr		header information
%   hc		head coordinate system information
%
% Due to a non-dislosure agreement between CTF and the F.C. Donders Centre, 
% the allowed use of this function is limited. Do not distribute this function.

% This program is provided to users of CTF MEG systems as a courtesy only. Please
% do not redistribute it without permission from CTF Systems Inc.
% This program has no warranty whatsoever.

% Original author: Jim McKay November 1999
% Copyright (C) 1999-2000 CTF Systems Inc. All Rights Reserved.
%
% modifications Copyright (C) 2002, Ole Jensen
% modifications Copyright (C) 2003, Robert Oostenveld
%
% construct the filenames for this dataset
[path, name, ext] = fileparts(fname);
fnameRes = fullfile(path, [name '.ds' filesep name '.res4']);
fnameMeg = fullfile(path, [name '.ds' filesep name '.meg4']);
fnameHc  = fullfile(path, [name '.ds' filesep name '.hc']);

% read header information
hdr = read_ctf_res4(fnameRes);

% read head location information
hc = read_ctf_hc(fnameHc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the actual data
% this is from Ole's getTrialCTF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(fnameMeg,'r','ieee-be');

if fid == -1
    errMsg = strcat('Could not open data file:',fnameMeg);
    error(errMsg); 
end

CTFformat =setstr(fread(fid,8,'char'))';
if (strcmp(CTFformat(1,1:7),'MEG41CP')==0),
    fclose(fid)
    error('datafile is not in CTF MEG4 format')
end 

% assign trial information to output
meg.Trial = trial;

for i=1:length(trial)
  Trial = trial(i);
  fseek(fid, hdr.nSamples*hdr.nChans*4*(Trial-1)+8, 'bof');
  B = fread(fid,[hdr.nSamples,hdr.nChans],'int32');
  for j=1:hdr.nChans    
      B(:,j) = hdr.gainV(j)*B(:,j);
  end
  % assign channel data to output
  meg.data(i,:,:) = B';
end  

fclose(fid);

