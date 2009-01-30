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
% $Log: not supported by cvs2svn $
% Revision 1.1  2005/12/06 06:24:23  psdlw
% Alternative functions from the FieldTrip package, which is now released under GPL (so I assume these functions can be committed to the sourceforge cvs)
%
% Revision 1.6  2003/04/22 08:50:10  roberto
% fixed bug with array indexing
%
% Revision 1.5  2003/04/22 08:43:00  roberto
% fixed bug for trials that did not start at 1
%
% Revision 1.4  2003/03/24 12:34:33  roberto
% minor changes
%
% Revision 1.3  2003/03/14 10:45:12  roberto
% changed copyright from GPL to "commercial confident" to adhere to NDA
%
% Revision 1.2  2003/03/13 14:27:40  roberto
% moved res4 header part into separate function
%
% Revision 1.1  2003/03/12 16:20:47  roberto
% new implementation based on Ole Jensens code
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

