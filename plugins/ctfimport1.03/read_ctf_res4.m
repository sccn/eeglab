function hdr = read_ctf_res4(fname)

% READ_CTF_RES4 reads the header in RES4 format from a CTF dataset
%
% hdr = read_ctf_res4(filename)
%
% Due to a non-dislosure agreement between CTF and the F.C. Donders Centre, 
% the allwoed use of this function is limited. Do not distribute this function.

% This program is provided to users of CTF MEG systems as a courtesy only. Please
% do not redistribute it without permission from CTF Systems Inc.
% This program has no warranty whatsoever.

% Author(s): Jim McKay November 1999
% Last revision: Jim McKay
% Copyright (c) 1999-2000 CTF Systems Inc. All Rights Reserved.
%
% modifications Copyright (C) 2002, Ole Jensen
% modifications Copyright (C) 2003, Robert Oostenveld
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read header information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(fname,'r','b');

% Check if header file exist
if fid == -1
    errMsg = strcat('Could not open header file:',fname);
    error(errMsg); 
end

% first 8 bytes contain filetype

% Check is fileformat is correct
r_head=setstr(fread(fid,8,'char'))';
if (strcmp(r_head(1,1:7),'MEG41RS')==0),
  fclose(fid)
  errMsg = strcat('Resource file is not in CTF MEG4 format in file:',fname)
  error(errMsg); 
end %if

% Read the initial parameters 
appName       = setstr(fread(fid,256,'char'))' ;
dataOrigin    = setstr(fread(fid,256,'char'))' ;
dataDescrip   = setstr(fread(fid,256,'char'))' ;
no_trial_avgd = fread(fid,1,'int16')          ;
data_time     = setstr(fread(fid,255,'char'))';
data_date     = setstr(fread(fid,255,'char'))';

fseek(fid,1288,'bof');
% Read the general recording parameters
no_samples  = fread(fid,1,'int32');
no_channels = fread(fid,1,'int16');
fseek(fid,2,'cof');			% hole of 2 bytes due to improper alignment
sample_rate = fread(fid,1,'double');
epoch       = fread(fid,1,'double');
no_trials   = fread(fid,1,'int16');
fseek(fid,2,'cof');			% hole of 2 bytes due to improper alignment
preTrigpts=fread(fid,1,'int32');

fseek(fid,1360,'bof');
% read in the meg4Filesetup structure
run_name     = setstr(fread(fid,32,'char')');
run_title    = setstr(fread(fid,256,'char')');
instruments  = setstr(fread(fid,32,'char')');
coll_desc    = setstr(fread(fid,32,'char')');
subj_id      = setstr(fread(fid,32,'char')');
operator     = setstr(fread(fid,32,'char')') ;
sensFilename = setstr(fread(fid,60,'char')') ;

fseek(fid,1839,'bof');
% Read in the run description length 
rd_len=fread(fid,1,'uint8');
% Go to the run description and read it in
fseek(fid,1844,'bof');
run_desc=setstr(fread(fid,rd_len,'char')');

% read in the filter information
temp=fread(fid,2,'uint8');
num_filt=temp(2);
for fi=0:(num_filt-1),
  filt_info=fread(fid,18,'uint8');
  num_fparm=filt_info(18);
  if num_fparm ~= 0,
    filt_parm=fread(fid,8*num_fparm,'uint8');
  end % if
end % for fi

% Read in the channel names
for i=1:no_channels,
    temp=fread(fid,32,'char')';
    temp(find(temp<32 )) = ' ';		% remove non-printable characters
    temp(find(temp>126)) = ' ';		% remove non-printable characters
    endstr = findstr(temp, '-'); temp(endstr:end) = ' ';	% cut off at '-'
    endstr = findstr(temp, ' '); temp(endstr:end) = ' ';	% cut off at ' '
    chan_name(i,:) = char(temp);		% as char array
    chan_label{i}  = deblank(char(temp));	% as cell array
end %for

% pre-allocate some memory space
sensGain = zeros([no_channels,1]);
qGain    = zeros([no_channels,1]);
ioGain   = zeros([no_channels,1]);
sensType = zeros([no_channels,1]);

% Read in the sensor information
fp = ftell(fid);
for chan=1:no_channels,
    fread(fid,1,'uint8');			% Read and ignore 1 byte from enum
    sensType(chan)=fread(fid,1,'uint8');	% Read sensor type
    fread(fid,2,'uint8');			% Read and ignore originalRunNum
    fread(fid,4,'uint8');			% Read and ignore coilShape
    sensGain(chan)=fread(fid,1,'double');	% Read sensor gain in Phi0/Tesla
    qGain(chan)=fread(fid,1,'double');		% Read qxx gain (usually 2^20 for Q20)
    ioGain(chan)=fread(fid,1,'double');		% Read i/o gain of special sensors (usually 1.0)
    ioOffset(chan)=fread(fid,1,'double');
    numCoils(chan)=fread(fid,1,'int16');
    grad_order_no(chan)=fread(fid,1,'int16');
    fread(fid,4,'char');

    % read the coil positions and orientations
    for i=1:8
      Chan(chan).coil(i).pos = fread(fid,3,'double')';	
      fread(fid,1,'double');
      Chan(chan).coil(i).ori = fread(fid,3,'double')';
      fread(fid,3,'double');
    end

    % read the coil positions and orientations in head coordinates(?)
    for i=1:8
      Chan(chan).coilHC(i).pos = fread(fid,3,'double')';	
      fread(fid,1,'double');
      Chan(chan).coilHC(i).ori = fread(fid,3,'double')';
      fread(fid,3,'double');
    end

    % jump to the next sensor info record
    fseek(fid, fp+chan*1328, 'bof');
end % for chan

% close the header file
fclose(fid);

% determine the different channel types
rowMEG  = [];
rowEEG  = [];
rowTRIG = [];
rowREF  = [];
rowALL  = 1:no_channels;
for k=1:no_channels
    if strcmp(char(chan_name(k,1)),'M')
        rowMEG = [rowMEG k];
    end
    if strcmp(char(chan_name(k,1:3)),'EEG')
        rowEEG = [rowEEG k];
    end
    if strcmp(char(chan_name(k,1:4)),'STIM')
        rowTRIG = [rowTRIG k];
    end
    if chan_name(k,1)=='B' | chan_name(k,1)=='G' | chan_name(k,1)=='P' | chan_name(k,1)=='Q' | chan_name(k,1)=='R'
        rowREF = [rowREF k];
    end
end

% assign all the variables that should be outputted as header information
hdr.Fs          = sample_rate;
hdr.nChans      = no_channels;
hdr.nSamples    = no_samples;
hdr.nSamplesPre = preTrigpts;
hdr.timeVec     = (1:no_samples)/sample_rate - preTrigpts/sample_rate;
hdr.nTrials     = no_trials;
hdr.gainV       = ioGain./(qGain.*sensGain);
hdr.label       = chan_label(:);
hdr.nameALL     = chan_name; 
hdr.rowMEG      = rowMEG;
hdr.rowEEG      = rowEEG;
hdr.rowTRIG     = rowTRIG;
hdr.rowREF      = rowREF;
hdr.Chan        = Chan;
% hdr.nameEEG     = [];
% hdr.nameMEG     = [];
% hdr.nameEOG     = [];
% hdr.trigV       = [];
% hdr.SwapData    = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the gradiometer system in DEWAR coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% combine the bottom and top coil of the MEG channels into hardware gradiometers
numMEG = length(hdr.rowMEG);
for i=1:numMEG
  coil1_pos(i,:) = hdr.Chan(hdr.rowMEG(i)).coil(1).pos;
  coil1_ori(i,:) = hdr.Chan(hdr.rowMEG(i)).coil(1).ori;
  coil2_pos(i,:) = hdr.Chan(hdr.rowMEG(i)).coil(2).pos;
  coil2_ori(i,:) = hdr.Chan(hdr.rowMEG(i)).coil(2).ori;
end

% apparently, some coils are oriented to the wrong side
tmp = coil2_pos - coil1_pos;
sel = find(dot(coil1_ori, tmp, 2)<1);
coil1_ori(sel,:) = -coil1_ori(sel,:);
coil2_ori(sel,:) = -coil2_ori(sel,:);

gradDEWAR.pnt = [coil1_pos; coil2_pos];
gradDEWAR.ori = [coil1_ori; coil2_ori];
gradDEWAR.tra = [eye(numMEG) -eye(numMEG)];
gradDEWAR.label = hdr.label(rowMEG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the gradiometer system in HEAD coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% combine the bottom and top coil of the MEG channels into hardware gradiometers
numMEG = length(hdr.rowMEG);
for i=1:numMEG
  coil1_pos(i,:) = hdr.Chan(hdr.rowMEG(i)).coilHC(1).pos;
  coil1_ori(i,:) = hdr.Chan(hdr.rowMEG(i)).coilHC(1).ori;
  coil2_pos(i,:) = hdr.Chan(hdr.rowMEG(i)).coilHC(2).pos;
  coil2_ori(i,:) = hdr.Chan(hdr.rowMEG(i)).coilHC(2).ori;
end

% apparently, some coils are oriented to the wrong side
tmp = coil2_pos - coil1_pos;
sel = find(dot(coil1_ori, tmp, 2)<1);
coil1_ori(sel,:) = -coil1_ori(sel,:);
coil2_ori(sel,:) = -coil2_ori(sel,:);

gradHC.pnt = [coil1_pos; coil2_pos];
gradHC.ori = [coil1_ori; coil2_ori];
gradHC.tra = [eye(numMEG) eye(numMEG)];
gradHC.label = hdr.label(rowMEG);

hdr.grad      = gradHC;		% default is in head-coordinates
hdr.grad.unit = 'cm';		% this is default in all CTF software

