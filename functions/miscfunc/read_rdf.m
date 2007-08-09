% read_rdf() - read RDF-formatted EEG files.
%
% Usage: 
%   >> [eeg,ev,header] = read_rdf(filename);
%
% Inputs:
%   filename - EEG data file in RDF format
%
% Outputs:
%   eeg - eeg data (array in size of [chan_no timepoint];
%   ev  - event structure
%      ev.sample_offset[] - event offsets in samples 
%                           from the first sample (0)
%      ev.event_code[]    - event codes (integers)
%   header - data structure for header information
%      header.ch_no     - number of channels
%      header.sample_no - number of samples

% Notes:
%
% Authors: Jeng-Ren Duann, CNL/Salk & INC/UCSD, 2002-12-12
%          with help from Andrey Vankov, creator of the RDF file format.

function [eeg,ev,header] = read_rdf(filename)
   
    if nargin < 1
        help read_rdf;
        return;
    end;
      
  eeg = [];
  ev = [];
  header = [];
  
  fp = fopen(filename,'rb','ieee-le');
  if fp == -1,
    disp('read_RDF(): Cannot open data file...!');
    return;
  end
  
  fseek(fp,6,-1);
  header.ch_no = fread(fp,1,'uint16');
  
  cnt = 0;
  ev_cnt = 0;
  while(~feof(fp)),
    tag = fread(fp,1,'uint32');
    if length(tag) == 0,
      break;
    end
    if tag == hex2dec('f0aa55'),
      cnt = cnt + 1;
      disp(['block ' num2str(cnt) ' found']);
      % read ch_no and block length
      fseek(fp,2,0);
      ch_no = fread(fp,1,'uint16');
      block_size = power(2,fread(fp,1,'uint16'));
      % read events
      fseek(fp,62,0);
      for i=1:110,
	samp_off = fread(fp,1,'char');
	cond_code = fread(fp,1,'char');
	ev_code = fread(fp,1,'uint16');
	if samp_off ~= 0,
	  ev_cnt = ev_cnt + 1;
	  ev(ev_cnt).sample_offset = samp_off + (cnt-1)*128;
	  ev(ev_cnt).event_code = ev_code;
	end
      end
      data = fread(fp,ch_no*block_size,'int16');
      data = reshape(data,ch_no,block_size);
      eeg = [eeg data];
    end
  end
  
  fclose(fp);
  header.sample_no = size(eeg,2);
