function [signal, accept, typeeeg, rt, response, chan_names, pnts, nsweeps, rate, xmin, xmax]=loadeeg( FILENAME, chanlist, TrialList, typerange, acceptype, rtrange, responsetype)
% eeg_load_scan_eeg - Load Neuroscan .EEG format
% 
% Usage: [signal, accept, typeeeg, rt, response, chan_names, pnts, nsweeps, rate, xmin, xmax]=loadeeg( FILENAME, chanlist, TrialList, typerange, acceptype, rtrange, responsetype)
%
%      FILENAME     input Neuroscan .avg file      
%      signal	    output signal	
%      variance     variance of the signal 
%      chan_names   array that represent the name of the electrodes
%
%      i.e. 
%	  [signal] = loadeeg( 'test.eeg' );     % load data into the array named 'signal'
%         plot( signal(1,:) );	  		% plot the signal for the first electrode of the first sweep
%
%      data are organised into an array of Number_of_electrode x (Number_of_points_per_trial*Number_of_sweeps)
%      for a file with 32 electrodes, 700 points per trial and 300 sweeps, the resulting array is 
%      of 32 collumn and 700*300 rows (300 consecutive blocs of 700 points) 	 	
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $

% Licence:  GNU GPL, no implied or express warranty
% History:  01/2001, arno_delorme@salk.edu
%	1062001		0.0		primitive version
%	1102001		1.0		fully working version
%	1112001		1.1		more parameters
%	1152001		1.2		fix bugs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 
	fprintf('Not enought arguments\n'); 
	help loadeeg 
	return;
end;
if nargin<2 CHAN='all'; end;
if nargin<3 TrialList='all'; end;
if nargin<4 typerange='all'; end;
if nargin<5 acceptype='all'; end;
if nargin<6 rtrange  ='all'; end;
if nargin<7 responsetype='all'; end;

% open file for reading
% ---------------------
fid=fopen(FILENAME,'r','ieee-le');
if fid<0
	fprintf(2,['Error LOADEEG: File ' FILENAME ' not found\n']);  
	return;
end;

% read general part of the erp header and set variables
% -----------------------------------------------------
erp = fread(fid, 362, 'uchar');	% skip the firsts 368 bytes
nsweeps = fread(fid, 1, 'ushort');	% number of sweeps
erp = fread(fid, 4, 'uchar'); 	% skip 4 bytes 
pnts= fread(fid, 1, 'ushort');	% number of point per waveform
chan= fread(fid, 1, 'ushort');  % number of channels
erp = fread(fid, 4, 'uchar'); 	% skip 4 bytes 
rate= fread(fid, 1, 'ushort');  % sample rate (Hz)
erp = fread(fid, 127, 'uchar');	% skip 125 bytes 
xmin= fread(fid, 1, 'float32'); % in s
xmax= fread(fid, 1, 'float32'); % in s
erp = fread(fid, 387, 'uchar');	% skip 387 bytes 

fprintf('number of channels         : %d\n', chan);
fprintf('number of points per trial : %d\n', pnts);
fprintf('sampling rate (Hz)         : %f\n', rate);
fprintf('xmin (s)                   : %f\n', xmin);
fprintf('xmax (s)                   : %f\n', xmax);

% read electrode configuration
% ----------------------------
fprintf('Electrode configuration\n');
for elec = 1:chan
   	channel_label_tmp = fread(fid, 10, 'uchar');
	chan_names(elec,:) = channel_label_tmp';
	for index = 2:9 if chan_names(elec,index) == 0 chan_names(elec,index)=' '; end; end;
	erp = fread(fid, 47-10, 'uchar');
	baseline(elec) = fread(fid, 1, 'ushort');
	erp = fread(fid, 10, 'uchar');
	sensitivity(elec) = fread(fid, 1, 'float32');
	erp = fread(fid, 8, 'uchar');
	calib(elec) = fread(fid, 1, 'float32');
	fprintf('%s: baseline: %d\tsensitivity: %f\tcalibration: %f\n', chan_names(elec,1:4), baseline(elec), sensitivity(elec), calib(elec));
	factor(elec) = calib(elec) * sensitivity(elec) / 204.8;
end;
%fprintf('Electrode configuration\n');
%for elec = 1:chan
%	erp = fread(fid, 47, 'uchar');
%	baseline(elec) = fread(fid, 1, 'ushort');
%	erp = fread(fid, 10, 'uchar');
%	sensitivity(elec) = fread(fid, 1, 'float32');
%	erp = fread(fid, 8, 'uchar');
%	calib(elec) = fread(fid, 1, 'float32');
%	fprintf('baseline: %d\tsensitivity: %f\tcalibration: %f\n', baseline(elec), sensitivity(elec), calib(elec));
%	factor(elec) = calib(elec) * sensitivity(elec) / 204.8;
%end;

xsize    = chan * pnts;
buf_size = chan * pnts ;			% size in shorts

% set tags for conditions
% -----------------------
if size(chanlist)  == size('all')	chanlist = [1:chan]; end;
if size(TrialList) == size('all')	trialtagI     = 1; else trialtagI     = 0; end;
if size(acceptype) == size('all')	acceptagI     = 1; else acceptagI     = 0; end;
if size(typerange) == size('all')	typetagI      = 1; else typetagI      = 0; end;
if size(responsetype) == size('all')	responsetagI  = 1; else responsetagI  = 0; end;
if size(rtrange)      == size('all')	rttagI        = 1; else rttagI        = 0; end;

count_selected = 1;
fprintf('Reserving array (can take some time)\n');
signal = zeros( chan, pnts*nsweeps);
fprintf('Array reserved, scanning file\n');

for sweep = 1:nsweeps

	% read sweeps header	
	% ------------------
	s_accept   = fread(fid, 1, 'uchar');
	s_type     = fread(fid, 1, 'ushort');
	s_correct  = fread(fid, 1, 'ushort');
	s_rt       = fread(fid, 1, 'float32');
	s_response = fread(fid, 1, 'ushort');
	s_reserved = fread(fid, 1, 'ushort');

	unreaded_buf = 1;

	% store the sweep or reject the sweep
	% -----------------------------------
	if trialtagI trialtag = 1;        else trialtag = ismember(sweep, TrialList); end;
	if acceptagI acceptag = 1;        else acceptag =  ismember(s_accept, acceptype); end;
	if typetagI  typetag  = 1; 	  else typetag  =  ismember(s_type, typerange); end;
	if responsetagI responsetag  = 1; else responsetag  = ismember(s_response, responsetype); end;
	if rttagI       rttag  = 1; 	  else rttag  =  ismember(s_rt, rtrange); end;

	if typetag
		if trialtag
			if acceptag
				if responsetag
					if rttag

						buf = fread(fid, [chan pnts], 'short');
						unreaded_buf = 0;

						% copy information to array
						% -------------------------
						accept(count_selected)   = s_accept;
						typeeeg(count_selected)  = s_type;
						rt(count_selected)       = s_rt;
						response(count_selected) = s_response;
		
						% demultiplex the data buffer and convert to microvolts
						% -----------------------------------------------------
						for elec = 1:chan
							buf(elec, :) = (buf(elec, :)-baseline(elec)-0.0)*factor(elec);
						end;
						signal(:,[((count_selected-1)*pnts+1):count_selected*pnts]) = buf;
						count_selected = count_selected + 1;
						if not(mod(count_selected,10)) fprintf('%d sweeps selected out of %d\n', count_selected-1, sweep); end;
					end;
				end;
			end;
		end;
	end;

	if unreaded_buf fseek(fid, buf_size*2, 'cof'); end;				
end;
nsweeps = count_selected-1;
fclose(fid);

% restrincting array
% ---------------------------------------
fprintf('rereservation of variables\n');
signal = signal(chanlist, 1:(count_selected-1)*pnts);
chan_names = chan_names(chanlist,:);

return;




% Frequency domain EEG File format
% 
% The frequency domain epoched EEG format has
% the same header as that  described in the appendix 
% of the SCAN manual (see the file sethead.h on the 
% download page). For each sweep of data, there is a 
% sweep header that is identical to that described 
% on page Headers-5.
% 
% At this point there is a difference.  After the sweep 
% header, the frequency domain data for each sweep has 
% the following format:
% 
% for each channel (erp.nchannels):
% 
% for each frequency bin (erp.pnts):
% 
% real value of the FFT stored as a 4-byte float;
% 
% for each frequency bin (erp.pnts):
% 
% imaginary value of the FFT stored as a 4-byte float;
% 
% The data have been scaled to microvolts prior to the FFT, 
% and it is these results which are stored to the file. 
