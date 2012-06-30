% loadavg() - loading eeg average data file from Neuroscan into
%             matlab.
% Usage:
%  >> [signal, variance, chan_names, ...
%               pnts, rate, xmin, xmax] = loadavg( filename, version );
% Input:
%      filename   - input Neuroscan .avg file   
%      version    - [1 or 2] function version. Default is 2 which scales
%                   the data properly.
%
% Output:
%      signal	  - output signal	
%      variance   - variance of the signal 
%      chan_names - array that represent the name of the electrodes
%      pnts       - number of data points
%      srate      - sampling rate
%      xmin       - ERP onset time
%      xmax       - ERP final time
%
% Example: 
%  % load data into the array named 'signal'
%  [signal]=loadavg( 'test.avg' );     
%  % plot the signal for the first electrode
%  plot( signal(1,:) );	  		
%
% Author: Arnaud Delorme (v1), Yang Zhang (v2)

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% Average binary file format
% Average data is stored as 4-byte floats in vectored format for each
% channel. Each channel has a 5-byte header that is no longer used. Thus,
% after the main file header, there is an unused 5-byte header followed by
% erp.pnts of 4-byte floating point numbers for the first channel; then a
% 5-byte header for channel two followed by erp.pnts*sizeof(float) bytes,
% etc. Therefore, the total number of bytes after the main header is:
% erp.nchannels * (5 + erp.pnts*sizeof(float)). To scale a data point to
% microvolts, multiply by the channel-specific calibration factor (i.e., for
% electrode j: channel[j]->calib) and divide by the number of sweeps in the
% average (i.e., channel[j]->n).

function [signal, variance, chan_names, pnts, rate, xmin, xmax]=loadavg( FILENAME,version)

if nargin<1 
	help loadavg 
	return; 
end;

if nargin < 2 || version == 2
    disp('Reading using method 2');
    [signal, variance, chan_names, pnts, rate, xmin, xmax]=loadavg_bcl( FILENAME );
    signal = signal';
    return;
end;
disp('Reading using method 1');

%if isempty(find(FILENAME=='.')) FILENAME=[FILENAME '.eeg']; end;

BOOL='int16';
ULONG='int32'; 
FLOAT='float32';
fid=fopen(FILENAME,'r','ieee-le');
if fid<0
	fprintf(2,['Error LOADEEG: File ' FILENAME ' not found\n']);  
	return;
end;

S_nsweeps_offset_total    = 362;
S_nsweeps_offset_accepted = 364;
S_pnts_offset 			  = 368;
S_nchans_offset 		  = 370;
S_variance_offset 		  = 375;
S_rate_offset 		  	  = 376;
S_xmin_offset 			  = 505;
S_xmax_offset 			  = 509;
packed_sizeof_SETUP       = 900;

% read general part of the erp header and set variables
% -----------------------------------------------------
%erp = fread(fid, 362, 'uchar');	% skip the firsts 368 bytes
%nsweeps = fread(fid, 1, 'ushort');	% number of sweeps
%erp = fread(fid, 4, 'uchar'); 	% skip 4 bytes 
%pnts= fread(fid, 1, 'ushort');	% number of point per waveform
%chan= fread(fid, 1, 'ushort');  % number of channels
%erp = fread(fid, 4, 'uchar'); 	% skip 4 bytes 
%rate= fread(fid, 1, 'ushort');  % sample rate (Hz)
%erp = fread(fid, 127, 'uchar');	% skip 125 bytes 
%xmin= fread(fid, 1, 'float32'); % in s
%xmax= fread(fid, 1, 'float32'); % in s
%erp = fread(fid, 387, 'uchar');	% skip 387 bytes 

% read # of channels, # of samples, variance flag, and real time bounds
% ---------------------------------------------------------------------
fseek(fid, S_nsweeps_offset_accepted, 'bof');  	nsweeps = fread(fid, 1, 'ushort');
if nsweeps == 0
    fseek(fid, S_nsweeps_offset_total, 'bof');  	nsweeps = fread(fid, 1, 'ushort');
end;
fseek(fid, S_pnts_offset, 'bof');  		pnts = fread(fid, 1, 'ushort');
fseek(fid, S_nchans_offset, 'bof');	  	chan = fread(fid, 1, 'ushort');
fseek(fid, S_variance_offset, 'bof');  	variance_flag = fread(fid, 1, 'uchar');
fseek(fid, S_rate_offset, 'bof');  		rate = fread(fid, 1, 'ushort');
fseek(fid, S_xmin_offset, 'bof');  		xmin = fread(fid, 1, 'float32');
fseek(fid, S_xmax_offset, 'bof');  		xmax = fread(fid, 1, 'float32');
fseek(fid, packed_sizeof_SETUP, 'bof');

fprintf('number of channels         : %d\n', chan);
fprintf('number of points per trial : %d\n', pnts);
fprintf('sampling rate (Hz)         : %f\n', rate);
fprintf('xmin (s)                   : %f\n', xmin);
fprintf('xmax (s)                   : %f\n', xmax);
fprintf('number of trials (s)       : %d\n', nsweeps);
if nsweeps == 0, nsweeps = 1; end;

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
	fprintf('%s: baseline: %d\tsensitivity: %f\tcalibration: %f\n', ...
            char(chan_names(elec,1:4)), baseline(elec), sensitivity(elec), calib(elec));
	factor(elec) = calib(elec) * sensitivity(elec) / 204.8;
end;

xsize    = chan * pnts;
buf_size = chan * pnts ;			% size in shorts

count_selected = 1;
fprintf('Reserving array (can take some time)\n');
signal = zeros( chan, pnts*nsweeps);
fprintf('Array reserved, scanning file\n');

signal   = zeros(pnts, chan);
variance = zeros(pnts, chan);

for elec = 1:chan

	% skip sweeps header and read data	
	% --------------------------------
	fseek(fid, 5, 'cof');
	signal(:, elec) = fread(fid, pnts, 'float32') * factor(elec) / nsweeps;
end;

if variance_flag
	for elec = 1:chan
		variance(:, elec) = fread(fid, pnts, 'float32');
	end;
else
	variance = 'novariance';
end;	

signal = signal';
variance = variance';

fclose(fid);
return;

function [signal,chan_names,variance, pnts, rate, xmin, xmax]=loadavg_bcl(FILENAME,chanNameList)
	%writed by allen zhang based on EEGLAB's loadavg.m
	%2009-11-25
	% NENU,CHANGCHUN,CHINA 
	% update by Yang Zhang
	%2011-2-28
	% NENU,Changchun, CHINA

	% revised by Yang Zhang
	% 2011-4-6
	%using automatic routine to identifiy the type of the avg file to get the true nsweeps parameter 

	if ~exist('chanNameList','var')
		chanNameList={'all'};
	end
	try
		fid=fopen(FILENAME,'r','ieee-le');
		% read general part of the erp header and set variables
		% -----------------------------------------------------
		fseek(fid, 362, 'cof');% skip the firsts 362 from BOF (368 bytes in orginal comment?)
		% disp(ftell(fid));% 360 bytes
		% hdr.nsweeps = fread(fid, 1, 'ushort');
		% disp(ftell(fid));% 362
		% hdr.compsweeps = fread(fid, 1, 'ushort');% the exact sweep numbers for eeg and single subject avg file| in grand avg file it represents the number of subjects
		% hdr.acceptcnt = fread(fid, 1, 'ushort');% number of accepted sweeps also the exact sweep numbers for grand avg file
		% hdr.rejectcnt = fread(fid, 1, 'ushort');%number of rejected sweeps
		% disp(ftell(fid));% 368

		compsweeps = fread(fid, 1, 'ushort');% the exact sweeps numbers for eeg file| in Grand avg file it represented number of subjects

		acceptcnt = fread(fid, 1, 'ushort');% number of accepted sweeps 
		rejectcnt = fread(fid, 1, 'ushort');%number of rejected sweeps

		% determine the type of avg file and choose the right value for nsweeps
		if (rejectcnt+acceptcnt)~=compsweeps
			disp('It''s a grand average file!!!');
			disp(['Subject number = ',num2str(compsweeps),'; sweeps = ',num2str(acceptcnt)]);
			nsweeps=compsweeps;
		else
			disp('It''s a single subject average file!!!');
			disp(['nsweeps = ',num2str(compsweeps),'; Accepted sweeps = ',num2str(acceptcnt),'; Rejected sweeps = ',num2str(rejectcnt)]);
			nsweeps=acceptcnt;
		end


		pnts=fread(fid, 1, 'ushort');	% number of point per waveform
		chan=fread(fid, 1, 'ushort');  % number of channels
		fseek(fid, 3, 'cof');	% skip 3 bytes
		variance_flag=fread(fid, 1, 'uchar');
		rate=fread(fid, 1, 'ushort');  % sample rate (Hz)
		fseek(fid, 127, 'cof');	% skip 127 bytes
		xmin=fread(fid, 1, 'float32'); % in s
		xmax=fread(fid, 1, 'float32'); % in s
		fseek(fid, 387, 'cof');	% skip 387 bytes

		% read electrode configuration
		% ----------------------------
		for elec = 1:chan
			channel_label_tmp = fread(fid, 10, 'uchar');
			electrodes(elec).tmp=channel_label_tmp;
			chan_names(elec,:) = channel_label_tmp'; %#ok<*AGROW>
			for index = 2:9
				if chan_names(elec,index) == 0
					chan_names(elec,index)=' ';
				end;
			end;
			fseek(fid, 61, 'cof');%skip 61 bytes
			electrodes(elec).calib= fread(fid, 1, 'float32');

			%     	erp = fread(fid, 47-10, 'uchar');
			% 	baseline(elec) = fread(fid, 1, 'ushort');
			% 	erp = fread(fid, 10, 'uchar');
			% 	sensitivity(elec) = fread(fid, 1, 'float32');
			% 	erp = fread(fid, 8, 'uchar');
			% 	    electrodes(elec).calib= fread(fid, 1, 'float32');
			%     
			% 	fprintf('%s: baseline: %d\tsensitivity: %f\tcalibration: %f\n', ...
			% char(chan_names(elec,1:4)), baseline(elec), sensitivity(elec), electrodes(elec).calib);
		end;

		signal   = zeros(pnts, chan);
		variance = zeros(pnts, chan);

		for elec = 1:chan
			% To scale a data point to
			% microvolts, multiply by the channel-specific calibration factor (i.e., for electrode j:
			% channel[j]->calib) and divide by the number of sweeps in the average (i.e.,
			% channel[j]->n);
			% skip sweeps header and read data
			% --------------------------------
			fseek(fid, 5, 'cof');
			signal(:, elec) =fread(fid, pnts, 'float32')*electrodes(elec).calib/nsweeps;

		end;

		if variance_flag
			for elec = 1:chan
				variance(:, elec) = fread(fid, pnts, 'float32')*electrodes(elec).calib/nsweeps;% not sure
			end;
		else
			variance = 'novariance';
		end;
		%%
		if ~strcmpi(chanNameList{1},'all')
			chanIDX=ones(1,chan);
			for ichanList=1:numel(chanNameList)
				for elec=1:chan
					if strcmpi(chanNameList{ichanList},char(chan_names(elec,1:numel(chanNameList{ichanList}))))
						chanIDX(elec)=0;
						break;
					end
				end
			end

			signal=signal(:,~chanIDX);
			if variance_flag
				variance=variance(:,~chanIDX);
			end
			chan_names=chan_names(~chanIDX,:);
		end
		% signal = signal';
		% variance = variance';
		fclose(fid);

	catch errorLOAD
		disp(FILENAME);
		rethrow(errorLOAD);
	end
	return;




