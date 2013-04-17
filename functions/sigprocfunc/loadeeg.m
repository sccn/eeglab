% loadeeg() - load a binary data file in Neuroscan .eeg file format.
%
% Usage:
% >> signal = loadeeg(filename);
% >> [signal, accept, typeeeg, rt, response, chan_names, pnts, ...
%      ntrials, srate, xmin, xmax] = loadeeg( filename, chanlist, ...
%        triallist, typerange, accepttype, rtrange, responsetype);
%
% Inputs:
%   filename     - [string] Input Neuroscan .eeg file      
%   chanlist     - [integer array] Only import selected channels
%                  Ex: [3,4:10] {default: 'all'} 
%   triallist    - [integer array] Only import selected trials {default: import all}
%   typerange    - [integer array] Only import trials of selected type
%                  {default: import all}
%   accepttype   - [integer array] Only import trials with the selected
%                  'accept' field values {default: import all}
%   rtrange      - [float array] [min max] (ms) Only import trials with subject
%                  reaction times in this range {default: import all}
%   responsetype - [integer array] Only import trials with selected 
%                  response type values {default: import all}
%   format       - ['short'|'int32'|'auto'] data format. Neuroscan v4.3+ assume 32-bit data
%                  while older versions assume 16-bit. {default: 'auto' = Determine}
%
% Outputs:
%   signal       - Output signal of size (trials, points)	
%   accept       - [1/0] vector of values for the accept field (one per trial)
%   typeeeg      - [Integer] values for the accept type (size 1,trials) 
%   rt           - [float] values for the accept rt (size trials) 
%   response     - [Integer] values for the accept response (size 1,trials) 
%   chan_names   - ['string' array] channel names 
%   pnts         - Number of time points per trial
%   ntrials      - Number of trials
%   srate        - Sampling rate (Hz)
%   xmin         - Trial start latency (ms)
%   xmax         - Trial end latency (ms)
%
% Example:
%   % Load .eeg data into an array named 'signal' 
%   >> [signal]=loadeeg( 'test.eeg' ); 
%   % Plot the signal in the first channel, first trial
%   >> plot( signal(1,:) );	  		
%
% Author: Arnaud Delorme, CNL, Salk Institute, 2001
%
% See also: pop_loadeeg(), eeglab()

% .eeg binary file format
%   data are organised into an array of Number_of_electrode x (Number_of_points_per_trial*Number_of_sweeps)
%   for a file with 32 electrodes, 700 points per trial and 300 sweeps, the resulting array is 
%   of 32 collumn and 700*300 rows (300 consecutive blocs of 700 points) 	 	

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

function [signal, accept, typeeeg, rt, response, chan_names, pnts, ...
	nsweeps, rate, xmin, xmax]=loadeeg( FILENAME, chanlist, ...
	TrialList, typerange, acceptype, rtrange, responsetype, format)
if nargin<1 
	fprintf('Not enought arguments\n'); 
	help loadeeg 
	return;
end;
if nargin<2 chanlist='all'; end;
if nargin<3 TrialList='all'; end;
if nargin<4 typerange='all'; end;
if nargin<5 acceptype='all'; end;
if nargin<6 rtrange  ='all'; end;
if nargin<7 responsetype='all'; end;
if nargin<8 format='auto'; end;
format = lower(format);

if ~strcmpi(format, 'short') && ~strcmpi(format, 'int32') && ~strcmpi(format, 'auto'), error('loadeeg: format error'); end;

% open file for reading
% ---------------------
fid=fopen(FILENAME,'r','ieee-le');
if fid<0
	fprintf('Error LOADEEG: File %s not found\n', FILENAME);  
	return;
end;

% determine the actual file format if auto is requested
%------------------------------------------------------
if strcmpi(format, 'auto')
    format = 'short';
    fseek(fid,12,'bof');
    nextfile = fread(fid,1,'long');
    if (nextfile > 0)
        fseek(fid,nextfile+52,'bof');
        is32bit = fread(fid,1,'char');       
        if (is32bit == 1)
            format = 'int32'
        end;
    end;
    frewind(fid);
end;

% read general part of the erp header and set variables
% -----------------------------------------------------
rev = char(fread(fid, 20, 'uchar'))'; % revision number
erp = fread(fid, 342, 'uchar');	% skip the firsts 362 from BOF (368 bytes in orginal comment?)
nsweeps = fread(fid, 1, 'ushort');	% number of sweeps
erp = fread(fid, 4, 'uchar'); 	% skip 4 bytes 
pnts= fread(fid, 1, 'ushort');	% number of point per waveform
chan= fread(fid, 1, 'ushort');  % number of channels
erp = fread(fid, 4, 'uchar'); 	% skip 4 bytes 
rate= fread(fid, 1, 'ushort');  % sample rate (Hz)
erp = fread(fid, 127, 'uchar');	% skip 127 bytes 
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
	fprintf('%s: baseline: %d\tsensitivity: %f\tcalibration: %f\n', ...
            char(chan_names(elec,1:4)), baseline(elec), sensitivity(elec), calib(elec));
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
if isstr(chanlist) && strcmpi(chanlist, 'all'), chanlist = [1:chan]; end;
if isstr(TrialList) && strcmpi(TrialList, 'all'),	trialtagI     = 1; else trialtagI     = 0; end;
if isstr(acceptype) && strcmpi(acceptype, 'all'),	acceptagI     = 1; else acceptagI     = 0; end;
if isstr(typerange) && strcmpi(typerange, 'all'), 	typetagI      = 1; else typetagI      = 0; end;
if isstr(responsetype) && strcmpi(responsetype, 'all'),	responsetagI  = 1; else responsetagI  = 0; end;
if isstr(rtrange) && strcmpi(rtrange, 'all'), 	rttagI        = 1; else rttagI        = 0; end;

count_selected = 1;
fprintf('Reserving array (can take some time)\n');
if isstr(TrialList)
    signal = zeros( chan, pnts*nsweeps);
else
    signal = zeros( chan, pnts*length(TrialList));
end;
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
	if trialtagI trialtag = 1;        else trialtag = ismember_bc(sweep, TrialList); end;
	if acceptagI acceptag = 1;        else acceptag =  ismember(s_accept, acceptype); end;
	if typetagI  typetag  = 1; 	  else typetag  =  ismember(s_type, typerange); end;
	if responsetagI responsetag  = 1; else responsetag  = ismember_bc(s_response, responsetype); end;
	if rttagI       rttag  = 1; 	  else rttag  =  ismember(s_rt, rtrange); end;

	if typetag
		if trialtag
			if acceptag
				if responsetag
					if rttag

						buf = fread(fid, [chan pnts], format);
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
                        try
                            signal(:,[((count_selected-1)*pnts+1):count_selected*pnts]) = buf;
                            count_selected = count_selected + 1;
                            if not(mod(count_selected,10)) fprintf('%d sweeps selected out of %d\n', count_selected-1, sweep); end;
                        catch,
                            disp('Warning: File truncated, aborting file read.');
                            count_selected = count_selected + 1;
                            break;
                        end;
					end;
				end;
			end;
		end;
	end;

	if unreaded_buf 
        if strcmpi(format, 'short'), fseek(fid, buf_size*2, 'cof'); 
        else                         fseek(fid, buf_size*4, 'cof'); 
        end;
    end;				
end;
nsweeps = count_selected-1;
fclose(fid);

% restrincting array
% ---------------------------------------
fprintf('rereservation of variables\n');
signal = signal(chanlist, 1:(count_selected-1)*pnts);
chan_names = chan_names(chanlist,:);

return;


















