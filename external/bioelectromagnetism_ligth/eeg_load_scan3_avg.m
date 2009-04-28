function [avg] = eeg_load_scan3_avg(FILENAME)

% eeg_load_scan3_avg - Read a NeuroScan 3.x AVG File
%
% USEAGE:  [avg] = eeg_load_scan3_avg(FILENAME)
%
%   FILENAME     input Neuroscan .avg file (version 3.x)
%   avg          output data structure, with fields:
%   
%   avg.signal      - ERP signal (uV, Npnts x Mchan)
%   avg.variance    - variance of the signal (Npnts x Mchan)
%   avg.chan_names  - electrode labels
%   avg.pnts        - number of points in ERP waveform
%   avg.rate        - sample rate (Hz)
%   avg.xmin        - prestimulus epoch start (e.g., -100 msec)
%   avg.xmax        - poststimulus epoch end (e.g., 900 msec)
%   avg.nsweeps     - number of accepted trials/sweeps in avg
%   
%   e.g.
%   avg = eeg_load_scan3_avg( 'test.avg' );
%   plot( avg.signal );
%


% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $


% This program is distributed under the GNU GPL; you can redistribute 
% it and/or modify it. This program is distributed in the hope that it 
% will be useful, but WITHOUT ANY WARRANTY; without even the implied 
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% Version 1.1, arno_delorme@salk.edu
% Version 1.2, Darren.Weber_at_radiology.ucsf.edu
%
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

% UPDATE    version remark
% 1062001   0.1     primitive version based on a C program 
% 1092001   1.0     fully working version based on loadegg.m that I programmed 
% 1102001   1.1     adding channel names loading
% 03/2002   1.2     Darren.Weber_at_radiology.ucsf.edu
%                   - S_nsweeps_offset from 362 to 364;
%                     so that it finds ACCEPTED not TOTAL SWEEPS, which
%                     has a great impact on conversion to uV.
%                   - modified xmin/xmax to msec from sec.
%                   - modified output to structure with fields
%                   - Changed the name of the function (from loadavg).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1,
  help eeg_load_scan3_avg; return;
end;

fid = fopen(FILENAME,'r','ieee-le');

if fid<0,
  msg = sprintf('EEG_LOAD_SCAN3_AVG: Cannot find:\n... %s\n', FILENAME);
  error(msg);
end;

BOOL='int16';
ULONG='int32'; 
FLOAT='float32';

% read # of channels, # of samples, variance flag, and time bounds
% ----------------------------------------------------------------
S_nsweeps_offset    = 364; % Darren Weber modified this from 362
S_pnts_offset       = 368;
S_nchans_offset     = 370;
S_variance_offset   = 375;
S_rate_offset       = 376;
S_xmin_offset       = 505;
S_xmax_offset       = 509;
packed_sizeof_SETUP = 900;

fseek(fid, S_nsweeps_offset, 'bof');    avg.nsweeps = fread(fid, 1, 'ushort');
fseek(fid, S_pnts_offset, 'bof');       avg.pnts = fread(fid, 1, 'ushort');
fseek(fid, S_nchans_offset, 'bof');     chan = fread(fid, 1, 'ushort');
fseek(fid, S_variance_offset, 'bof');   variance_flag = fread(fid, 1, 'uchar');
fseek(fid, S_rate_offset, 'bof');       avg.rate = fread(fid, 1, 'ushort');
fseek(fid, S_xmin_offset, 'bof');       avg.xmin = fread(fid, 1, 'float32') * 1000;
fseek(fid, S_xmax_offset, 'bof');       avg.xmax = fread(fid, 1, 'float32') * 1000;
fseek(fid, packed_sizeof_SETUP, 'bof');

fprintf('number of channels : %d\n', chan);
fprintf('number of points   : %d\n', avg.pnts);
fprintf('sampling rate (Hz) : %f\n', avg.rate);
fprintf('xmin (msec)        : %f\n', avg.xmin);
fprintf('xmax (msec)        : %f\n', avg.xmax);
fprintf('Accepted sweeps    : %d\n', avg.nsweeps);

% read electrode configuration
% ----------------------------
fprintf('Electrode configuration\n');
for elec = 1:chan,
  channel_label_tmp = fread(fid, 10, 'uchar');
  avg.chan_names(elec,:) = channel_label_tmp';
  for index = 2:9,
    if avg.chan_names(elec,index) == 0,
      avg.chan_names(elec,index) = ' ';
    end;
  end;
  erp = fread(fid, 47-10, 'uchar');
  baseline(elec) = fread(fid, 1, 'ushort');
  erp = fread(fid, 10, 'uchar');
  sensitivity(elec) = fread(fid, 1, 'float32');
  erp = fread(fid, 8, 'uchar');
  calib(elec) = fread(fid, 1, 'float32');
  fprintf('%s: baseline: %d\tsensitivity: %f\tcalibration: %f\n', avg.chan_names(elec,1:4), baseline(elec), sensitivity(elec), calib(elec));
  factor(elec) = calib(elec) * sensitivity(elec) / 204.8;
end;

% Read signal data (amplifier units)
signal = zeros(avg.pnts, chan);
for elec = 1:chan,
  fseek(fid, 5, 'cof'); % skip sweeps header
  signal(:, elec) = fread(fid, avg.pnts, 'float32');
end;

if variance_flag,
  variance = zeros(avg.pnts, chan);
  for elec = 1:chan,
    variance(:, elec) = fread(fid, avg.pnts, 'float32');
  end;
  avg.variance = variance';
else
  avg.variance = [];
end;

% Convert signal to microvolts
baseline = repmat(baseline,avg.pnts,1);
calib    = repmat(calib,   avg.pnts,1);

if avg.nsweeps,
  signal = (signal - baseline) .* calib ./ avg.nsweeps;
else
  signal = (signal - baseline) .* calib;
end

avg.signal = signal';

fclose(fid);
return;
