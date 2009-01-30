function [ctf] = ctf_read_meg4(folder,ctf,CHAN,TIME,TRIALS,COEFS);

% ctf_read_meg4 - read meg4 format data from a CTF .ds folder
%
% [ctf] = ctf_read_meg4([folder],[ctf],[CHAN],[TIME],[TRIALS]);
%
% This function reads all or select portions of the raw meg data matrix in
% the .meg4 file within any .ds folder.  It may call the ctf_read_res4
% function to identify the relevant parameters of the dataset.
%
% The .meg4 file contains the raw numbers sampled from the electronics. In
% this function, these raw analog2digial numbers are multiplied by the
% appropriate sensor gains, which are read from the .res4 file.  However,
% note that the data values returned can be very small (10^-12 Tesla),
% which may be a problem for some computations.
%
% INPUTS
%
% If you do not wish to specify an input option, use [], but keep the order
% of the input options as above.  Only specify as many input options as
% required.  With no input options, the function will prompt for a folder,
% call ctf_read_res4 and then read all of the data matrix.
%
% folder - the directory of the .ds data set to read.  By
% default, a gui prompts for the folder.
%
% ctf - a struct with setup, sensor and data fields.  If the setup field is
% missing or empty, this function calls ctf_read_res4.
%
% CHAN - a integer array of channel numbers to read.
%        eg, [30:35] reads channels 30 to 35.  Also
%        If CHAN = 'eeg', read only eeg channels/sensorIndices
%        If CHAN = 'meg', read only meg channels/sensorIndices
%        If CHAN = 'ref', read only reference channels/sensorIndices
%        If CHAN = 'other', read only the other channels/sensorIndices
%        If CHAN = 'megeeg', read meg and eeg channels/sensorIndices
%        If CHAN = 'eegmeg', read eeg and meg channels/sensorIndices
%
% TIME - eg. [0 5] - the desired time interval to read, in sec.
%        If TIME = 'all', all data is read (the default)
%
% TRIALS - If TRIALS = n, the nth trial will be read.
%          If TRIALS = [3,5,8], reads trials 3,5, and 8 such that
%          ctf.data(:,:,1) = data for trial 3,
%          ctf.data(:,:,2) = data for trial 5, and 
%          ctf.data(:,:,3) = data for trial 8.
%          If TRIALS = [3:7], reads trials 3 to 7
%          If TRIALS = 'all', reads all data (the default)
%
% OUTPUTS
%
% ctf.data - matrix of all the data read, such that data(x,y,z)
%            contains sample point x, channel y and trial z.  The
%            input options, CHAN, TIME, TRIALS can be used to 
%            select subsections of the .meg4 data matrix
%
% ctf.sensor - has the following fields:
%              .names       - cell array of sensor names
%              .location    - array of sensor locations for plotting
%              .orientation - array of sensor orientations
%
% ctf.res4 - has the following fields
%            .file   - the .res4 file path and file name
%            .header - the format of the .res4 file
%
% ctf.meg4 - has the following fields
%            .file   - the .meg4 file path and file name
%            .header - the format of the .meg4 file
%
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %
%      <                                                       > %  
%      <                      DISCLAIMER:                      > %
%      <                                                       > %
%      < THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY.  > %
%      < THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR    > %
%      <                     OFFICIAL USE.                     > %
%      <                                                       > %
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<> %
%


% $Revision: 1.1 $ $Date: 2009-01-30 03:49:27 $

% Copyright (C) 2003  Darren L. Weber
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Modified: 11/2003, Darren.Weber_at_radiology.ucsf.edu
%                    - modified from NIH code simply to allocate data into
%                    one large struct (ctf)
%                    - modified channel selection section at the end so
%                    that it doesn't try to get orientation information for
%                    EEG channels
%                    - changed ctf.data into a 3D matrix, rather than a
%                    cell array of matrices
% Modified: 07/2004, Arnaud Delorme, fixed reading time interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------
% check the function input parameters and assign any defaults

if ~exist('CHAN','var'),   CHAN   = 'all'; end
if ~exist('TIME','var'),   TIME   = 'all'; end
if ~exist('TRIALS','var'), TRIALS = 'all'; end

if isempty(CHAN),   CHAN   = 'all'; end
if isempty(TIME),   TIME   = 'all'; end
if isempty(TRIALS), TRIALS = 'all'; end

CHAN = ctf_channel_select(ctf,CHAN);

%-----------------------------------------------
% ensure we have the data parameters

if ~exist('COEFS','var'),
    COEFS = false;
end

if ~exist('folder','var'),
  if ~exist('ctf','var'),
    ctf = ctf_folder;
  else
    ctf = ctf_folder([],ctf);
  end
else
  if ~exist('ctf','var'),
    ctf = ctf_folder(folder);
  else
    ctf = ctf_folder(folder,ctf);
  end
end

if ~isfield(ctf,'setup'),
  ctf = ctf_read_res4(ctf.folder,1,COEFS);
end




%--------------------------------------------------------------
ver = '$Revision: 1.1 $';
fprintf('\nCTF_READ_MEG4 [v %s]\n',ver(11:15)); tic;


%----------------------------------------------------------------
% open the data file

[folderPath,folderName,folderExt] = fileparts(ctf.folder);
ctf.meg4.file = findmeg4file( ctf.folder );
[fid,message] = fopen(ctf.meg4.file,'rb','s');
if fid < 0, error('cannot open .meg4 file'); end


%----------------------------------------------------------------
% Read the header

% The data file consists of a header and the raw samples from the
% electronics. The header is the 8-byte character sequence: MEG41CP+NULL.

header_bytes = 8;

ctf.meg4.header = char(fread(fid,[1,header_bytes],'char'));

% check the format
if strmatch('MEG41CP',ctf.meg4.header),
  % OK, we can handle this format
else
  msg = sprintf('May not read "%s" format correctly',ctf.meg4.header);
  warning(msg);
end



%-------------------------------------------------------------------
% double check the input parameters

switch num2str(TIME),
  case 'all',
    TIME = ctf.setup.time_sec;
    TIME_index = 1:ctf.setup.number_samples;
  otherwise
    % assume the input is a range of times in sec
    % check the range
    if TIME(1) < ctf.setup.time_sec(1),
      fprintf('...setting TIME(1) = ctf.setup.time_sec(1)\n');
      TIME(1) = ctf.setup.time_sec(1);
    end
    if TIME(end) > ctf.setup.time_sec(end),
      fprintf('...setting TIME(end) = ctf.setup.time_sec(end)\n');
      TIME(end) = ctf.setup.time_sec(end);
    end
    % now find the nearest indices into the samples matrix
    TIME_index = intersect(find(ctf.setup.time_sec >= TIME(1)), find(ctf.setup.time_sec <= TIME(end)));
    % TIME_index = interp1(ctf.setup.time_sec,1:ctf.setup.number_samples,TIME,'nearest');
    % now ensure that the TIME array is consistent with ctf.setup.time_sec
    TIME = ctf.setup.time_sec(TIME_index);
end
TIME = sort(TIME);

% check the duration
duration = TIME(end) - TIME(1);
if duration > ctf.setup.duration_trial,
  fprintf('...TIME input too large for trial\n');
  fprintf('...setting TIME = %g seconds (ctf.setup.duration_trial)',ctf.setup.duration_trial);
  duration = ctf.setup.duration_trial;
end
if duration <= 0,
  fprintf('...TIME(end) - TIME(1) is <= 0, quitting now!\n');
  return
end


% calculate the number of samples selected
number_samples = round((duration) * ctf.setup.sample_rate) + 1;



switch num2str(TRIALS),
  case 'all',
    TRIALS = 1:ctf.setup.number_trials;
  otherwise
    % assume the input is an array of trials
end
TRIALS = unique(sort(TRIALS));




%----------------------------------------------------------------
% Calculate sensor gains

megIndex =   ctf.sensor.index.meg_sens;
refIndex =   ctf.sensor.index.meg_ref;
eegIndex =   ctf.sensor.index.eeg_sens;
otherIndex = ctf.sensor.index.other;

channel_gain = zeros(1,ctf.setup.number_channels);

channel_gain(megIndex)   = [ctf.sensor.info(megIndex).proper_gain] .* [ctf.sensor.info(megIndex).q_gain];
channel_gain(refIndex)   = [ctf.sensor.info(refIndex).proper_gain] .* [ctf.sensor.info(refIndex).q_gain];
channel_gain(eegIndex)   = [ctf.sensor.info(eegIndex).q_gain];
channel_gain(otherIndex) = [ctf.sensor.info(otherIndex).q_gain];


%-------------------------------------------------------------------------
% Read trial data from .meg4 file

% The data is stored as a sequence of (signed) 4-byte integers, starting
% with the first trial and first channel, then the first trial and second
% channel, etc. The number of channels per trial and the number of samples
% in every trialchannel block are constant per dataset. The constants are
% found in the general resources stored in the resource file, see 'The
% Resource File Format'. The numbers stored in the data file are the raw
% numbers collected from the electronics. For these numbers to be useful,
% the various gains, stored in the sensor resources, must be applied.

number_trials = length(TRIALS);
number_channels = length(CHAN);

dataSize = [number_samples, number_channels, number_trials];
v = version('-release');
if strmatch(v, {'11','12','13'}),
  ctf.data = zeros(dataSize);
else,
  ctf.data = repmat(single(0), dataSize);
end

% Calculate trial byte size
trial_bytes = 4 * ctf.setup.number_samples * ctf.setup.number_channels;

trial_count = 0;
for trial = 1:length(TRIALS),
  
  trial_number = TRIALS(trial);
  
  if trial > 1,
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
  end
  fprintf('...reading %4d of %4d trials\n', trial, ctf.setup.number_trials);
  
  % calculate the byte offset in the file for this trial
  trial_offset = header_bytes + ( trial_bytes * ( trial_number - 1 ) );
  
  for channel = 1:length(CHAN),
    
    channel_number = CHAN(channel);
    
    % calculate the channel offset in the current trial
    channel_offset = 4 * ctf.setup.number_samples * ( channel_number - 1 );
    
    % seek to the trial offset, relative to the beginning of the file
    fseek(fid,trial_offset,-1);
    
    % now seek to the channel offset, in the current trial
    fseek(fid,channel_offset,0);
    
    % read the entire set of samples for this channel
    channel_samples = fread(fid, [ctf.setup.number_samples,1], 'int32');
    
    % extract just the selected time array
    channel_samples = channel_samples(TIME_index);
    
    if channel_gain(channel_number),
        channel_samples2tesla = channel_samples ./ channel_gain(channel_number);
    else
        channel_samples2tesla = channel_samples;
    end
    
    % assign the selected time samples into the ctf.data matrix
    ctf.data(:,channel,trial) = channel_samples2tesla;
    
  end
end

fclose(fid);










%-------------------------------------------------------------------------
% assign sensor locations and orientations for selected channels, this
% section will simplify the data allocated by ctf_read_res4

fprintf('...sorting %d from %d sensors\n',number_channels, ctf.setup.number_channels);

ctf.sensor.location = zeros(3,number_channels);
ctf.sensor.orientation = zeros(3,number_channels);

ctf.sensor.label = [];
ctf.sensor.location = [];
ctf.sensor.orientation = [];

for c = 1:length(CHAN),
    
    channel = CHAN(c);
    
    % All channels have a label
    ctf.sensor.label{1,c} = ctf.sensor.info(channel).label;
    
    % All channels have a location
    
    % EEG channels do not have any orientation
    
    switch ctf.sensor.info(channel).index,
        
        case {ctf.sensor.type.meg_sens, ctf.sensor.type.meg_ref},
            %0=Reference Channels, 
            %1=More Reference Channels, 
            %5=MEG Channels
            
            % MEG channels are radial gradiometers, so they have an inner (1) and
            % an outer (2) location - it might be better to take the average of
            % their locations
            if ~isempty(ctf.sensor.info(channel).location),
                ctf.sensor.location(:,c) = ctf.sensor.info(channel).location(:,1);
            end
            
            if ~isempty(ctf.sensor.info(channel).orientation),
                ctf.sensor.orientation(:,c) = ctf.sensor.info(channel).orientation(:,1);
            end
            
        case ctf.sensor.type.eeg_sens,
            %9=EEG Channels
            
            if ~isempty(ctf.sensor.info(channel).location),
                ctf.sensor.location(:,c) = ctf.sensor.info(channel).location(:,1);
            end
            
    end
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEED TO CHECK ctf.setup parameters here, to adjust for any changes
% required by the CHAN, TIME, TRIALS inputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% modify the setup parameters so they correspond with the data selected
ctf.setup.number_samples  = number_samples;
ctf.setup.number_channels = number_channels;
ctf.setup.number_trials   = number_trials;

if ctf.setup.number_samples ~= size(ctf.data,1),
  error('ctf.setup.number_samples ~= size(ctf.data,1)');
end
if ctf.setup.number_channels ~= size(ctf.data,2),
  error('ctf.setup.number_channels ~= size(ctf.data,2)');
end
if ctf.setup.number_trials ~= size(ctf.data,3),
  error('ctf.setup.number_trials ~= size(ctf.data,3)');
end


t = toc; fprintf('...done (%6.2f sec)\n\n',t);

return



% find file name if truncated or with uppercase extension
% added by Arnaud Delorme June 15, 2004
% -------------------------------------------------------
function meg4name = findmeg4file( folder )

meg4name = dir([ folder filesep '*.meg4' ]);
if isempty(meg4name)
    meg4name = dir([ folder filesep '*.MEG4' ]);
end;

if isempty(meg4name)
    error('No file with extension .meg4 or .MEG4 in selected folder');
else
    meg4name = [ folder filesep meg4name.name ];
end;

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------
function sensorName = parse_sensor_label(temp)

% sensorName = parse_sensor_label(temp)
% parse sensor label names

temp(temp>127) = 0;
temp(temp<0) = 0;
temp = strtok(temp,char(0));
temp = strtok(temp,'-');
sensorName = char(temp)';

return
