function [ctf] = ctf_read_res4(folder,VERBOSE,COEFS);
  
% ctf_read_res4 - Read a CTF .res4 file
%
% ctf = ctf_read_res4( [folder], [verbose], [coefs])
% 
% This function reads the resource information from a CTF .ds folder.  This
% resource information must be read before reading the .meg4 data file.
% All input arguments are optional.
%
% INPUTS
%
%   folder: the .ds directory containing the data to be read.  With
%           no input, it will prompt with a gui folder locator 
%           (called by ctf_folder).
%
%   verbose: If verbose = 1, display 'ctf.setup' structure (default)
%            If verbose = 0, do not display 'ctf.setup' structure
%
%   coefs: an option to read the MEG sensor and reference coefficients, 
%          which give the weights for calculation of synthetic 2nd or 3rd 
%          order gradiometers.
%          If coefs = 1, read the sensor coefficients
%          If coefs = 0, do not read the sensor coefficients (this is the
%          default because it is assumed that data preprocessing with
%          CTF tools has already applied AND saved a synthetic
%          gradiometer transformation for the meg4 data file).  The
%          ctf.sensor.info(x).grad_order_no will indicate the gradient 
%          order of the data (this only applies to MEG sensors, other
%          channels have zero values).
% 
% OUTPUTS
%
%   ctf.folder - path of the .ds folder
%
%   ctf.res4.file - data file path/name
%   ctf.res4.header - data format header
%
%   ctf.setup - a header structure consisting of date, time, run name, run
%   title, subject, run description, operator, number of channels, number
%   of samples, sample rate, number of trials, duration, pretrigger_samples,
%   sensor filename, head zeroing, and number of filters.
%   
%   ctf.sensor.index - a sensor structure consisting of EEG sensors, MEG
%   sensors, reference sensors, and other sensors.
%   
%   ctf.sensor.info - a structure with gain and offset information consisting
%   of proper gain, Q gain, io gain, io offset, and index.
% 
% The proper gain is the channel-wide constant quotient of a value in raw
% units and a value in user units. (For EEG this is one because there is
% no conversion; i.e., the values are in volts, for MEG there is conversion
% between phi0 and teslas).
% proper gain = raw value / user value
% 
% The Q-gain is used to convert the 32-bit integers produced by the
% electronics to the raw values stored in the data file.
% raw value = 32 bit value / Q gain
% 
% The IO-gain is essentially a correction factor that is applied to the
% raw values before they are stored to disk.
% 
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%      <                                                      >
%      <                    DISCLAIMER:                       >
%      <                                                      >
%      < THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. >
%      < THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   >
%      <                    OFFICIAL USE.                     >
%      <                                                      >
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>
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
%                    - modified from NIH code readresfile.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $';
fprintf('\nCTF_READ_RES4 [v %s]\n',ver(11:15)); tic;

if ~exist('folder','var'),
    ctf = ctf_folder;
else
    ctf = ctf_folder(folder);
end
%[folderPath,folderName,folderExt] = fileparts(ctf.folder);
ctf.res4.file = findres4file(ctf.folder);

if ~exist('VERBOSE','var'), VERBOSE = 1; end
if ~exist('COEFS','var'), COEFS = 0; end

%----------------------------------------------------------------
% open the data file

[fid,message] = fopen(ctf.res4.file,'rb','ieee-be.l64');
if fid < 0, error('cannot open .res4 file'); end


%-------------------------------------------------------------
% READ HEADER

fseek(fid,0,-1);
ctf.res4.header = char(fread(fid,8,'char'))';

% check for the right format
if strmatch('MEG41RS',ctf.res4.header),
    % OK, we can handle this format
else
    msg = sprintf('This function is designed to read MEG41RS format.\nIt may not read "%s" format correctly',ctf.res4.header);
    warning(msg);
end


%-------------------------------------------------------------
% READ SETUP


%---DATE/TIME
fseek(fid, 778,-1);
ctf.setup.time = char(fread(fid,255,'char'))';
fseek(fid,1033,-1);
ctf.setup.date = char(fread(fid,255,'char'))';

%---NUMBER OF SAMPLES
fseek(fid,1288,-1);
ctf.setup.number_samples = (fread(fid,1,'int32')');

%---NUMBER OF CHANNELS
fseek(fid,1292,-1);
ctf.setup.number_channels = (fread(fid,1,'int16')');

%---SAMPLE RATE
fseek(fid,1296,-1);
ctf.setup.sample_rate = fread(fid,1,'double')';
ctf.setup.sample_msec = 1000 / ctf.setup.sample_rate;
ctf.setup.sample_sec  =    1 / ctf.setup.sample_rate;

%---NUMBER OF TRIALS
fseek(fid,1312,-1);
ctf.setup.number_trials = fread(fid,1,'int16')';
fseek(fid, 776,-1);
ctf.setup.number_trials_averaged = fread(fid,1,'int16');

%---DURATION
fseek(fid,1304,-1);
ctf.setup.duration_total = fread(fid,1,'double');
ctf.setup.duration_trial = ctf.setup.duration_total / ctf.setup.number_trials;

%---PRE_TRIG POINTS
fseek(fid,1316,-1);
ctf.setup.pretrigger_samples = fread(fid,1,'int32');
ctf.setup.pretrigger_msec = (ctf.setup.pretrigger_samples / ctf.setup.sample_rate) * 1000;

%---HEAD ZEROING
fseek(fid,1348,-1);
h_zero = fread(fid,1,'int32')';
no_yes = {'no','yes'};
ctf.setup.head_zero = no_yes{h_zero+1};

%---RUN NAME
fseek(fid,1360,-1);
ctf.setup.run_name = char(fread(fid,32,'char'))';

%---RUN TITLE
fseek(fid,1392,-1);
ctf.setup.run_title = char(fread(fid,256,'char'))';

%---SUBJECT
fseek(fid,1712,-1);
ctf.setup.subject = char(fread(fid,32,'char'))';

%---OPERATOR
fseek(fid,1744,-1);
ctf.setup.operator = char(fread(fid,32,'char'))';

%---SENSOR FILE NAME
fseek(fid,1776,-1);
ctf.setup.sensor_file_name = char(fread(fid,60,'char'))';

%---RUN DESCRIPTION & FILTERS
fseek(fid,1836,-1);
run_size = fread(fid,1,'int32');
fseek(fid,1844,-1);
ctf.setup.run_description = char(fread(fid,run_size,'char'));

ctf.setup.number_filters = fread(fid,1,'int16');

for i = 1:ctf.setup.number_filters,
    ctf.setup.filters(i).freq     = fread(fid,1,'double');
    ctf.setup.filters(i).class    = fread(fid,1,'int32');
    ctf.setup.filters(i).type     = fread(fid,1,'int32');
    ctf.setup.filters(i).numparam = fread(fid,1,'int16');
    ctf.setup.filters(i).params   = fread(fid,ctf.setup.filters(i).numparam,'double');
end

if(COEFS)
    b = 1846 + run_size;
    if(ctf.setup.number_filters == 0)
        np = 0;
    else
        np = sum([ctf.setup.filters.numparam]);
        warning('3rd gradient + hardware filter parameters not fully tested! let''s see what happens... :)');
    end
    nf = ctf.setup.number_filters;
    f = ( nf * 18 ) + ( np * 8 );
    offset = b + f + ctf.setup.number_channels * 1360;
end


%-------------------------------------------------------------
% CREATE TIME ARRAYS

% the time arrays must be based on increments of the sample_msec
ctf.setup.time_msec = [0:ctf.setup.number_samples - 1]' * ctf.setup.sample_msec;
ctf.setup.time_msec = ctf.setup.time_msec - ctf.setup.pretrigger_msec;

% adjust the sample point closest to zero so that it is zero, if it
% is reasonably close to zero, say within 3 decimal places for msec timing
zero_index = find(abs(ctf.setup.time_msec) == min(abs(ctf.setup.time_msec)));
zero_value = ctf.setup.time_msec(zero_index);
if (-0.0001 < zero_value) & (zero_value < 0.0001),
    ctf.setup.time_msec(zero_index) = 0;
end

ctf.setup.start_msec = ctf.setup.time_msec(1);
ctf.setup.end_msec   = ctf.setup.time_msec(end);

ctf.setup.time_sec  = ctf.setup.time_msec / 1000;
ctf.setup.start_sec = ctf.setup.time_sec(1);
ctf.setup.end_sec   = ctf.setup.time_sec(end);


%-------------------------------------------------------------
% PRINT SETUP
if VERBOSE,
    ctf_print_setup(ctf);
end





%-------------------------------------------------------------
% READ SENSOR INFORMATION

ctf.sensor.info = struct(...
    'proper_gain',[],...
    'q_gain',[],...
    'io_gain',[],...
    'io_offset',[],...
    'index',[],...
    'extra',[],...
    'label',[],...
    'grad_order_no',[]);

% read channel names
for chan = 1:ctf.setup.number_channels,
    temp = fread(fid,32,'char');
    ctf.sensor.info(chan).label = parse_sensor_label(temp);
end

for chan = 1:ctf.setup.number_channels,
    
    
    % The proper gain is the channel-wide constant quotient of a value in raw
    % units and a value in user units. (For EEG this is one because there is
    % no conversion; i.e., the values are in volts, for MEG there is conversion
    % between phi0 and teslas).
    % proper gain = raw value / user value
    % 
    % The Q-gain is used to convert the 32-bit integers produced by the
    % electronics to the raw values stored in the data file.
    % raw value = 32 bit value / Q gain
    % 
    % The IO-gain is essentially a correction factor that is applied to the
    % raw values before they are stored to disk.
    
    
    %ftell(fid);
    
    ctf.sensor.info(chan).index = fread(fid,1,'int16');
    ctf.sensor.info(chan).extra = fread(fid,1,'int16');
    id = fread(fid,1,'int32')+1;
    ctf.sensor.info(chan).proper_gain = fread(fid,1,'double');
    ctf.sensor.info(chan).q_gain = fread(fid,1,'double');
    ctf.sensor.info(chan).io_gain = fread(fid,1,'double');
    ctf.sensor.info(chan).io_offset = fread(fid,1,'double');
    fread(fid,1,'int16');
    ctf.sensor.info(chan).grad_order_no = fread(fid,1,'int16');
    fread(fid,1,'int32');
    
    %fseek(fid,ftell(fid)+6,0);
    
    for pos = 1:8,
        ctf.sensor.info(chan).coil(pos).position.x = fread(fid,1,'double');
        ctf.sensor.info(chan).coil(pos).position.y = fread(fid,1,'double');
        ctf.sensor.info(chan).coil(pos).position.z = fread(fid,1,'double');
        fread(fid,1,'double');
        ctf.sensor.info(chan).coil(pos).orient.x = fread(fid,1,'double');
        ctf.sensor.info(chan).coil(pos).orient.y = fread(fid,1,'double');
        ctf.sensor.info(chan).coil(pos).orient.z = fread(fid,1,'double');
        fread(fid,1,'double');
        fread(fid,1,'int16');
        fread(fid,1,'int32');
        fread(fid,1,'int16');
        fread(fid,1,'double');
        
        %fseek(fid,ftell(fid)+56,0);
        %fseek(fid,ftell(fid)-80,0);
    end
    
    for pos = 1:8,
        ctf.sensor.info(chan).hcoil(pos).position.x = fread(fid,1,'double');
        ctf.sensor.info(chan).hcoil(pos).position.y = fread(fid,1,'double');
        ctf.sensor.info(chan).hcoil(pos).position.z = fread(fid,1,'double');
        fread(fid,1,'double');
        ctf.sensor.info(chan).hcoil(pos).orient.x = fread(fid,1,'double');
        ctf.sensor.info(chan).hcoil(pos).orient.y = fread(fid,1,'double');
        ctf.sensor.info(chan).hcoil(pos).orient.z = fread(fid,1,'double');
        fread(fid,1,'double');
        fread(fid,1,'int16');
        fread(fid,1,'int32');
        fread(fid,1,'int16');
        fread(fid,1,'double');
        
        %fseek(fid,ftell(fid)+56,0);
        %fseek(fid,ftell(fid)+80,0);
    end
    %fseek(fid,ftell(fid)+1288,-1);
end




%-------------------------------------------------------------
% Find channel types and define channel sets, see the
% System Administrators .pdf, 'Channel Sets Configuration'

ctf = ctf_channel_sets(ctf);



%-------------------------------------------------------------
% Channel coordinates, in centimeters, in subject head space

for chan = 1:ctf.setup.number_channels,
    
    switch ctf.sensor.info(chan).index,
        
        case {0,1,5},
            %0=Reference Magnetometers
            %1=Reference Gradiometers
            %5=MEG Channels
            
            coord = [ctf.sensor.info(chan).hcoil(1:2).position];
            ctf.sensor.info(chan).location = [coord.x; coord.y; coord.z];
            
            orient = [ctf.sensor.info(chan).hcoil(1).orient];
            ctf.sensor.info(chan).orientation = [orient.x; orient.y; orient.z];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This ensures that the orientation of the sensor is away from the
            % center of a sphere.  It uses the sign of the dot product between
            % the orientation vector and the location vector.
            tmp = ctf.sensor.info(chan).orientation' * ctf.sensor.info(chan).location;
            tmp = sign(tmp(1));
            ctf.sensor.info(chan).orientation = tmp * ctf.sensor.info(chan).orientation;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        case 9,
            %EEG Channels    
            
            coord = [ctf.sensor.info(chan).hcoil(1:2).position];
            ctf.sensor.info(chan).location = [coord.x; coord.y; coord.z];
            ctf.sensor.info(chan).orientation = [];
            
    end
end


%-------------------------------------------------------------------------
% assign sensor locations and orientations
% simplify information in the ctf.sensor.info struct

ctf.sensor.location = zeros(3,ctf.setup.number_channels);
ctf.sensor.orientation = zeros(3,ctf.setup.number_channels);

ctf.sensor.label = [];
ctf.sensor.location = [];
ctf.sensor.orientation = [];

for c = 1:ctf.setup.number_channels,
    
    % All channels have a label
    ctf.sensor.label{1,c} = strtok(ctf.sensor.info(c).label,'-');
    
    % All channels have a location
    % EEG channels do not have any orientation
    if length(ctf.sensor.type.meg_ref) > 1, ind2 = 2; else ind2 = 1; end;
    
    switch ctf.sensor.info(c).index,
        
        case {ctf.sensor.type.meg_sens, ctf.sensor.type.meg_ref(1), ctf.sensor.type.meg_ref(ind2)},
        % modification above based on bug 438 of EEGLAB - Arnaud Delorme, August 2007
            %0=Reference Channels, 
            %1=More Reference Channels, 
            %5=MEG Channels
            
            % MEG channels are radial gradiometers, so they have an inner (1) and
            % an outer (2) location - it might be better to take the average of
            % their locations
            if ~isempty(ctf.sensor.info(c).location),
                ctf.sensor.location(:,c) = ctf.sensor.info(c).location(:,1);
            end
            
            if ~isempty(ctf.sensor.info(c).orientation),
                ctf.sensor.orientation(:,c) = ctf.sensor.info(c).orientation(:,1);
            end
            
        case ctf.sensor.type.eeg_sens,
            %9=EEG Channels
            
            if ~isempty(ctf.sensor.info(c).location),
                ctf.sensor.location(:,c) = ctf.sensor.info(c).location(:,1);
            end
            
    end
end






%%%%%% Coefficient reading appended by SSD %%%%%%

if(COEFS)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NUMBER OF COEFFICIENTS, byte offset = b+f+nc*1360, byte size = 1
    
    fseek(fid,offset,-1);
    ctf.res4.numberCoefficients = fread(fid,1,'int16');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SENSOR COEFFICIENT RECORDS, byte offset = b+f+nc*1360+2, byte size = 1992
    
    if VERBOSE & ctf.res4.numberCoefficients,
        fprintf('...reading %d coefficients\n',ctf.res4.numberCoefficients);
    end
    
    SENSOR_LABEL  = 31;
    MAX_NUM_COEFS = 50;
    MAX_BALANCING = MAX_NUM_COEFS;
    
    hexadef = {'00000000','47314252','47324252','47334252','47324f49','47334f49'};
    strdef = {'NOGRAD','G1BR','G2BR','G3BR','G2OI','G3OI'};
    
    for i = 1:ctf.res4.numberCoefficients,
        
        % read the sensor name (channel label)
        temp = fread(fid,[1,32],'char');
        sensorName = parse_sensor_label(temp);
        sensorIndex = strmatch( sensorName, {ctf.sensor.info.label} );
        
        % read the coefficient type
        coefType = fread(fid,1,'bit32');
        
        padding = fread(fid,1,'int32'); % not sure why this is needed???
        
        % read the coefficient record
        numberCoefs = fread(fid,1,'int16');
        
        if numberCoefs > MAX_NUM_COEFS,
            msg = sprintf('numberCoefs > MAX_NUM_COEFS\n');
            warning(msg);
        end
        
        sensor_list = char(fread(fid,[SENSOR_LABEL,MAX_BALANCING],'uchar')');
        % clean-up the sensor_list
        sensor_list = sensor_list(1:numberCoefs,:);
        for j=1:numberCoefs,
            temp = strtok(sensor_list(j,:),char(0));
            temp = strtok(temp,'-');
            % check if this sensor is a reference
            refLabels = {ctf.sensor.info(ctf.sensor.index.meg_ref).label};
            refIndex = strmatch(temp,refLabels);
            if refIndex,
                % ensure this one has the same label
                temp = refLabels{refIndex};
            end
            
            new_sensor_list(j) = refIndex;
        end
        sensor_list = ctf.sensor.index.meg_ref(:,new_sensor_list)';
        
        coefs_list = fread(fid,MAX_BALANCING,'double');
        % clean-up the coefs_list
        coefs_list = coefs_list(1:numberCoefs,:)';
        
        % allocate the coefficient parameters into the ctf struct
        ctf.res4.sensorCoef(i).sensorName = sensorName;
        ctf.res4.sensorCoef(i).coefType   = coefType;
        ctf.res4.sensorCoef(i).coefRec.numberCoefs = numberCoefs;
        ctf.res4.sensorCoef(i).coefRec.sensor_list = sensor_list;
        ctf.res4.sensorCoef(i).coefRec.coefs_list  = coefs_list;
        
        
        
        % DLW:
        % This is a brainstorm variable, note the use of coefType
        % Not clear why this is checked and allocated as such
        coefType = find( hex2dec(hexadef) == coefType );
        if coefType,
            CoefInfo{sensorIndex,coefType-1}.numberCoefs = numberCoefs;
            CoefInfo{sensorIndex,coefType-1}.sensor_list = sensor_list;
            CoefInfo{sensorIndex,coefType-1}.coefs = coefs_list;
        end
    end	
    
    
    
    % Channel Gains
    gain_chan = zeros(size(ctf.setup.number_channels,1),1);
    gain_chan(ctf.sensor.index.meg_sens) = ([ctf.sensor.info(ctf.sensor.index.meg_sens).proper_gain]'.*[ctf.sensor.info(ctf.sensor.index.meg_sens).q_gain]');
    gain_chan(ctf.sensor.index.meg_ref) = ([ctf.sensor.info(ctf.sensor.index.meg_ref).proper_gain]'.*[ctf.sensor.info(ctf.sensor.index.meg_ref).q_gain]');
    %  gain_chan(ieegsens) = 1./([SensorRes(ieegsens).qGain]'*1e-6);
    %     gain_chan(ieegsens) = 1./([SensorRes(ieegsens).qGain]');
    %     gain_chan(iothersens) = ([SensorRes(iothersens).qGain]'); % Don't know exactly which gain to apply here
    
    
    
    % Calculus of the matrix for nth-order gradient correction
    % Coefficients for unused reference channels are weigthed by zeros in
    % the correction matrix.
    Gcoef = zeros(length(ctf.sensor.index.meg_sens),length(min(ctf.sensor.index.meg_ref):max(ctf.sensor.index.meg_ref)));
    grad_order_no = 3*ones(306,1);
    for k = 1:length(ctf.sensor.index.meg_sens)
        
        % Reference coils for channel k
        if grad_order_no(ctf.sensor.index.meg_sens(k)) == 0
            %Data is saved as RAW
            %Save 3rd order gradient sensor-list for subsequent correction if requested later by the user
            [refs] = (CoefInfo{ctf.sensor.index.meg_sens(k),3}.sensor_list);
            Gcoef(k,refs-min(ctf.sensor.index.meg_ref)+1) = CoefInfo{ctf.sensor.index.meg_sens(k),3}.coefs ... 
                .* gain_chan(refs)'/gain_chan(ctf.sensor.index.meg_sens(k)); 
        else
            [refs] = (CoefInfo{ctf.sensor.index.meg_sens(k),grad_order_no(ctf.sensor.index.meg_sens(k))}.sensor_list);
            Gcoef(k,refs-min(ctf.sensor.index.meg_ref)+1) = CoefInfo{ctf.sensor.index.meg_sens(k),grad_order_no(ctf.sensor.index.meg_sens(k))}.coefs ... 
                .* gain_chan(refs)/gain_chan(ctf.sensor.index.meg_sens(k)); 
        end
        ctf.sensor.info(ctf.sensor.index.meg_sens(k)).Gcoef = Gcoef(k,:);
    end
end   %% end COEF block


fclose(fid);

t = toc; fprintf('...done (%6.2f sec)\n\n',t);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------
function res4name = findres4file( folder )

% find file name if truncated or with uppercase extension
% added by Arnaud Delorme June 15, 2004

res4name = dir([ folder filesep '*.res4' ]);
if isempty(res4name)
    res4name = dir([ folder filesep '*.RES4' ]);
end

if isempty(res4name)
    error('No file with extension .res4 or .RES4 in selected folder');
else
    res4name = [ folder filesep res4name(1).name ];
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
