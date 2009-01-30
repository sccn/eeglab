function [ctf] = ctf_channel_sets(ctf);

% ctf_channel_sets - Define CTF MEG sensor regions
%
% ctf = ctf_channel_sets(ctf)
% 
% This function parses the resource information in ctf.sensor,
% which is returned from ctf_read_res4
%
% INPUTS
%
%   ctf - the struct created by ctf_read_res4
%
% OUTPUTS
%
%   ctf.sensor.index - various regional sensor indices
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

% $Revision: 1.1 $ $Date: 2009-01-30 03:49:26 $

% Copyright (C) 2004  Darren L. Weber
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

% Modified: 02/2004, Darren.Weber_at_radiology.ucsf.edu
%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ver = '$Revision: 1.1 $';
%fprintf('\nCTF_CHANNEL_SETS [v %s]\n',ver(11:15)); tic;

if ~exist('ctf','var'),
  ctf = ctf_folder;
end

if ~isfield(ctf,'setup'),
  ctf = ctf_read_res4(ctf.folder);
end

%-------------------------------------------------------------
% Find channel types and define channel sets, see the
% System Administrators .pdf, 'Channel Sets Configuration'

%types = unique([ctf.sensor.info.index])
% 0     1     5     9    11    17
%for typeValue = types,
%  index  = find([ctf.sensor.info.index] == typeValue);
%  labels = ctf.sensor.info(index).label;
%end



% ctf.SensorNames = {...
%     'Ref Magnetometer',... % Sensor Type Index of 0
%     'Ref Gradiometer' ,... % Index of 1
%     ''  ,...               % 2
%     '' ,...                % 3
%     '' ,...                % 4
%     'MEG Sensor',...       % 5
%     '' ,...                % 6
%     '',...                 % 7
%     '',...                 % 8
%     'EEG Sensor',...       % 9
%     'ADC Input Current',...% 10 ADC Input Current (Amps)
%     'Stimulation input',...% 11
%     'Video Time',...       % 12
%     '',...                 % 13
%     '',...                 % 14
%     'SAM Sensor',...       % 15
%     'Virtual Channel',...  % 16
%     'System Clock',...     % 17 System Time Ref
%     'ADC Input Voltage',...% 18 ADC Input Voltage (Volts)
%   };



% find the indices of the MEG reference sensors

ctf.sensor.type.meg_ref = [0 1];
ctf.sensor.type.refMagnetometers = 0; % include B*
ctf.sensor.type.refGradiometers  = 1; % include G*,P*,Q*,R*

magref  = find([ctf.sensor.info.index] == 0);
gradref = find([ctf.sensor.info.index] == 1);

ctf.sensor.index.meg_ref_mag  = magref;
ctf.sensor.index.meg_ref_grad = gradref;
ctf.sensor.index.meg_ref  = [magref,gradref];


% find the indices of the MEG head sensors

ctf.sensor.type.meg_sens = 5; % include B*
ctf.sensor.index.meg_sens = find([ctf.sensor.info.index] == 5);


% find the indices of the EEG head sensors

ctf.sensor.type.eeg_sens = 9; % include EEG*
ctf.sensor.index.eeg_sens = find([ctf.sensor.info.index] == 9);


% find the indices of the ADC channels

ctf.sensor.type.adc = 10;
ctf.sensor.index.adc = find([ctf.sensor.info.index] == 10);


% find the indices of the STIM channels

ctf.sensor.type.stim_ref = 11;
ctf.sensor.index.stim_ref = find([ctf.sensor.info.index] == 11);


% find the indices of the STIM channels

ctf.sensor.type.video_time = 12;
ctf.sensor.index.video_time = find([ctf.sensor.info.index] == 12);


% find the indices of the SAM channels (Synthetic Aperture Magnetometry)

ctf.sensor.type.sam = 15;
ctf.sensor.index.sam = find([ctf.sensor.info.index] == 15);


% find the indices of the virtual channels

ctf.sensor.type.virtual_channels = 16;
ctf.sensor.index.virtual_channels = find([ctf.sensor.info.index] == 16);


% find the indices of the system clock

ctf.sensor.type.sclk_ref = 17; % include SCLK*
ctf.sensor.index.sclk_ref = find([ctf.sensor.info.index] == 17);


% combine eeg and meg channels

ctf.sensor.index.all_sens = [ ctf.sensor.index.meg_sens, ctf.sensor.index.eeg_sens ];

all = [ ctf.sensor.index.all_sens, ctf.sensor.index.meg_ref, ...
        ctf.sensor.index.stim_ref, ctf.sensor.index.sclk_ref ];

ctf.sensor.index.other = setdiff([1:ctf.setup.number_channels],all);


% define regions of MEG sensors
labels = {ctf.sensor.info.label};
ctf.sensor.index.meg_left_central    = strmatch('MLC',labels)';
ctf.sensor.index.meg_left_frontal    = strmatch('MLF',labels)';
ctf.sensor.index.meg_left_parietal   = strmatch('MLP',labels)';
ctf.sensor.index.meg_left_occipital  = strmatch('MLO',labels)';
ctf.sensor.index.meg_left_parietal   = strmatch('MLP',labels)';
ctf.sensor.index.meg_left_temporal   = strmatch('MLT',labels)';
ctf.sensor.index.meg_right_central   = strmatch('MRC',labels)';
ctf.sensor.index.meg_right_frontal   = strmatch('MRF',labels)';
ctf.sensor.index.meg_right_occipital = strmatch('MRO',labels)';
ctf.sensor.index.meg_right_parietal  = strmatch('MRP',labels)';
ctf.sensor.index.meg_right_temporal  = strmatch('MRT',labels)';
ctf.sensor.index.meg_mid_central     = strmatch('MZC',labels)';
ctf.sensor.index.meg_mid_frontal     = strmatch('MZF',labels)';
ctf.sensor.index.meg_mid_occipital   = strmatch('MZO',labels)';
ctf.sensor.index.meg_mid_parietal    = strmatch('MZP',labels)';

ctf.sensor.index.meg_left  = [ ...
    ctf.sensor.index.meg_left_central, ...
    ctf.sensor.index.meg_left_frontal, ...
    ctf.sensor.index.meg_left_occipital, ...
    ctf.sensor.index.meg_left_parietal, ...
    ctf.sensor.index.meg_left_temporal ];

ctf.sensor.index.meg_right  = [ ...
    ctf.sensor.index.meg_right_central, ...
    ctf.sensor.index.meg_right_frontal, ...
    ctf.sensor.index.meg_right_occipital, ...
    ctf.sensor.index.meg_right_parietal, ...
    ctf.sensor.index.meg_right_temporal ];

ctf.sensor.index.meg_central = [ ...
    ctf.sensor.index.meg_left_central, ...
    ctf.sensor.index.meg_right_central, ...
    ctf.sensor.index.meg_mid_central ];

ctf.sensor.index.meg_frontal = [ ...
    ctf.sensor.index.meg_left_frontal, ...
    ctf.sensor.index.meg_right_frontal, ...
    ctf.sensor.index.meg_mid_frontal ];

ctf.sensor.index.meg_occipital = [ ...
    ctf.sensor.index.meg_left_occipital, ...
    ctf.sensor.index.meg_right_occipital, ...
    ctf.sensor.index.meg_mid_occipital ];

ctf.sensor.index.meg_parietal = [ ...
    ctf.sensor.index.meg_left_parietal, ...
    ctf.sensor.index.meg_right_parietal, ...
    ctf.sensor.index.meg_mid_parietal ];

ctf.sensor.index.meg_temporal = [ ...
    ctf.sensor.index.meg_left_temporal, ...
    ctf.sensor.index.meg_right_temporal ];


%t = toc; fprintf('...done (%6.2f sec)\n\n',t);

return
