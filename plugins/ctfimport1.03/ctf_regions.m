function [ctf] = ctf_regions(ctf);

% ctf_channel_sets - Define CTF MEG sensor regions
%
% ctf = ctf_read_res4(ctf)
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

% $Revision: 1.1 $ $Date: 2009-01-30 03:49:27 $

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

ver = '$Revision: 1.1 $';
fprintf('\nCTF_CHANNEL_SETS [v %s]\n',ver(11:15)); tic;

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


% find the indices of the MEG reference sensors

ctf.sensor.type.meg_ref = [0 1];
ctf.sensor.type.refMagnetometers = 0; % include B*
ctf.sensor.type.refGradiometers  = 1; % include G*,P*,Q*,R*

magref  = find([ctf.sensor.info.index] == 0);
gradref = find([ctf.sensor.info.index] == 1);

ctf.sensor.index.ref_mag  = magref;
ctf.sensor.index.ref_grad = gradref;
ctf.sensor.index.meg_ref  = [magref,gradref];


% find the indices of the MEG head sensors

ctf.sensor.type.meg_sens = 5; % include B*
ctf.sensor.index.meg_sens = find([ctf.sensor.info.index] == 5);


% find the indices of the EEG head sensors

ctf.sensor.type.eeg_sens = 9; % include EEG*
ctf.sensor.index.eeg_sens = find([ctf.sensor.info.index] == 9);


% find the indices of the STIM channels

ctf.sensor.type.stim_ref = 11; % include STIM*
ctf.sensor.index.stim_ref = find([ctf.sensor.info.index] == 11);


% find the indices of the system clock

ctf.sensor.type.sclk_ref = 17; % include SCLK*

ctf.sensor.index.sclk_ref = find([ctf.sensor.info.index] == 17);


% find the indices of virtual channels
%ctf.sensor.type.vc  = ??;
%vcsens = find([ctf.sensor.info.index] == ??);


ctf.sensor.index.all_sens = [ ctf.sensor.index.meg_sens, ctf.sensor.index.eeg_sens ];

all = [ ctf.sensor.index.all_sens, ctf.sensor.index.meg_ref, ...
        ctf.sensor.index.stim_ref, ctf.sensor.index.sclk_ref ];

ctf.sensor.index.other = setdiff([1:ctf.setup.number_channels],all);


% define regions of MEG sensors


ctf.sensor.index.meg_left_central    = [];
ctf.sensor.index.meg_left_frontal    = [];
ctf.sensor.index.meg_left_occipital  = [];
ctf.sensor.index.meg_left_parietal   = [];
ctf.sensor.index.meg_left_temporal   = [];
ctf.sensor.index.meg_right_central   = [];
ctf.sensor.index.meg_right_frontal   = [];
ctf.sensor.index.meg_right_occipital = [];
ctf.sensor.index.meg_right_parietal  = [];
ctf.sensor.index.meg_right_temporal  = [];
ctf.sensor.index.meg_mid_central     = [];
ctf.sensor.index.meg_mid_frontal     = [];
ctf.sensor.index.meg_mid_occipital   = [];
ctf.sensor.index.meg_mid_parietal    = [];


for megIndex = ctf.sensor.index.meg_sens,
  megLabel = ctf.sensor.info(megIndex).label;
  if findstr('MLC',megLabel),
    % these are MEG left central
    ctf.sensor.index.meg_left_central(end+1) = megIndex;
  end
  if findstr('MLF',megLabel),
    % these are MEG left frontal
    ctf.sensor.index.meg_left_frontal(end+1) = megIndex;
  end
  if findstr('MLO',megLabel),
    % these are MEG left occipital
    ctf.sensor.index.meg_left_occipital(end+1) = megIndex;
  end
  if findstr('MLP',megLabel),
    % these are MEG left parietal
    ctf.sensor.index.meg_left_parietal(end+1) = megIndex;
  end
  if findstr('MLT',megLabel),
    % these are MEG left temporal
    ctf.sensor.index.meg_left_temporal(end+1) = megIndex;
  end
  if findstr('MRC',megLabel),
    % these are MEG right central
    ctf.sensor.index.meg_right_central(end+1) = megIndex;
  end
  if findstr('MRF',megLabel),
    % these are MEG right frontal
    ctf.sensor.index.meg_right_frontal(end+1) = megIndex;
  end
  if findstr('MRO',megLabel),
    % these are MEG right occipital
    ctf.sensor.index.meg_right_occipital(end+1) = megIndex;
  end
  if findstr('MRP',megLabel),
    % these are MEG right parietal
    ctf.sensor.index.meg_right_parietal(end+1) = megIndex;
  end
  if findstr('MRT',megLabel),
    % these are MEG right temporal
    ctf.sensor.index.meg_right_temporal(end+1) = megIndex;
  end
  if findstr('MZC',megLabel),
    % these are MEG mid central
    ctf.sensor.index.meg_mid_central(end+1) = megIndex;
  end
  if findstr('MZF',megLabel),
    % these are MEG mid frontal
    ctf.sensor.index.meg_mid_frontal(end+1) = megIndex;
  end
  if findstr('MZO',megLabel),
    % these are MEG mid occipital
    ctf.sensor.index.meg_mid_occipital(end+1) = megIndex;
  end
  if findstr('MZP',megLabel),
    % these are MEG mid parietal
    ctf.sensor.index.meg_mid_parietal(end+1) = megIndex;
  end
end

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


t = toc; fprintf('...done (%6.2f sec)\n\n',t);

return
