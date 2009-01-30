function [CHAN,type] = ctf_channel_select(ctf,CHAN)

% ctf_channel_select - select channels using character codes
% 
% [CHAN,type] = ctf_channel_select(ctf,CHAN)
%
% where CHAN input is an array of channel numbers or a string that
% corresponds to a channel set below (these channel sets are defined by
% ctf_channel_sets):
%
% 'all', CHAN = [1:ctf.setup.number_channels];  % default
% 'ref', CHAN = ctf.sensor.index.meg_ref;
% 'meg', CHAN = ctf.sensor.index.meg_sens;
% 'eeg', CHAN = ctf.sensor.index.eeg_sens;
% 'other', CHAN = ctf.sensor.index.other;
% {'megeeg','eegmeg'}, CHAN = [ ctf.sensor.index.eeg_sens ctf.sensor.index.meg_sens ];
% 'lc', CHAN = ctf.sensor.index.meg_left_central;
% 'lf', CHAN = ctf.sensor.index.meg_left_frontal;
% 'lo', CHAN = ctf.sensor.index.meg_left_occipital;
% 'lp', CHAN = ctf.sensor.index.meg_left_parietal;
% 'lt', CHAN = ctf.sensor.index.meg_left_temporal;
% 'rc', CHAN = ctf.sensor.index.meg_right_central;
% 'rf', CHAN = ctf.sensor.index.meg_right_frontal;
% 'ro', CHAN = ctf.sensor.index.meg_right_occipital;
% 'rp', CHAN = ctf.sensor.index.meg_right_parietal;
% 'rt', CHAN = ctf.sensor.index.meg_right_temporal;
% 'mc', CHAN = ctf.sensor.index.meg_mid_central;
% 'mf', CHAN = ctf.sensor.index.meg_mid_frontal;
% 'mo', CHAN = ctf.sensor.index.meg_mid_occipital;
% 'mp', CHAN = ctf.sensor.index.meg_mid_parietal;
% 'l',  CHAN = ctf.sensor.index.meg_left;
% 'r',  CHAN = ctf.sensor.index.meg_right;
% 'c',  CHAN = ctf.sensor.index.meg_central;
% 'f',  CHAN = ctf.sensor.index.meg_frontal;
% 'o',  CHAN = ctf.sensor.index.meg_occipital;
% 'p',  CHAN = ctf.sensor.index.meg_parietal;
% 't',  CHAN = ctf.sensor.index.meg_temporal;
%
% CHAN output is a numeric array of channel indices
% type output is the string type of the channel indices
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
%           02/2005, removed unique(CHAN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ver = '$Revision: 1.1 $ $Date: 2009-01-30 03:49:26 $';
fprintf('\nCTF_CHANNEL_SELECT [v %s]\n',ver(11:15));


if ~exist('CHAN','var'), CHAN = 'all'; end
if ~exist('BAD','var'), BAD = 1; end

if ischar(CHAN),
    fprintf('...selecting %s channels\n', CHAN);
else
    fprintf('...selecting %d channels\n', length(CHAN));
end

type = num2str(lower(CHAN));

switch type,
 case 'ref',
  CHAN = ctf.sensor.index.meg_ref;
 case 'meg',
  CHAN = ctf.sensor.index.meg_sens;
 case 'eeg',
  CHAN = ctf.sensor.index.eeg_sens;
 case 'other',
  CHAN = ctf.sensor.index.other;
 case {'megeeg','eegmeg'},
  CHAN = [ ctf.sensor.index.eeg_sens ctf.sensor.index.meg_sens ];
 case 'all',
  CHAN = [1:ctf.setup.number_channels];
 case 'lc',
  CHAN = ctf.sensor.index.meg_left_central;
 case 'lf',
  CHAN = ctf.sensor.index.meg_left_frontal;
 case 'lo',
  CHAN = ctf.sensor.index.meg_left_occipital;
 case 'lp',
  CHAN = ctf.sensor.index.meg_left_parietal;
 case 'lt',
  CHAN = ctf.sensor.index.meg_left_temporal;
 case 'rc',
  CHAN = ctf.sensor.index.meg_right_central;
 case 'rf',
  CHAN = ctf.sensor.index.meg_right_frontal;
 case 'ro',
  CHAN = ctf.sensor.index.meg_right_occipital;
 case 'rp',
  CHAN = ctf.sensor.index.meg_right_parietal;
 case 'rt',
  CHAN = ctf.sensor.index.meg_right_temporal;
 case 'mc',
  CHAN = ctf.sensor.index.meg_mid_central;
 case 'mf',
  CHAN = ctf.sensor.index.meg_mid_frontal;
 case 'mo',
  CHAN = ctf.sensor.index.meg_mid_occipital;
 case 'mp',
  CHAN = ctf.sensor.index.meg_mid_parietal;
 case 'l',
  CHAN = ctf.sensor.index.meg_left;
 case 'r',
  CHAN = ctf.sensor.index.meg_right;
 case 'c',
  CHAN = ctf.sensor.index.meg_central;
 case 'f',
  CHAN = ctf.sensor.index.meg_frontal;
 case 'o',
  CHAN = ctf.sensor.index.meg_occipital;
 case 'p',
  CHAN = ctf.sensor.index.meg_parietal;
 case 't',
  CHAN = ctf.sensor.index.meg_temporal;
 otherwise
  type = 'numeric';
  % assume the input is an array of channel numbers
end

%CHAN = unique(CHAN);


% remove bad channels
% if BAD,
%     fprintf('...removing bad channels\n');
%     ctf = ctf_read_badchannels([],ctf);
%     badChanNames = ctf.sensor.bad;
%     chanNames = ctf.sensor.label(CHAN);
%     [newChanNames,index] = setdiff(chanNames,badChanNames);
%     index = sort(index); % to undo the sorting of setdiff
%     CHAN = CHAN(index);
% end


return
