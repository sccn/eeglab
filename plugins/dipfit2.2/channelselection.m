function [channel] = channelselection(channel, datachannel)

% CHANNELSELECTION for EEG and MEG labels
%
% This function translates the user-specified list of channels into channel
% labels as they occur in the data. This channel selection procedure can
% be used throughout fieldtrip.
%
% You can specify a mixture of real channel labels and of special strings,
% or index numbers that will be replaced by the corresponding channel
% labels. Channels that are not present in the raw datafile are
% automatically removed from the channel list.
%
% E.g.
%  'gui'     a graphical user interface will pop up to select the channels 
%  'all'     is replaced by all channels in the datafile
%  'MEG'     is replaced by all channels in the CTF datafile starting with 'M'
%  'EEG'     is replaced by all channels in the CTF datafile starting with 'EEG'
%  'EEG1020' is replaced by 'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', ...
%  'EOG'     is replaced by all recognized EOG channels
%  'EMG'     is replaced by all channels in the datafile starting with 'EMG'
%  'lfp'     is replaced by all channels in the datafile starting with 'lfp'
%  'mua'     is replaced by all channels in the datafile starting with 'mua'
%  'spike'   is replaced by all channels in the datafile starting with 'spike'
%  10        is replaced by the 10th channel in the datafile
%
% Other channel groups are
%   'EEG1010'    with approximately 90 electrodes
%   'EEG1005'    with approximately 350 electrodes
%   'EEGCHWILLA' for Dorothee Chwilla's electrode caps (used at the NICI)
%   'EEGBHAM'    for the 128 channel EEG system used in Birmingham
%   'EEGREF'     for mastoid and ear electrodes (M1, M2, LM, RM, A1, A2)
%   'MZ'         for MEG central
%   'ML'         for MEG left
%   'MR'         for MEG right
%   'MLx', 'MRx' and 'MZx' with x=C,F,O,P,T for left/right central, frontal,
%   occipital, parietal and temporal
%
% You can also exclude channels or channel groups using the following syntax
%   {'all', '-POz', '-Fp1', -EOG'}

% Copyright (C) 2003-2007, Robert Oostenveld
%
% $Log: channelselection.m,v $
% Revision 1.1  2009/01/30 04:02:01  arno
% *** empty log message ***
%
% Revision 1.19  2007/05/30 13:22:01  roboos
% concatenate indices as column instead of row vector
%
% Revision 1.18  2007/01/22 10:32:28  roboos
% added 'gui' option for graphical user interface, thanks to Vladimir
%
% Revision 1.17  2006/06/06 16:28:54  ingnie
% added option channel is numeric; channel index replaced by channel name
%
% Revision 1.16  2006/05/03 15:08:56  roboos
% already remove double channels prior to the translation of the channel groups
%
% Revision 1.15  2006/03/23 22:26:18  roboos
% added channel groups for lfp, mua and spike
%
% Revision 1.14  2005/12/16 14:03:20  roboos
% added two VEOG channels
% added the EEGBHAM channel group
%
% Revision 1.13  2005/09/14 12:32:24  roboos
% small change in help
%
% Revision 1.12  2005/05/23 09:33:43  roboos
% return immediately if input is empty
%
% Revision 1.11  2004/10/22 16:10:07  roboos
% replaced all occurences of strmatch with strncmp, which runs much faster
%
% Revision 1.10  2004/03/10 15:03:08  roberto
% made selection of EOG channels more general
%
% Revision 1.9  2004/02/24 17:21:34  roberto
% slightly different implementation to undo the channel name sorting
%
% Revision 1.8  2004/02/24 16:53:44  roberto
% added 2 lines that unbdo the alphabetical sorting, the labels are now
% sorted according to their occurence in the data
%
% Revision 1.7  2004/01/26 11:54:47  roberto
% fixed recursion bug in bad channels
% added a line to remove double occurences of channels
%
% Revision 1.6  2004/01/22 21:41:12  roberto
% added support for excluding bad channels or channel groups
%
% Revision 1.5  2004/01/15 17:01:35  roberto
% added channel groups for MZx, with x=f,c,p,o
%
% Revision 1.4  2004/01/09 15:17:33  roberto
% added EEGCHWILLA as channel group (for NICI/Eric Maris)
% added LM and RM to to reference electrodes
%
% Revision 1.3  2003/12/08 12:32:10  roberto
% added 2 channels for eog
%
% Revision 1.2  2003/11/12 07:50:14  roberto
% added separate group for EEG reference electrodes (mastoid, ear)
%
% Revision 1.1  2003/10/28 15:09:27  roberto
% previously known under misc/translate_channel_list.m
% added label to EOG group
%
% Revision 1.2  2003/09/11 21:55:32  roberto
% added labels for 10-10 and for 10-5 (5%) electrode system
%
% Revision 1.1  2003/06/16 15:32:58  roberto
% new implementation, starting from code in preprocessing
%

if any(size(channel) == 0)
  % there is nothing to do if it is empty
  return
end

if isnumeric(channel)
  % change index into channelname
  channel = datachannel(channel);
  return
end

if ~iscell(channel)
  % ensure that a single input argument like 'all' also works
  channel = {channel};
end

% ensure that both inputs are column vectors
channel     = channel(:);
datachannel = datachannel(:);

% remove channels that occur more than once, this sorts the channels alphabetically
[channel, indx] = unique(channel);
% undo the sorting, make the order identical to that of the data channels
[dum, indx] = sort(indx);
channel = channel(indx);

% define the known groups with channel labels
labelall  = datachannel;
label1020 = {'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8',  'T7', 'C3', 'Cz', 'C4', 'T8',  'P7', 'P3', 'Pz', 'P4', 'P8', 'O1', 'Oz', 'O2'}';
label1010 = {'Fp1', 'Fpz', 'Fp2', 'AF9', 'AF7', 'AF5', 'AF3', 'AF1', 'AFz', 'AF2', 'AF4', 'AF6', 'AF8', 'AF10', 'F9', 'F7', 'F5', 'F3', 'F1', 'Fz', 'F2', 'F4', 'F6', 'F8', 'F10', 'FT9', 'FT7', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'FT8', 'FT10', 'T9', 'T7', 'C5', 'C3', 'C1', 'Cz', 'C2', 'C4', 'C6', 'T8', 'T10', 'TP9', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'CP6', 'TP8', 'TP10', 'P9', 'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10', 'PO9', 'PO7', 'PO5', 'PO3', 'PO1', 'POz', 'PO2', 'PO4', 'PO6', 'PO8', 'PO10', 'O1', 'Oz', 'O2', 'I1', 'Iz', 'I2'}';
label1005 = {'Fp1', 'Fpz', 'Fp2', 'AF9', 'AF7', 'AF5', 'AF3', 'AF1', 'AFz', 'AF2', 'AF4', 'AF6', 'AF8', 'AF10', 'F9', 'F7', 'F5', 'F3', 'F1', 'Fz', 'F2', 'F4', 'F6', 'F8', 'F10', 'FT9', 'FT7', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'FT8', 'FT10', 'T9', 'T7', 'C5', 'C3', 'C1', 'Cz', 'C2', 'C4', 'C6', 'T8', 'T10', 'TP9', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'CP6', 'TP8', 'TP10', 'P9', 'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10', 'PO9', 'PO7', 'PO5', 'PO3', 'PO1', 'POz', 'PO2', 'PO4', 'PO6', 'PO8', 'PO10', 'O1', 'Oz', 'O2', 'I1', 'Iz', 'I2', 'AFp9h', 'AFp7h', 'AFp5h', 'AFp3h', 'AFp1h', 'AFp2h', 'AFp4h', 'AFp6h', 'AFp8h', 'AFp10h', 'AFF9h', 'AFF7h', 'AFF5h', 'AFF3h', 'AFF1h', 'AFF2h', 'AFF4h', 'AFF6h', 'AFF8h', 'AFF10h', 'FFT9h', 'FFT7h', 'FFC5h', 'FFC3h', 'FFC1h', 'FFC2h', 'FFC4h', 'FFC6h', 'FFT8h', 'FFT10h', 'FTT9h', 'FTT7h', 'FCC5h', 'FCC3h', 'FCC1h', 'FCC2h', 'FCC4h', 'FCC6h', 'FTT8h', 'FTT10h', 'TTP9h', 'TTP7h', 'CCP5h', 'CCP3h', 'CCP1h', 'CCP2h', 'CCP4h', 'CCP6h', 'TTP8h', 'TTP10h', 'TPP9h', 'TPP7h', 'CPP5h', 'CPP3h', 'CPP1h', 'CPP2h', 'CPP4h', 'CPP6h', 'TPP8h', 'TPP10h', 'PPO9h', 'PPO7h', 'PPO5h', 'PPO3h', 'PPO1h', 'PPO2h', 'PPO4h', 'PPO6h', 'PPO8h', 'PPO10h', 'POO9h', 'POO7h', 'POO5h', 'POO3h', 'POO1h', 'POO2h', 'POO4h', 'POO6h', 'POO8h', 'POO10h', 'OI1h', 'OI2h', 'Fp1h', 'Fp2h', 'AF9h', 'AF7h', 'AF5h', 'AF3h', 'AF1h', 'AF2h', 'AF4h', 'AF6h', 'AF8h', 'AF10h', 'F9h', 'F7h', 'F5h', 'F3h', 'F1h', 'F2h', 'F4h', 'F6h', 'F8h', 'F10h', 'FT9h', 'FT7h', 'FC5h', 'FC3h', 'FC1h', 'FC2h', 'FC4h', 'FC6h', 'FT8h', 'FT10h', 'T9h', 'T7h', 'C5h', 'C3h', 'C1h', 'C2h', 'C4h', 'C6h', 'T8h', 'T10h', 'TP9h', 'TP7h', 'CP5h', 'CP3h', 'CP1h', 'CP2h', 'CP4h', 'CP6h', 'TP8h', 'TP10h', 'P9h', 'P7h', 'P5h', 'P3h', 'P1h', 'P2h', 'P4h', 'P6h', 'P8h', 'P10h', 'PO9h', 'PO7h', 'PO5h', 'PO3h', 'PO1h', 'PO2h', 'PO4h', 'PO6h', 'PO8h', 'PO10h', 'O1h', 'O2h', 'I1h', 'I2h', 'AFp9', 'AFp7', 'AFp5', 'AFp3', 'AFp1', 'AFpz', 'AFp2', 'AFp4', 'AFp6', 'AFp8', 'AFp10', 'AFF9', 'AFF7', 'AFF5', 'AFF3', 'AFF1', 'AFFz', 'AFF2', 'AFF4', 'AFF6', 'AFF8', 'AFF10', 'FFT9', 'FFT7', 'FFC5', 'FFC3', 'FFC1', 'FFCz', 'FFC2', 'FFC4', 'FFC6', 'FFT8', 'FFT10', 'FTT9', 'FTT7', 'FCC5', 'FCC3', 'FCC1', 'FCCz', 'FCC2', 'FCC4', 'FCC6', 'FTT8', 'FTT10', 'TTP9', 'TTP7', 'CCP5', 'CCP3', 'CCP1', 'CCPz', 'CCP2', 'CCP4', 'CCP6', 'TTP8', 'TTP10', 'TPP9', 'TPP7', 'CPP5', 'CPP3', 'CPP1', 'CPPz', 'CPP2', 'CPP4', 'CPP6', 'TPP8', 'TPP10', 'PPO9', 'PPO7', 'PPO5', 'PPO3', 'PPO1', 'PPOz', 'PPO2', 'PPO4', 'PPO6', 'PPO8', 'PPO10', 'POO9', 'POO7', 'POO5', 'POO3', 'POO1', 'POOz', 'POO2', 'POO4', 'POO6', 'POO8', 'POO10', 'OI1', 'OIz', 'OI2'}';
labelchwilla = {'Fz', 'Cz', 'Pz', 'F7', 'F8', 'LAT', 'RAT', 'LT', 'RT', 'LTP', 'RTP', 'OL', 'OR', 'FzA', 'Oz', 'F7A', 'F8A', 'F3A', 'F4A', 'F3', 'F4', 'P3', 'P4', 'T5', 'T6', 'P3P', 'P4P'}';
labelbham = {'P9', 'PPO9h', 'PO7', 'PPO5h', 'PPO3h', 'PO5h', 'POO9h', 'PO9', 'I1', 'OI1h', 'O1', 'POO1', 'PO3h', 'PPO1h', 'PPO2h', 'POz', 'Oz', 'Iz', 'I2', 'OI2h', 'O2', 'POO2', 'PO4h', 'PPO4h', 'PO6h', 'POO10h', 'PO10', 'PO8', 'PPO6h', 'PPO10h', 'P10', 'P8', 'TPP9h', 'TP7', 'TTP7h', 'CP5', 'TPP7h', 'P7', 'P5', 'CPP5h', 'CCP5h', 'CP3', 'P3', 'CPP3h', 'CCP3h', 'CP1', 'P1', 'Pz', 'CPP1h', 'CPz', 'CPP2h', 'P2', 'CPP4h', 'CP2', 'CCP4h', 'CP4', 'P4', 'P6', 'CPP6h', 'CCP6h', 'CP6', 'TPP8h', 'TP8', 'TPP10h', 'T7', 'FTT7h', 'FT7', 'FC5', 'FCC5h', 'C5', 'C3', 'FCC3h', 'FC3', 'FC1', 'C1', 'CCP1h', 'Cz', 'FCC1h', 'FCz', 'FFC1h', 'Fz', 'FFC2h', 'FC2', 'FCC2h', 'CCP2h', 'C2', 'C4', 'FCC4h', 'FC4', 'FC6', 'FCC6h', 'C6', 'TTP8h', 'T8', 'FTT8h', 'FT8', 'FT9', 'FFT9h', 'F7', 'FFT7h', 'FFC5h', 'F5', 'AFF7h', 'AF7', 'AF5h', 'AFF5h', 'F3', 'FFC3h', 'F1', 'AF3h', 'Fp1', 'Fpz', 'Fp2', 'AFz', 'AF4h', 'F2', 'FFC4h', 'F4', 'AFF6h', 'AF6h', 'AF8', 'AFF8h', 'F6', 'FFC6h', 'FFT8h', 'F8', 'FFT10h', 'FT10'};
labelref  = {'M1', 'M2', 'LM', 'RM', 'A1', 'A2'}';
labeleog  = datachannel(strncmp('EOG', datachannel, length('EOG')));  % anything that starts with EOG
labeleog  = {labeleog{:} 'HEOG', 'VEOG', 'VEOG-L', 'VEOG-R'}';        % or any of these
labelemg  = datachannel(strncmp('EMG', datachannel, length('EMG')));
labeleeg  = datachannel(strncmp('EEG', datachannel, length('EEG')));
labelmeg  = datachannel(strncmp('M'  , datachannel, length('M'  )));	% all MEG channels start with "M"
labelmz   = datachannel(strncmp('MZ' , datachannel, length('MZ' )));	% central MEG channels
labelml   = datachannel(strncmp('ML' , datachannel, length('ML' )));	% left    MEG channels
labelmr   = datachannel(strncmp('MR' , datachannel, length('MR' )));	% right   MEG channels
labelmlc  = datachannel(strncmp('MLC', datachannel, length('MLC')));
labelmlf  = datachannel(strncmp('MLF', datachannel, length('MLF')));
labelmlo  = datachannel(strncmp('MLO', datachannel, length('MLO')));
labelmlp  = datachannel(strncmp('MLP', datachannel, length('MLP')));
labelmlt  = datachannel(strncmp('MLT', datachannel, length('MLT')));
labelmrc  = datachannel(strncmp('MRC', datachannel, length('MRC')));
labelmrf  = datachannel(strncmp('MRF', datachannel, length('MRF')));
labelmro  = datachannel(strncmp('MRO', datachannel, length('MRO')));
labelmrp  = datachannel(strncmp('MRP', datachannel, length('MRP')));
labelmrt  = datachannel(strncmp('MRT', datachannel, length('MRT')));
labelmzc  = datachannel(strncmp('MZC', datachannel, length('MZC')));
labelmzf  = datachannel(strncmp('MZF', datachannel, length('MZF')));
labelmzo  = datachannel(strncmp('MZO', datachannel, length('MZO')));
labelmzp  = datachannel(strncmp('MZP', datachannel, length('MZP')));
labellfp  = datachannel(strncmp('lfp', datachannel, length('lfp')));
labelmua  = datachannel(strncmp('mua', datachannel, length('mua')));
labelspike  = datachannel(strncmp('spike', datachannel, length('spike')));

% figure out if there are bad channels or channel groups that should be excluded
findbadchannel = strncmp('-', channel, length('-'));      % bad channels start with '-'
badchannel = channel(findbadchannel);
if ~isempty(badchannel)
  for i=1:length(badchannel)
    badchannel{i} = badchannel{i}(2:end);                 % remove the '-' from the channel label
  end
  badchannel = channelselection(badchannel, datachannel); % support exclusion of channel groups
  channel(findbadchannel) = [];                           % remove them from the channels to be processed
end

% determine if any of the known groups is mentioned in the channel list
findall        = find(strcmp(channel, 'all'));
findmeg        = find(strcmp(channel, 'MEG'));
findemg        = find(strcmp(channel, 'EMG'));
findeeg        = find(strcmp(channel, 'EEG'));
findeeg1020    = find(strcmp(channel, 'EEG1020'));
findeeg1010    = find(strcmp(channel, 'EEG1010'));
findeeg1005    = find(strcmp(channel, 'EEG1005'));
findeegchwilla = find(strcmp(channel, 'EEGCHWILLA'));
findeegbham    = find(strcmp(channel, 'EEGBHAM'));
findeegref     = find(strcmp(channel, 'EEGREF'));
findeog        = find(strcmp(channel, 'EOG'));
findmz         = find(strcmp(channel, 'MZ' ));
findml         = find(strcmp(channel, 'ML' ));
findmr         = find(strcmp(channel, 'MR' ));
findmlc        = find(strcmp(channel, 'MLC'));
findmlf        = find(strcmp(channel, 'MLF'));
findmlo        = find(strcmp(channel, 'MLO'));
findmlp        = find(strcmp(channel, 'MLP'));
findmlt        = find(strcmp(channel, 'MLT'));
findmrc        = find(strcmp(channel, 'MRC'));
findmrf        = find(strcmp(channel, 'MRF'));
findmro        = find(strcmp(channel, 'MRO'));
findmrp        = find(strcmp(channel, 'MRP'));
findmrt        = find(strcmp(channel, 'MRT'));
findmzc        = find(strcmp(channel, 'MZC'));
findmzf        = find(strcmp(channel, 'MZF'));
findmzo        = find(strcmp(channel, 'MZO'));
findmzp        = find(strcmp(channel, 'MZP'));
findlfp        = find(strcmp(channel, 'lfp'));
findmua        = find(strcmp(channel, 'mua'));
findspike      = find(strcmp(channel, 'spike'));
findgui        = find(strcmp(channel, 'gui'));

% remove any occurences of groups in the channel list
channel([
  findall
  findmeg
  findemg
  findeeg
  findeeg1020
  findeeg1010
  findeeg1005
  findeegchwilla
  findeegbham
  findeegref
  findeog
  findmz
  findml
  findmr
  findmlc
  findmlf
  findmlo
  findmlp
  findmlt
  findmrc
  findmrf
  findmro
  findmrp
  findmrt
  findmzc
  findmzf
  findmzo
  findmzp
  findlfp
  findmua
  findspike
  findgui
]) = [];

% add the full channel labels to the channel list
if findall,        channel = [channel; labelall]; end
if findmeg,        channel = [channel; labelmeg]; end
if findemg,        channel = [channel; labelemg]; end
if findeeg,        channel = [channel; labeleeg]; end
if findeeg1020,    channel = [channel; label1020]; end
if findeeg1010,    channel = [channel; label1010]; end
if findeeg1005,    channel = [channel; label1005]; end
if findeegchwilla, channel = [channel; labelchwilla]; end
if findeegbham,    channel = [channel; labelbham]; end
if findeegref,     channel = [channel; labelref]; end
if findeog,        channel = [channel; labeleog]; end
if findmz ,        channel = [channel; labelmz ]; end
if findml ,        channel = [channel; labelml ]; end
if findmr ,        channel = [channel; labelmr ]; end
if findmlc,        channel = [channel; labelmlc]; end
if findmlf,        channel = [channel; labelmlf]; end
if findmlo,        channel = [channel; labelmlo]; end
if findmlp,        channel = [channel; labelmlp]; end
if findmlt,        channel = [channel; labelmlt]; end
if findmrc,        channel = [channel; labelmrc]; end
if findmrf,        channel = [channel; labelmrf]; end
if findmro,        channel = [channel; labelmro]; end
if findmrp,        channel = [channel; labelmrp]; end
if findmrt,        channel = [channel; labelmrt]; end
if findmzc,        channel = [channel; labelmzc]; end
if findmzf,        channel = [channel; labelmzf]; end
if findmzo,        channel = [channel; labelmzo]; end
if findmzp,        channel = [channel; labelmzp]; end
if findlfp,        channel = [channel; labellfp]; end
if findmua,        channel = [channel; labelmua]; end
if findspike,      channel = [channel; labelspike]; end

% remove channel labels that have been excluded by the user
badindx = match_str(channel, badchannel);
channel(badindx) = [];

% remove channel labels that are not present in the data
chanindx = match_str(channel, datachannel);

channel  = channel(chanindx);

if findgui
    indx = select_channel_list(datachannel, match_str(datachannel, channel), 'Select channels');
    channel = datachannel(indx);
end

% remove channels that occur more than once, this sorts the channels alphabetically
[channel, indx] = unique(channel);

% undo the sorting, make the order identical to that of the data channels
[dum, indx] = sort(indx);
channel = channel(indx);
