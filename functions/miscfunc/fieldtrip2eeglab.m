% load data file ('dataf') preprocessed with fieldtrip
% and show in eeglab viewer
%
% This function is provided as is. It only works for some specific type of
% data. This is a simple function to help the developer and by no mean
% an all purpose function.

function [EEG] = fieldtrip2eeglab(dataf)

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

if exist(dataf,'file')
  load(dataf)
end

% load chanlocs.mat
% EEG.chanlocs = chanlocs;
EEG.chanlocs = [];

for i=1:size(data.trial,2)
  EEG.data(:,:,i) = single(data.trial{i});
end

EEG.setname    = dataf; %data.cfg.dataset;
EEG.filename   = '';
EEG.filepath   = '';
EEG.subject    = '';
EEG.group      = '';
EEG.condition  = '';
EEG.session    = [];
EEG.comments   = 'preprocessed with fieldtrip';
EEG.nbchan     = size(data.trial{1},1);
EEG.trials     = size(data.trial,2);
EEG.pnts       = size(data.trial{1},2);
EEG.srate      = data.fsample;
EEG.xmin       = data.time{1}(1);
EEG.xmax       = data.time{1}(end);
EEG.times      = data.time{1};
EEG.ref        = []; %'common';
EEG.event      = [];
EEG.epoch      = [];
EEG.icawinv    = [];
EEG.icasphere  = [];
EEG.icaweights = [];
EEG.icaact     = [];
EEG.saved      = 'no';

[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
eeglab redraw
pop_eegplot( EEG, 1, 1, 1);

 