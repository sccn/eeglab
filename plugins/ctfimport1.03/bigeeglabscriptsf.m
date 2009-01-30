clear all;

folder='sampledataset.ds';

[data]=readepochs(folder);
%wind = [-.1 .15];
%[data]=readepochs(folder,'markers',{'click'},'window',wind);
%[data]=readepochs(folder,'window',wind,'trials',[1:10],'megsens',[2:7]);
if ~exist('wind','var')
    wind = data.setup.start_sec;
end
dat=permute(data.epochs{1},[2,1,3]);
size(dat)
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

EEG = pop_importdata( 'nbchan', size(dat,1), 'dataformat', 'array', 'data', 'dat', ...
  'pnts', size(dat,2), 'srate', data.setup.sample_rate, 'xmin', wind(1));

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname','test');
eeglab redraw;

for i=1:size(dat,1)
     EEG.chanlocs(i).labels = char(data.sensor.label(i));
     EEG.chanlocs(i).X      = data.sensor.location(1,i);
     EEG.chanlocs(i).Y      = data.sensor.location(2,i);
     EEG.chanlocs(i).Z      = data.sensor.location(3,i);
end
EEG.chanlocs = pop_chanedit(EEG.chanlocs, 'convert', 'cart2topo', 'convert', 'cart2sph');

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

eeglab redraw;
