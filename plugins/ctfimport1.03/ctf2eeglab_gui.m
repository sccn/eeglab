%GUI script for reading ctf data into eeglab.
%to run type ctf2eeglab on the command line.

clear all;
folder = uigetdir('*.ds','Data Set Finder');
[ctf2] = ctf_read_res4(folder,1);
sensors = {'meg','eeg','ref','vc'};
sensloc = 'ctf2.sensor.index';
sensorlen = '1:';
num2str(length(ctf2.sensor.index.meg_sens));

k = menu('Choose the channels you would like to use:','MEG','EEG','Reference','Virtual');
sensloc = strcat(sensloc,sensors(k));
sensorlen = strcat(sensorlen, num2str(length(eval(sensloc{1}))));  
prompt   = {'Enter Sensor Numbers to Read (ie 2:7 or 1 3 7 23 65)'};
title    = 'Input for Sensor Range (Default is all sensors)';
lines = 1;
def     = {sensorlen};
sensnum   = inputdlg(prompt,title,lines,def);
for i = str2num(sensnum{1});
    if ~ismember(str2num(def{1}),i)
        errordlg('There is no such sensor','Sensor Range Error')
        exit(1);
    end
end

trialnums = '1:';
trialnums = strcat(trialnums, num2str(ctf2.setup.number_trials));
prompt2 = {'Enter the trials you would like to use(ie 1:3 or 2 4 5):'};
def2 = {trialnums};
trials = inputdlg(prompt2, 'Input for trial range', 1,def2);
for i = str2num(trials{1})
    if ~ismember(str2num(def2{1}),i)
        errordlg('Trial not in dataset','Trial Error')
        exit(1);
    end
end
button = questdlg('Would you like to use markers?','Use Markers?','Yes','No','Yes');
if strcmp(button,'Yes')
    [marker_info] = readmarkerfile(folder);
    markers = marker_info.marker_names;
    m = menu('Please choose which marker you would like to use', markers);
end
prompt   = {'Enter the start time (sec):', 'Enter the end time:'};
title    = 'Time Window';
lines = 1;
def     = {'0', '.15'};
def(1) = {num2str(ctf2.setup.start_sec)};
def(2) = {num2str(ctf2.setup.end_sec)};
markeranswer = inputdlg(prompt,title,lines,def);
wind = [0 0];
wind(1) = str2num(markeranswer{1});
wind(2) = str2num(markeranswer{2});

if strcmp(button,'Yes')
    [data]=readepochs(folder,'marker_info',marker_info,'ctf',ctf2,sensors{k}, str2num(sensnum{1}),'trials',str2num(trials{1}), 'markers',markers(m),'window', wind);
else
    [data]=readepochs(folder,'ctf',ctf2,sensors{k}, str2num(sensnum{1}),'trials',str2num(trials{1}), 'window', wind);
end


%wind = [-.1 .15];
%[data]=readepochs(folder,'markers',{'click'},'window',wind);
%[data]=readepochs(folder,'window',wind,'trials',[1:10],'megsens',[2:7]);

dat=permute(data.epochs{1},[2,1,3]);
size(dat)
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

EEG = pop_importdata( 'nbchan', size(dat,1), 'dataformat', 'array', 'data', 'dat', 'pnts', size(dat,2), 'srate', data.setup.sample_rate, 'xmin', wind(1));
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,  'setname',data.setup.subject);
eeglab redraw;

for i=1:size(dat,1)
    if ~isempty(data.sensor.location)
        EEG.chanlocs(i).labels=char(data.sensor.label(i));
        EEG.chanlocs(i).X=data.sensor.location(1,i);
        EEG.chanlocs(i).Y=data.sensor.location(2,i);
        EEG.chanlocs(i).Z=data.sensor.location(3,i);
    end
end
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'convert', 'cart2topo', 'convert', 'cart2sph');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw;
clear dat data ctf2 marker_info markers sensors;
