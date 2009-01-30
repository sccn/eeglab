% ctftopo
% EEGLAB script for plotting topo maps of CTF data
% Thomas Ferree @ UCSF
% Revised 2/17/2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

% toggle saving eeglab .set, 0 = no, 1 = yes
saveSet = 0;

% set machine: 1 = TiBoard, 2 = Glutamate/Acetylcholine/Brainmap
machine = 2;

% set input and output paths
if machine == 1
	pathdata = '/Users/tferree/Documents/Data/Attention/HighN/';
	pathersp = '/Users/tferree/Documents/Research/Attention/HighN/Numerical/';
	pathtemp = '/Users/tferree/Documents/Research/Attention/HighN/Temp/';
else
	addpath /home/tferree/EEGLAB42/;
	addpath /home/tferree/Attention/HighN/ctf2eeglab/;
	pathdata = '/home/cdale/HighN/';
	pathersp = '/home/tferree/Attention/HighN/Numerical/';
	pathtemp = '/home/tferree/Attention/HighN/Temp/';
end

% set subject and conditions

%subject = 'HM'
%date = '20040126';

%subject = 'GN'
%date = '20040109';

subject = 'CJ'
date = '20040116';

cue = 'CueL'

% load ERSP results
load([pathersp 'ERSP_' subject '_' cue '.mat']);

% call Darren's scripts
%ctffolder = [pathdata subject '_DNL_' date '_HiN_' cue 'InvNT' '-fEOG.ds'];
ctffolder = [pathdata subject '_epoch' filesep 'processed' filesep subject '_DNL_' date '_HiN_' cue 'InvNT' '-fEOG.ds'];

ctfchannels = 'meg';
ctftime = 'all';
ctftrials = [1];
ctf = ctf_read(ctffolder,ctfchannels,ctftime,ctftrials);
ctf2eeglab;

% make topo maps
freqs = allfreqs(1,:,1);
times = alltimes(1,:,1);    
maxcolor = 5;

for plottime = 0:100:1800
    
    if plottime <= 9
        timestring = ['000' num2str(plottime)];
    elseif plottime >= 10 & plottime <= 99
        timestring = ['00' num2str(plottime)];
    elseif plottime >= 100 & plottime < 1000
        timestring = ['0' num2str(plottime)];
    else
        timestring = num2str(plottime);
    end
    
    timefreqlist = [plottime 5; plottime 7; plottime 9; plottime 11; plottime 13; plottime 15; plottime 17];
    h3 = figure;
    tftopo(allersp,times,freqs,'mode','ave','signifs', allerspboot,'limits',[min(times) max(times) 1 20 -maxcolor, maxcolor],'vert',[0 1000],'chanlocs',EEG.chanlocs,'timefreqs',timefreqlist,'electrodes','off','maplimits',[-maxcolor,maxcolor]);
    saveas(h3,[pathtemp 'Topo_' subject '_' cue '_' timestring 'msec.jpg'],'jpg');
    close(h3);

end
