% ctfcomp
% EEGLAB script for joint time-frequency analysis for MEG data
% Saves single-channel ERP and JTF 
% Saves all outputs in .mat file
% Thomas Ferree @ UCSF
% Revised 1/28/2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;




% toggle saving eeglab .set, 0 = no, 1 = yes
saveSet = 0;

% set machine: 1 = TiBoard, 2 = Glutamate/Acetylcholine/Brainmap
machine = 2 ;

% user-defined paths
if machine == 1
	pathdata = '/Users/tferree/Documents/Data/Attention/HighN/';
	pathersp = '/Users/tferree/Documents/Research/Attention/Temp/';
	pathfigs = '/Users/tferree/Documents/Research/Attention/Temp/';
else
	addpath /netopt/lib/matlab/eeglab/ -end;
	addpath /home/tferree/Attention/HighN/ctf2eeglab/ -end;
	pathdata = '/home/cdale/HighN/HM_epoch/processed/';
	pathersp = '/home/tferree/Attention/HighN/Temp/';
	pathfigs = '/home/tferree/Attention/HighN/Temp/';
end

% call Darren's scripts
ctffolder = [pathdata 'HM_DNL_20040126_HiN_CueLValNT-fEOG.ds'];
ctfchannels = 'meg';
ctftime = 'all';
%ctftrials = [1:5:100];
ctftrials = 'all';
ctf = ctf_read(ctffolder,ctfchannels,ctftime,ctftrials);
ctf2eeglab;

% define JTFA parameters
winsize = 128; 
padratio = 2; % frequency resolution
bootalpha = 0.05; % already two-tailed
bootnaccu = 200; % 200-20000

% scale data upward to order unity
EEG.data = EEG.data * 10^12;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop on channels
Nchan = EEG.nbchan;
for ichan = 1:Nchan
        
    % channel string for saving figures to files with sequential file names
    if ichan < 10
        chanstring = ['00' num2str(ichan)];
    end
    if ichan > 9 & chan < 100
        chanstring = ['0' num2str(ichan)];
    end
    if ichan > 99
        chanstring = num2str(ichan);
    end
    
    % make single-channel ERP image
    h1 = figure(1000+ichan); 
    avewidth = 1;
    decimate = 1;
    pop_erpimage(EEG,1, [ichan],[],['Channel ' num2str(ichan)],avewidth,decimate,{},[],'' ,'erp','cbar');
    saveas(h1,[pathfigs 'ERP_Ch' chanstring '.jpg'],'jpg');
    close(h1);
    
    % make single-channel JTF image
    h2 = figure(2000+ichan);
    [ersp,itc,powbase,times,freqs,erspboot,itcboot]=pop_timef( EEG, 1, ichan, [-1200 2000], 0 ,'type', 'phasecoher', 'title',['Channel ' num2str(ichan)],'padratio', padratio, 'plotphase', 'off','winsize',winsize,'alpha',bootalpha,'naccu',bootnaccu,'baseline',-1200);
    saveas(h2,[pathfigs 'JTF_Ch' chanstring '.jpg'],'jpg');
    close(h2);
    
    % store all timef results in 'all' array, referring to all channels
    allersp(:,:,ichan) = ersp;
    allitc(:,:,ichan) = itc;
    allpowbase(:,:,ichan) = powbase;
    alltimes(:,:,ichan) = times;
    allfreqs(:,:,ichan) = freqs;        
    allerspboot(:,:,ichan) = erspboot;
    allitcboot(:,:,ichan) = itcboot;

end % end loop over channels

% save results to file - avoiding overwrite with 'old' designation
%if exist([pathersp 'ERSP.mat'])
%	mv [pathersp 'ERSP.mat'] [pathersp 'ERSP_old.mat'];
%end
folderpath = [pathersp 'ERSP.mat'];
save(folderpath,'allersp','allitc','allpowbase','alltimes','allfreqs','allerspboot','allitcboot');

