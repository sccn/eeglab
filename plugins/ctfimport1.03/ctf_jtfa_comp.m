% ctfcomp
% EEGLAB script for joint time-frequency analysis for MEG data
% Reads one channel at a time to minimize memory requirements
% Saves all single-channel ERSP figures including ITC
% Saves all numerical results in .mat file
% Thomas Ferree @ UCSF
% Revised 2/16/2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

% set machine: 1 = TiBoard, 2 = Glutamate/Acetylcholine/Brainmap, Corby
machine = 3;

% set CTF parameters
ctftrials = 'all';
%ctftrials = [1:10];
ctftime = 'all';

% set JTFA parameters
winsize = 128; % time window width (->256, Scott 4/21/2004)
padratio = 4; % frequency resolution (->2?, Scott 4/21/2004)
bootalpha = 0.05; % two-tailed p-value
bootnaccu = 2000; % 200-20000

% user-defined paths
if machine == 1
	pathdata = '/Users/tferree/Documents/Data/Attention/HighN/';
	pathtemp = '/Users/tferree/Documents/Research/Attention/HighN/Temp/';
end
if machine >= 2
	addpath /home/tferree/EEGLAB42/;
	addpath /home/tferree/Attention/HighN/ctf2eeglab/;
	pathdata = '/home/cdale/HighN/';
end
if machine == 2
	pathtemp = '/home/tferree/Attention/HighN/Temp/';
elseif machine == 3	
	pathtemp = '/home/cdale/Temp/';
end

% toggle saving eeglab .set, 0 = no, 1 = yes
saveSet = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%subject = 'HM'
%date = '20040126';

%subject = 'GN'
%date = '20040109';

subject = 'CJ'
date = '20040116';

% loop on cues
for icue = 1:2

    % define cue condition
    if icue == 1
        cue = 'CueL'
    else icue == 2
        cue = 'CueR'
    end
    
    % loop on channels
	for ichan = 1:275

        % initialize data vector
        datavector = [];
		ntrials = 0;
        
		% loop over four targets with correct response + one for concatination
		for itarget = 1:4
		
			% define target condition
			if itarget == 1
                target = 'InvNT'
			elseif itarget == 2 
                target = 'InvTarg'
			elseif itarget == 3
                target = 'ValNT'
			elseif itarget == 4 
                target = 'ValTarg'
			else
                error('Invalid index for target condition.')
			end
			
			% create output folders
			
			if ~exist([pathtemp 'ERSP_' subject(1:2) '_' cue])
                 mkdir(pathtemp,['ERSP_' subject(1:2) '_' cue]);
			end
			
% 			if ~exist([pathtemp subject(1:2) '_' cue filesep 'ERP'])
%                  mkdir(pathtemp, [subject(1:2) '_' cue filesep 'ERP']);
% 			end
% 			
% 			if ~exist([pathtemp subject(1:2) '_' cue filesep 'JTF'])
%                  mkdir(pathtemp, [subject(1:2) '_' cue filesep 'JTF']);
% 			end
			
			% load CTF data
			ctffolder = [pathdata subject '_epoch' filesep 'processed' filesep subject '_DNL_' date '_HiN_' cue target '-fEOG.ds'];
            ctf = ctf_read_res4(ctffolder);
            ctfchannels = ctf.sensor.index.meg;

			ctf = ctf_read(ctffolder,ctfchannels(ichan),ctftime,ctftrials);
			ctf2eeglab;
            
            % concatinate trials over target conditions
            ntrials = ntrials + EEG.trials;
            for trial = 1:EEG.trials
                datavector = [datavector EEG.data(1,:,trial)];
            end
            
        end % end loop over targets
    
        % scale data upward to order unity
        datavector = datavector * 10^12;
            
        % channel string for saving figures to files with sequential file names
        if ichan < 10
            chanstring = ['00' num2str(ichan)];
        end
        if ichan > 9 & ichan < 100
            chanstring = ['0' num2str(ichan)];
        end
        if ichan > 99
            chanstring = num2str(ichan);
        end
        
%         % make single-channel ERP image
%         h1 = figure(1000+ichan); 
%         avewidth = 1;
%         decimate = 1;
%         pop_erpimage(EEG,1, [ichan],[],['Channel ' num2str(ichan)],avewidth,decimate,{},[],'' ,'erp','cbar','yerplabel','pT');
%         saveas(h1,[pathtemp subject '_' cue filesep 'Ch' chanstring '.jpg'],'jpg');
%         close(h1);
        
        % make single-channel JTF image
        h2 = figure(2000+ichan);
        %[ersp,itc,powbase,times,freqs,erspboot,itcboot] = pop_timef( EEG, 1, ichan, [-1200 2000], 0 ,'type', 'phasecoher','baseline',-300, 'title',['Channel ' num2str(ichan)],'padratio', padratio, 'plotphase', 'off','winsize',winsize,'alpha',bootalpha,'naccu',bootnaccu);
        [ersp,itc,powbase,times,freqs,erspboot,itcboot] = timef(datavector, EEG.pnts, [-1200 2000], EEG.srate, 0,'type', 'phasecoher','baseline',-300, 'title',['Channel ' num2str(ichan)],'padratio', padratio, 'plotphase', 'off','winsize',winsize,'alpha',bootalpha,'naccu',bootnaccu);
        saveas(h2,[pathtemp 'ERSP_' subject '_' cue filesep 'Ch' chanstring '.jpg'],'jpg');
        close(h2);
        
        % store all timef results in 'all' arrays, referring to all channels
        allersp(:,:,ichan) = ersp;
        allitc(:,:,ichan) = itc;
        allpowbase(:,:,ichan) = powbase;
        alltimes(:,:,ichan) = times;
        allfreqs(:,:,ichan) = freqs;        
        allerspboot(:,:,ichan) = erspboot;
        allitcboot(:,:,ichan) = itcboot;
	
	end % end loop over channels
	
	folderpath = [pathtemp 'ERSP_' subject '_' cue '.mat'];
	save(folderpath,'allersp','allitc','allpowbase','alltimes','allfreqs','allerspboot','allitcboot','ntrials');
	
end % end loop over cues
