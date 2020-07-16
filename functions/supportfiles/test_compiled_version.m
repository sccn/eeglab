% EEGLAB history file generated on the 08-Jul-2020
% ------------------------------------------------
EEG = pop_loadset( 'filename', 'eeglab_data.set', 'filepath', '/data/matlab/eeglab/sample_data/');
EEG = eeg_checkset( EEG );
EEG = pop_loadset('filename','eeglab_data.set','filepath','/data/matlab/eeglab/sample_data/');
EEG = eeg_checkset( EEG );
EEG.etc.eeglabvers = 'development head'; % this tracks which version of EEGLAB is being used, you may ignore it
EEG = pop_loadset('filename','eeglab_data.set','filepath','/data/matlab/eeglab/sample_data/');
EEG = eeg_checkset( EEG );
EEG=pop_chanedit(EEG, 'lookup','/data/matlab/eeglab/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp','load',{'/data/matlab/eeglab/sample_data/eeglab_chan32.locs' 'filetype' 'autodetect'});
EEG = eeg_checkset( EEG );
EEG.etc.eeglabvers = '2020.0'; % this tracks which version of EEGLAB is being used, you may ignore it
EEG = pop_loadset('filename','eeglab_data.set','filepath','/data/matlab/eeglab/sample_data/');
EEG = eeg_checkset( EEG );
EEG=pop_chanedit(EEG, 'lookup','/data/matlab/eeglab/plugins/dipfit/standard_BEM/elec/standard_1005.elc');
EEG = eeg_checkset( EEG );
EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on', 'pca', 10);
EEG = eeg_checkset( EEG );
EEG = pop_iclabel(EEG, 'default');
EEG = eeg_checkset( EEG );
EEG = pop_dipfit_settings( EEG, 'hdmfile','/data/matlab/eeglab/plugins/dipfit/standard_BEM/standard_vol.mat','coordformat','MNI','mrifile','/data/matlab/eeglab/plugins/dipfit/standard_BEM/standard_mri.mat','chanfile','/data/matlab/eeglab/plugins/dipfit/standard_BEM/elec/standard_1005.elc','coord_transform',[0 0 0 0 0 -1.5708 1 1 1] ,'chansel',[1:31] );
EEG = eeg_checkset( EEG );
EEG = pop_dipfit_gridsearch(EEG, [1:31] ,[-85     -77.6087     -70.2174     -62.8261     -55.4348     -48.0435     -40.6522     -33.2609     -25.8696     -18.4783      -11.087     -3.69565      3.69565       11.087      18.4783      25.8696      33.2609      40.6522      48.0435      55.4348      62.8261      70.2174      77.6087           85] ,[-85     -77.6087     -70.2174     -62.8261     -55.4348     -48.0435     -40.6522     -33.2609     -25.8696     -18.4783      -11.087     -3.69565      3.69565       11.087      18.4783      25.8696      33.2609      40.6522      48.0435      55.4348      62.8261      70.2174      77.6087           85] ,[0      7.72727      15.4545      23.1818      30.9091      38.6364      46.3636      54.0909      61.8182      69.5455      77.2727           85] ,0.4);
EEG = eeg_checkset( EEG );
EEG = pop_headplot(EEG, 0, [1:3] , 'Components of dataset: Continuous EEG Data', [2  2], 'setup',{'/Users/arno/eeglab_data.spl','meshfile','mheadnew.mat','transform',[-1.136 7.7523 11.4527 -0.027117 0.015531 -1.5455 0.91234 0.93161 0.80698] });
EEG = eeg_checkset( EEG );

% STUDY
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
pop_editoptions( 'option_storedisk', 1);
[STUDY ALLEEG] = pop_loadstudy('filename', 'stern3s.study', 'filepath', '/data/data/STUDIES/STERN');
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
[STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG, {},'savetrials','on','interp','on','recompute','on','erp','on');

STUDY = std_erpplot(STUDY,ALLEEG,'channels',{'LEYE'}, 'design', 3);
STUDY = pop_statparams(STUDY, 'mode','fieldtrip');
STUDY = std_erpplot(STUDY,ALLEEG,'channels',{'LEYE'}, 'design', 3);
STUDY = std_makedesign(STUDY, ALLEEG, 3, 'name','Design 3','delfiles','off','defaultdesign','off','subjselect',{'S01','S02','S03'},'filepath','/data/data/STUDIES/STERN');
STUDY = std_selectdesign(STUDY, ALLEEG, 2);
STUDY = std_erpplot(STUDY,ALLEEG,'channels',{'LEYE'}, 'design', 2);
eeglab redraw;
