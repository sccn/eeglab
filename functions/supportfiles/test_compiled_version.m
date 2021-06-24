function test_compiled_version(varargin)
global CURRENTSTUDY CURRENTSET ALLEEG EEG STUDY

res = questdlg2('Are you seing the command line outputs?', 'Compiled version', 'No', 'Yes', 'Yes');
if strcmpi(res, 'No')
    res = questdlg2( [ 'Users need to be able to see command line' 10 ...
        'output even for the compiled version. This' 10 ...
        'require fixing and recompiling. Continue anyway?' ], 'Compiled version', 'No', 'Yes', 'Yes');
    if strcmpi(res, 'No')
        return
    end
end

try
    % EEGLAB history file generated on the 08-Jul-2020
    % ------------------------------------------------
    EEG = pop_loadset( 'filename', 'eeglab_data.set', 'filepath', 'sample_data');
    EEG = eeg_checkset( EEG );
    EEG = pop_loadset('filename','eeglab_data.set','filepath','sample_data');
    EEG = eeg_checkset( EEG );
    EEG.etc.eeglabvers = 'development head'; % this tracks which version of EEGLAB is being used, you may ignore it
    EEG = pop_loadset('filename','eeglab_data.set','filepath','sample_data');
    EEG = eeg_checkset( EEG );
    EEG=pop_chanedit(EEG, 'lookup',fullfile('plugins','dipfit','standard_BESA','standard-10-5-cap385.elp'),'load',{fullfile('sample_data','eeglab_chan32.locs') 'filetype' 'autodetect'});
    EEG = eeg_checkset( EEG );
    EEG.etc.eeglabvers = '2020.0'; % this tracks which version of EEGLAB is being used, you may ignore it
    EEG = pop_loadset('filename','eeglab_data.set','filepath','sample_data');
    EEG = eeg_checkset( EEG );
    EEG=pop_chanedit(EEG, 'lookup',fullfile('plugins','dipfit','standard_BEM','elec', 'standard_1005.elc'));
    EEG = eeg_checkset( EEG );
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
    EEG = eeg_checkset( EEG );
    EEG = eeg_checkset( EEG );
    %EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on', 'pca', 10);
    EEG = pop_runica(EEG, 'icatype', 'picard', 'pca', 10);
    EEG = eeg_checkset( EEG );
    EEG = pop_iclabel(EEG, 'default');
    EEG = eeg_checkset( EEG );
    EEG = pop_dipfit_settings( EEG, 'hdmfile',fullfile('plugins','dipfit','standard_BEM','standard_vol.mat'),'coordformat','MNI','mrifile',fullfile('plugins','dipfit','standard_BEM','standard_mri.mat'),'chanfile',fullfile('plugins','dipfit','standard_BEM','elec', 'standard_1005.elc'),'coord_transform',[0 0 0 0 0 -1.5708 1 1 1] ,'chansel',[1:31] );
    EEG = eeg_checkset( EEG );
    EEG = pop_multifit(EEG, 1,'threshold',100,'plotopt',{'normlen','on'});
    EEG = eeg_checkset( EEG );
    EEG = pop_headplot(EEG, 0, [1:3] , 'Components of dataset: Continuous EEG Data', [2  2], 'setup',{'eeglab_data.spl','meshfile','mheadnew.mat','transform',[-1.136 7.7523 11.4527 -0.027117 0.015531 -1.5455 0.91234 0.93161 0.80698] });
    close;
    EEG = pop_loadset( 'filename', 'eeglab_data_epochs_ica.set', 'filepath', 'sample_data');
    EEG = pop_dipfit_settings( EEG, 'hdmfile',fullfile('plugins','dipfit','standard_BEM','standard_vol.mat'),'coordformat','MNI','mrifile',fullfile('plugins','dipfit','standard_BEM','standard_mri.mat'),'chanfile',fullfile('plugins','dipfit','standard_BEM','elec', 'standard_1005.elc'),'coord_transform',[0.83215 -15.6287 2.4114 0.081214 0.00093739 -1.5732 1.1742 1.0601 1.1485],'chansel',[1:31] );
    pop_dipfit_loreta(EEG, 6);
    close
    
    % STUDY
    % -----
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    pop_editoptions( 'option_storedisk', 1);
    
    commands = {};
    EEG = pop_loadset('filename','eeglab_data_epochs_ica.set','filepath','sample_data');
    for iS = 1:10
        fileName1 = sprintf('eeglab_data_epochs_ica_c1_%d.set', iS);
        EEG2 = EEG; EEG2.data = EEG2.data+rand(EEG.nbchan, EEG.pnts, EEG.trials)*20-10;
        pop_saveset(EEG2, 'filename',fileName1);
        
        fileName2 = sprintf('eeglab_data_epochs_ica_c2_%d.set', iS);
        EEG2 = EEG; EEG2.data = EEG2.data+rand(EEG.nbchan, EEG.pnts, EEG.trials)*20-10;
        EEG2.data(:, 129:end,:) = EEG2.data(:, 129:end,:)+rand(EEG.nbchan, 256, EEG.trials)*20-5;
        pop_saveset(EEG2, 'filename',fileName2);
        
        commands = { commands{:} ...
            { 'index' iS*2-1 'load' fileName1 'subject' sprintf('S%d', iS) 'condition' '1' } ...
            { 'index' iS*2   'load' fileName2 'subject' sprintf('S%d', iS) 'condition' '2' } };
    end
    
    [STUDY, ALLEEG] = std_editset(STUDY, ALLEEG, 'commands', commands);
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
    STUDY = std_makedesign(STUDY, ALLEEG);
    
    [STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG, {},'savetrials','on','interp','on','recompute','on','erp','on');
    STUDY = std_erpplot(STUDY,ALLEEG,'channels',{'Fz'}, 'design', 1);
    close;
    
    STUDY = pop_statparams(STUDY, 'mode','fieldtrip', 'condstats', 'on');
    STUDY = std_erpplot(STUDY,ALLEEG,'channels',{'Fz'}, 'design', 1);
    close;
    
    % cleaning up
    delete('eeglab_data_epochs_ica_*.*');
    delete('*.daterp');
    delete('eeglab_data.spl');
    
catch lasterr
    % show error on command line
    disp(lasterr.getReport())
    warndlg2( [ 'Script failed with error below (see also command line). Fix and try again.' 10 10 ...
        lasterr.message ], 'Compiled version');
    return
end

warndlg2( [ 'Script run successfully.' ], 'Compiled version');

end
