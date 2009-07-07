% demonstration of the combined use of EEG and MEG data

% clear data and fix the random number generator (RNG)

fclose('all');
close all
clear all

format long

% fixing the RNG in order to reproduce the experiment
rand('twister',1); randn('state',1);

% specify transformations of spectral data

% Here, we may preselect only a few channels, times, and/or frequencies,
% and we may average over channels, time, or frequencies.

cfg.frequency   = [5 30]; %[1 150]; % a range in data{i}.freq
cfg.avgoverchan = 'no'; % default
cfg.avgoverfreq = 'no'; % default
cfg.avgovertime = 'no'; % default

%% specify classification procedure

% myproc = clfproc({ ...
%     preprocessor('prefun',@(x)(log10(x))) ...
%     standardizer('bstd',false) ...
%     combiner('procedure',clfproc({da()}),'combination','concatenate') ...    %combiner('procedure',clfproc({filterer('verbose',true) nb()}),'combination','concatenate') ... % combines classification of multiple datasets
%     lr() ... 
% });

% myproc = clfproc({ ...
%     preprocessor('prefun',@(x)(log10(x))) ...
%     standardizer() ...
%     combiner('procedure',clfproc({da()}),'combination','product') ...    %combiner('procedure',clfproc({filterer('verbose',true) nb()}),'combination','concatenate') ... % combines classification of multiple datasets    
% });

myproc = clfproc({ ...
    preprocessor('prefun',@(x)(log10(x))) ...
    standardizer() ...
    combiner('procedure',clfproc({filterer('maxfeatures',100,'verbose',true,'validator',crossvalidator('procedure',clfproc({nb()}),'cvfolds',0.8)) nb()}),'combination','concatenate') ... % combines classification of multiple datasets    
    lr() ...
});

% validation method
cv = crossvalidator('procedure',myproc,'cvfolds',10,'randomize',true,'verbose',true);

%% load data and perform FieldTrip transformations

datas = cell(1,2);
designs = cell(1,2);
for d=1:2

    % load experimental data for a subject
    if d==1
        data = load('~/data/christian/eegcharlotte/freqftrs_eeg_hjorth_test2.mat','dat'); data = data.dat;
  %      cfg.channel     = {'C3' 'C4' 'Cz'};

    else
        data = load('~/data/christian/megcharlotte/freqftrs_meg_sensor_test2.mat','dat'); data = data.dat;
 %       cfg.channel     = {'MLF42' 'MRF42'};

    end

    data = { data{2} data{4}}; % select conditions from {Rx, Ri, Lx, Li, No}

    % view times as repetitions

    for c=1:length(data)
        data{c}.time = 1;
        data{c}.dimord = 'rpt_chan_freq';
        sz = size(data{c}.powspctrm);
        biol = zeros(sz(1)*4,sz(2),sz(3));
        for j=1:4
            biol(((j-1)*sz(1)+1):(j*sz(1)),:,:) = squeeze(data{c}.powspctrm(:,:,:,j));
        end
        data{c}.powspctrm = biol;
    end

    % possible other transformations may be specified here
    % (baseline correction etc)

    [cfg2,data] = prepare_timefreq_data(cfg,data{:});

    % make data suitable for classification

    % create design matrix; this specifies at least the class labels (1,2,3) for each example
    design = [];
    for c=1:length(data.biol)
        design = [design; c*ones(size(data.biol{c},1),1)];
    end
   
    % concatenate the data into one dataset
    datas{d} = cat(1,data.biol{:});
    
end

%% compute crossvalidation results

cv = cv.validate(datas(:),design);

%% get average classification accuracy

evaluate(cv.post,cv.design,'metric','accuracy')














