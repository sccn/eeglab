% This Matlab script demonstrates how to use the neuroML toolbox
% together with FieldTrip in order to classify EEG/MEG data based on
% single trials. The dataset is already preprocessed and freqanalysed and
% consists of a cell array representing left motor execution (Lx), left motor
% imagery (Li), right motor execution (Rx), right motor imagery (Ri), and no motor
% response (No) classes.

%% clear data and fix the random number generator (RNG)

fclose('all');
close all
clear all

format long

% fixing the RNG in order to reproduce the experiment
rand('twister',1); randn('state',1);

%% specify transformations of spectral data

% Here, we may preselect only a few channels, times, and/or frequencies,
% and we may average over channels, time, or frequencies.

cfg.channel     = {'C3' 'C4'};
cfg.frequency   = [8 14]; %[1 150]; % a range in data{i}.freq
cfg.latency     = [0 4]; % should not be touched for this experiment
cfg.avgoverchan = 'no'; % default
%cfg.avgoverfreq = 'yes'; % default
%cfg.avgovertime = 'yes'; % default

%% specify classification procedure

% myproc = clfproc({ ...
%     preprocessor('prefun',@(x)(log10(x))) ...
%     standardizer() ...
%     one_against_one('procedure',clfproc({da()}))
%    });

myproc = clfproc({ ...
    preprocessor('prefun',@(x)(log10(x))) ...
    standardizer() ...
    one_against_one('procedure',clfproc({kernelmethod()}),'combination','majority')
   });

% validation method
cv = crossvalidator('procedure',myproc,'cvfolds',10,'randomize',true);

%% load data and perform FieldTrip transformations

% load experimental data for a subject
data = load('~/data/christian/eegcharlotte/freqftrs_eeg_hjorth_test2.mat','dat'); data = data.dat;
data = { data{2} data{4} data{5}}; % select conditions from {Rx, Ri, Lx, Li, No}

% possible other transformations may be specified here 
% (baseline correction etc)
[cfg,data] = prepare_timefreq_data(cfg,data{:});

%% make data suitable for classification

% create design matrix; this specifies at least the class labels (1,2,3) for each example
design = [];
for c=1:length(data.biol)
    design = [design; c*ones(size(data.biol{c},1),1)];
end

% concatenate the data into one dataset
data = cat(1,data.biol{:});

%% compute crossvalidation results

cv = cv.validate(data,design);

%% get average classification accuracy

evaluate(cv.post,cv.design,'metric','accuracy')














