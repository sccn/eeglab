%% Using statistics_crossvalidate together with FieldTrip data
% This example demonstrates how to use neuroimaging data obtained from FieldTrip
% together with the statistics_crossvalidate function. In the example, we make use of covert
% attention data of one subject that has already been frequency analyzed.
%
% The data consists of 7 different frequencies at 274 channels at time
% points [-0.5 0 0.5 1 1.5 2 2.5]. We can expect evoked response after the
% cue and alpha modulation after about 1 second.
%
% Copyright (C) 2008  Marcel van Gerven
%

%% Classification examples using FieldTrip's statistics_crossvalidate
function statistics_crossvalidate_demo()

%% 
% Some initialization

fclose('all');
close all
clear all
format long

% fix the RNG in order to reproduce the experiment
rand('twister',1); randn('state',1);

%% 
% load data
load ~/code/classification/toolboxes/gerven/bayesbrain/examples/freqli;
load ~/code/classification/toolboxes/gerven/bayesbrain/examples/freqri;
data = {freqLI freqRI}; clear freqLI; clear freqRI;

%%
% initialize parameters

cfg = [];
cfg.channel     = {'MLO' 'MRO'};
cfg.frequency   = [8 30];
cfg.latency     = [0.5 2.5];
cfg.avgovertime = 'yes';

%%
% create design matrix; this specifies at least the class labels
% for each example
cfg.design = [];
for c=1:length(data)
    cfg.design = [cfg.design; c*ones(size(data{c}.powspctrm,1),1)];
end

%% 
% specify crossvalidate as statistical procedure and use discriminant
% analysis as a classifier
cfg.method = 'crossvalidate';
cfg.clfproc = clfproc({standardizer() svmmethod()});

%% 
% call crossvalidation through freqstatistics
stat = freqstatistics(cfg,data{:});

%%
% output classification rate

stat.prob

%% Classification using a different classifier

% logistic regression instead of discriminant analysis

cfg.clfproc = clfproc({standardizer() lr()});
stat = freqstatistics(cfg,data{:});

stat.prob

%% Using a different number of folds and a different metric

cfg.metric = 'auc'; % area under the ROC curve
cfg.cvfolds = 0.8; % use 80% for training and 20% for testing
stat = freqstatistics(cfg,data{:});
stat.prob

%%
% for this classifier we can project parameters back onto the scalp

cfg.layout = 'CTF275.lay';
cfg.zparam = 'model1';
topoplotTFR(cfg,stat);

end