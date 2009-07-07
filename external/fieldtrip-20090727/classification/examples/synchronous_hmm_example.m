% This Matlab script demonstrates how to use the neuroML toolbox
% together with FieldTrip in order to classify EEG/MEG data based on
% single trials. 

% here we demonstrate hidden markov models for synchronous experiments.
% i.e., we have finite horizon temporal data. Note that this is typically
% the case for neuroimaging data. I.e., we have a fixed number of trials of
% finite horizon

% Synchronous temporal data has dimensions N x (M x T) with N examples, M features, and 
% T time slices. In a synchronous experiment, each example is a finite horizon dynamic
% model of which the underlying state is assumed to be constant. In
% general, the data should have time as the last dimension

% In an asynchronous experiment, data has dimensions T x M. In this case, 
% examples should not be randomized during validation!
% the rationale for this distinction is that in synchronous experiments one
% example extends over time (requiring a finite horizon model) whereas for
% asynchronous experiments, one example is identical to one point in time
% (i.e., an example does not extend over time, but subsquent examples are
% coupled).

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
cfg.frequency   = [8 12]; %[1 150]; % a range in data{i}.freq
cfg.avgoverchan = 'no'; % default
cfg.avgoverfreq = 'yes'; % default
cfg.avgovertime = 'no'; % default

%% specify classification procedure

% does this preprocessing lead to multiplication interactions due to log ?????

% for finite horizon models we may also define a baseline procedure

myproc = clfproc({ ...    
    preprocessor('prefun',@(x)(log10(x))) ...
    standardizer() ...
    hmm('horizon',5,'coupled',false,'ar',false,'mixture',1,'ie',@canonical_jtree_ie) ...
    });

% validation method
cv = crossvalidator('procedure',myproc,'cvfolds',0.9);

%% load data and perform FieldTrip transformations

% load experimental data for a subject
load ~/data/christian/freqbcicomp/freq9.mat
data = freqtraindata{1};

% % time as examples
% data.time = 1;
% data.dimord = 'rpt_chan_freq';
% sz = size(data.powspctrm);
% biol = zeros(sz(1)*sz(4),sz(2),sz(3));
% for j=1:sz(4)
%     biol(((j-1)*sz(1)+1):(j*sz(1)),:,:) = squeeze(data.powspctrm(:,:,:,j));
% end
% data.powspctrm = biol;

% the associated design
design = trainlabels{1};

% replicate for the labels
% design = zeros(size(adesign,1)*sz(4),1);
% for j=1:size(adesign,1)
%     design(((j-1)*sz(4)+1):j*sz(4)) = adesign(j);
% end

% possible other transformations may be specified here
% (baseline correction etc)
[cfg,data] = prepare_timefreq_data(cfg,data);
data = data.biol{:};

% select a number of time slices
data = data(:,:,:,10:14);

%% compute crossvalidation results

cv = cv.validate(data,design);

%% get average classification accuracy

cv.evaluate('metric','accuracy')














