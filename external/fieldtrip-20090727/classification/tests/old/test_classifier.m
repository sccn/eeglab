%% this script can be used to test your classifier

fclose('all');
close all
clear all
format long

% fix the RNG in order to reproduce the experiment
rand('twister',1); randn('state',1);

%% 
% load data
load ~/data/fieldtrip/covattfrq1
data = {left right}; 

%% 
% specify transformations of spectral data.
% Here, we may preselect only a few channels, times, and/or frequencies,
% and we may average over channels, time, or frequencies.
% we call a FieldTrip private function (this should be replaced at some
% point)
cfg.channel     = {'MLO34' 'MRO34'};
cfg.frequency   = [8 14];
cfg.latency     = [0.5 2.5];
cfg.avgoverfreq = 'no';
cfg.avgovertime = 'no';
[cfg,data] = prepare_timefreq_data(cfg,data{:});

nchan = length(data.label);

%%
% create design matrix
design = [];
for c=1:length(data.biol)
    design = [design; c*ones(size(data.biol{c},1),1)];
end

% concatenate the data into one dataset
data = cat(1,data.biol{:});

%% 
% specify classification procedure: insert your method here!

myproc = { ...    
    standardizer() ...
    svmmethod() ...
    };

%%
% validation method; randomize trial order and give verbose output

cv = crossvalidator('balanced',false,'procedure',myproc,'cvfolds',0.8,'randomize',true,'verbose',true);

%% 
% compute crossvalidation results
cv = cv.validate(data,design);

%% 
% get average classification accuracy

cv.evaluate('metric','accuracy')

%% 
% significance compared with a baseline classifier
cv.significance

