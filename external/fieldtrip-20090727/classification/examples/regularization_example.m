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
cfg.frequency   = [11 13]; %[1 150]; % a range in data{i}.freq
cfg.avgoverchan = 'no'; % default
cfg.avgoverfreq = 'no'; % default
cfg.avgovertime = 'no'; % default

myproc = clfproc({ ...
    preprocessor('prefun',@(x)(log10(x))) ...
    standardizer() ...
    gslr('p',2) ...
    });

cv = crossvalidator('procedure',myproc,'cvfolds',1,'randomize',true,'verbose',true);


% load data and perform FieldTrip transformations

% load experimental data for a subject
data = load('~/data/christian/eegcharlotte/freqftrs_eeg_hjorth_test2.mat','dat'); data = data.dat;

data = { data{2} data{4} data{5}}; % select conditions from {Rx, Ri, Lx, Li, No}

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
[cfg,data] = prepare_timefreq_data(cfg,data{:});

% make data suitable for classification

% create design matrix; this specifies at least the class labels (1,2,3) for each example
design = [];
for c=1:length(data.biol)
    design = [design; c*ones(size(data.biol{c},1),1)];
end
%design = [design (1:size(design,1))'];

% concatenate the data into one dataset
data = cat(1,data.biol{:});

% no grouping

cv = cv.validate(data,design);
evaluate(cv.post,cv.design,'metric','accuracy')

% create plot

path = cv.procedure.clfmethods{3}.diagnostics.path(2:end);
lambdas = cv.procedure.clfmethods{3}.diagnostics.lambdas(2:end);
lambdas = full(lambdas);
lambdas = lambdas(lambdas > 0);

paths = zeros(length(path),size(path{1},1)*(size(path{1},2)-1));
for j=1:length(path)
    p = path{j}(:,1:(end-1));
    paths(j,:) = p(:);
end

stairs(lambdas,paths);
set(gca, 'XDir', 'reverse')
xlim([min(lambdas) max(lambdas)]);
ylim(1.1*[-max(abs(paths(:))) max(abs(paths(:)))]);
xlabel('\lambda');
ylabel('\theta');





