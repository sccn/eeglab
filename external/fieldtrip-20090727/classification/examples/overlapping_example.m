% demonstrate groupwise regularization using overlapping groups

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

cfg.channel     = {'C4'};
cfg.frequency   = [8 20]; %[1 150]; % a range in data{i}.freq
cfg.avgoverchan = 'no'; % default
cfg.avgoverfreq = 'no'; % default
cfg.avgovertime = 'no'; % default


% specify groups

% for j=1:18 % channels
% 
%     groups{g} = [];
%     for f=1:46
%         for k=1:3 % classes
%             groups{g} = [groups{g}; [k (f-1)*18+j]; ];
%         end
%     end
% 
%     g = g + 1;
% end

% mf = cfg.frequency(2) - cfg.frequency(1);
% groups = cell(1,mf);
% g=1;
% for f=2:mf
%     groups{g} = [ 1 f-1; 1 f; 1 f+1 ];
%     g = g+1;
% end
% groups{g} = [1 mf+2]; % bias term

mf = cfg.frequency(2) - cfg.frequency(1) + 1;
groups = cell(1,mf);
g=1;
for f=2:mf
    groups{g} = [ 1 f-1; 1 f ];
    g = g+1;
end
groups{g} = [1 mf+1]; % bias term


% regularization
isregularized = true(1,g);
isregularized(g) = false; % bias is unregularized

% specify classification procedure

lambda = [80.066517871922301  68.714000369253313  59.800557055773922  37.019967470982671  18.015480234313330  12.130584814946991   6.449333824048325   4.394067695411302   4.284177177252106 3.45448267579721 2.888874626090915 0.753533591155653];

myproc = clfproc({ ...
    preprocessor('prefun',@(x)(log10(x))) ...
    standardizer() ...
    gslr('p',2,'overlap',true,'isregularized',isregularized,'groups',groups) ...
    });

% validation method
cv = crossvalidator('procedure',myproc,'cvfolds',0.9,'randomize',true,'verbose',true);

% load data and perform FieldTrip transformations

% load experimental data for a subject
data = load('~/data/christian/eegcharlotte/freqftrs_eeg_hjorth_test2.mat','dat'); data = data.dat;

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
[cfg,data] = prepare_timefreq_data(cfg,data{:});

% make data suitable for classification

% create design matrix; this specifies at least the class labels (1,2,3) for each example
design = [];
for c=1:length(data.biol)
    design = [design; c*ones(size(data.biol{c},1),1)];
end

% concatenate the data into one dataset
data = cat(1,data.biol{:});

% compute crossvalidation results

%[post,design] = cv.validate({data data},{design design});
cv = cv.validate(data,design);

% get average classification accuracy

cv.evaluate('metric','accuracy')

% show groups

% imagesc(abs(cv.procedure.clfmethods{3}.model(1,1:(end-1))))
% 
% figure

bar(full(abs(cv.procedure.clfmethods{3}.model(1,1:(end-1)))))








