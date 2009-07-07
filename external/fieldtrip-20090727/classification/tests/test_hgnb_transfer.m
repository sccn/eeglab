% This Matlab script demonstrates how to use the neuroML toolbox
% together with FieldTrip in order to classify EEG/MEG data based on
% single trials. The dataset is already preprocessed and freqanalysed and
% consists of a cell array representing left motor execution (Lx), left motor
% imagery (Li), right motor execution (Rx), right motor imagery (Ri), and no motor
% response (No) classes.

%% clear data and fix the random number generator (RNG)

fclose('all');
%close all
clear all

format long

% fixing the RNG in order to reproduce the experiment
rand('twister',1); randn('state',1);

%% specify transformations of spectral data

% Here, we may preselect only a few channels, times, and/or frequencies,
% and we may average over channels, time, or frequencies.

cfg.channel     = {'C3' 'C4' 'Cz'};
cfg.frequency   = [8 14]; %[1 150]; % a range in data{i}.freq
cfg.avgoverchan = 'no'; % default
cfg.avgoverfreq = 'no'; % default
cfg.avgovertime = 'no'; % default

%% specify classification procedure

myproc = clfproc({ ...    
    preprocessor('prefun',@(x)(log10(x))) ...
    standardizer('bstd',false) ... 
    hgnb_transfer() ...
    });


% validation method
cv = crossvalidator('procedure',myproc,'cvfolds',10,'randomize',true,'verbose',true);

%% load data and perform FieldTrip transformations

% load experimental data for all subjects
subjects = [2 1];

datas = cell(1,length(subjects));
designs = cell(1,length(subjects));

for s=1:length(subjects)

    data = load(sprintf('~/data/bci/eegcharlotte/freqftrs_eeg_hjorth_test%d.mat',subjects(s)),'dat'); data = data.dat;

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
        data{c}.powspctrm = biol(:,:,:,:);
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
    datas{s} = cat(1,data.biol{:});
    designs{s} = design;
end
 
rand('twister',1); randn('state',1);

prm1 = randperm(length(designs{1}));
prm2 = randperm(length(designs{2}));

N = 100;

acctrans  = zeros(2,N);
acc       = zeros(2,N);
for j=20:10:N
 
  disp(j);
  
  data    = {datas{1}(prm1(:),:) datas{2}(prm2(1:j),:)};
  design  = {designs{1}(prm1(:),:) designs{2}(prm2(1:j),:)};
  
  % compute crossvalidation results
  rand('twister',1); randn('state',1);
  cv = cv.validate(data,design);

  acctrans(:,j) = transpose(cv.evaluate);

  % validation method
  cv1 = crossvalidator('procedure',clfproc({preprocessor('prefun',@(x)(log10(x))) standardizer('bstd',false) nb() }),'cvfolds',10,'randomize',true,'verbose',true);
  rand('twister',1); randn('state',1);
  cv1 = cv1.validate(data{1},design{1});
  acc(1,j) = cv1.evaluate;

  cv2 = crossvalidator('procedure',clfproc({preprocessor('prefun',@(x)(log10(x))) standardizer('bstd',false) nb() }),'cvfolds',10,'randomize',true,'verbose',true);
  rand('twister',1); randn('state',1);
  cv2 = cv2.validate(data{2},design{2});
  acc(2,j) = cv2.evaluate;  

end

% check if means are closer to each other!

m   = cv.getmodel(1);
m1  = cv1.getmodel(1);
m2  = cv2.getmodel(1);

bar([m1; m2; m{1}; m{2}]')
legend('model 1','model 2','transfer model 1','transfer model 2');

idx = find(acc(1,:)~=0);
plot(idx,[acc(:,idx); acctrans(:,idx)]')
legend('model 1','model 2','transfer model 1','transfer model 2');
