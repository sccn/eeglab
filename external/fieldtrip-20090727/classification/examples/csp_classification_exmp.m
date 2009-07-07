% csp_classification_exmp.m
% 
% This script depicts a dichotomous classification process with the use of CSP features
% (CSPPROCESSOR object is used).
%

% Pawel Herman, 2009

%addpath(genpath('h:\common\matlab\fieldtrip\classification'));

datapath = 'D:\new4Pawel\JM\';
load(cat(2,datapath,'preprocmissesSimple.mat'));
load(cat(2,datapath,'preprochitsSimple.mat'));

HITdata  = preprochits;      clear preprochits;
MISSdata = preprocmisses;    clear preprocmisses;

disp(' ');
disp(sprintf('%3d hits',length(HITdata.trial)));
disp(sprintf('%3d misses',length(MISSdata.trial)));

cfg = [];
cfg.channel = {'MRO','MLO','MZO','MLT','MRT','MLP','MRP','MZP'};%,'MLC','MRC','MZC'};
cfg.bpfilter = 'yes';
cfg.bpfiltdir = 'twopass'; 
cfg.bpfilttype = 'but'; 
cfg.detrend = 'yes'; 
cfg.datatype = 'continuous'; 
cfg.keeptrials = 'yes'; 
cfg.bpfreq = [8 12];

% filtering in the predefined alpha band
filtHIT = preprocessing(cfg,HITdata);
filtMISS = preprocessing(cfg,MISSdata);

clear HITdata;
clear MISSdata;

assert( size(filtHIT.trial{1},1) == size(filtMISS.trial{1},1) );
nchan = size(filtHIT.trial{1},1);

disp(sprintf('Processing %d channels',nchan));

% extracting timelockanalysis structures
cfg = [];
cfg.keeptrials = 'yes';
tlHIT = timelockanalysis(cfg,filtHIT);
tlMISS = timelockanalysis(cfg,filtMISS);

clear filtHIT;
clear filtMISS;

% data (pre-filtered trials) and design
orig_data = [tlHIT.trial; tlMISS.trial];
design = [1*ones(size(tlHIT.trial,1),1); 2*ones(size(tlMISS.trial,1),1)];

clear tlHIT;
clear tlMISS;

% CSPPROCESSOR could be invoked independently of any CLFPROC pipeline:
% csp_analysis = cspprocessor('numpatterns',3,'outputdatatype','logpowcsp','filttype','biosig_CSP0','numchan',nchan);
% csp_analysis = csp_analysis.train(orig_data,design);
% data = csp_analysis.test(orig_data);

data = orig_data;
clear orig_data;

% ----------------   LDA classifier with CSP (logpow) features  --------------------
myproc = clfproc({ ... 
    cspprocessor('numchan',nchan,'numpatterns',3,'outputdatatype','logpowcsp','filttype','CSP0')...
    standardizer() ...
    da() ...
    });

% TRAINING error (or rather classification accuracy)
tr = crossvalidator('procedure',myproc,'cvfolds',1,'randomize',false,'verbose',true);
tr = tr.validate(data,design);
tracc_lda = tr.evaluate('metric','accuracy');

% 5-fold CROSSVALIDATION error (or rather classification accuracy)
cv = crossvalidator('procedure',myproc,'cvfolds',5,'randomize',true,'verbose',true);
M=1;  % for mutliple n-fold crossvalidation (here M=1 implies single partitioning)
acc = zeros(M,1);
for i=1:M
    cv = cv.validate(data,design);
    acc(i) = cv.evaluate('metric','accuracy');
end
disp('--------------   LDA results  -------------');
disp(sprintf('Crossvalidation ACC with = %5.4f +- %5.4f',mean(acc),std(acc)));
disp(sprintf('Training ACC with LDA    = %5.4f ',tracc_lda));
disp('-------------------------------------------');



%-------------   SVM classifier with CSP (logpow) features  ----------------
res = [];
tracc = [];

for C_svm = [ 1 3 5 10 20:20:100 200]   %different C settings are tried out in the crossvalidation analysis
    myproc = clfproc({ ...
        cspprocessor('numchan',nchan,'numpatterns',3,'outputdatatype','logpowcsp','filttype','CSP0')...
        standardizer() ...
        svmmethod('method',@svm_km_l2,'C',C_svm,'ratio4estplatt',0.0,'kernel','linear','kerparam',1) ...
        });

    % 5-fold CROSSVALIDATION error (or rather classification accuracy)
    cv = crossvalidator('procedure',myproc,'cvfolds',5,'randomize',true,'verbose',true);

    M=1;  %as above - 1 x 5-fold crossvalidation
    acc = zeros(M,1);
    for i=1:M
        cv = cv.validate(data,design);
        acc(i) = cv.evaluate('metric','accuracy');
    end

    disp(sprintf('Crossvalidation ACC = %5.4f +- %5.4f',mean(acc),std(acc)));

    % contains information about classification accuracy obtained in
    % crossvalidation for various values of hyperparameter C
    res = [res mean(acc)];

    % TRAINING error (or rather classification accuracy)
    tr = crossvalidator('procedure',myproc,'cvfolds',1,'randomize',false,'verbose',true);
    tr = tr.validate(data,design);
    tracc = [tracc tr.evaluate('metric','accuracy')];
    disp(sprintf('Training ACC = %5.4f  for C = %3.2f',tracc(end),C_svm));
end




