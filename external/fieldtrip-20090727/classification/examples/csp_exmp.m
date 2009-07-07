% csp_examp.m

% This script demonstrates how CSP patterns can be extracted from two-class (two-condition) data
% outside the scope of classification toolbox (using typical Fieldtrip syntax without an OOP approach)

%preprochits and preprocmisses are two data sets

%%
cfg.channel = 'all';
cfg.bpfilter = 'yes';
cfg.bpfiltdir = 'twopass'; 
cfg.bpfilttype = 'but'; 
cfg.detrend = 'yes'; 
cfg.datatype = 'continuous'; 
cfg.keeptrials = 'yes'; 
cfg.bpfreq = [8 12];

% preprocessing to filter the signals within a specified freq band (here 8-12 Hz)
filt_data1 = preprocessing(cfg,preprochits);
filt_data2 = preprocessing(cfg,preprocmisses);

% design matrix containing labels
design = [1 * ones(length(filt_data1.trial),1); ...
          2 * ones(length(filt_data2.trial),1)];
%% 

%%
%  TIMELOCKANALYSIS for VARIANTS II and III
cfg = [];
cfg.keeptrials = 'yes';
tlpreprochits = timelockanalysis(cfg, filt_data1);
tlpreprocmisses = timelockanalysis(cfg, filt_data2);
%%


%%
%  TIMELOCKANALYSIS on appended data: VARIANT I 
tlboth = timelockanalysis(cfg, appenddata([], filt_data1, filt_data2));

[filters1,d1] = csp_train(tlboth,design,2,'fieldtrip_cov');      
[csp_proj1,csp_pow1] = csp_test(tlboth,filters1);
%%

%%
% Explicit appending the data sets obtained using TIMELOCKANALYSIS - VARIANT II
data_both = appenddata([],tlpreprochits,tlpreprocmisses);

[filters2,d2] = csp_train(data_both,design,2,'biosig_CSP0'); 
[csp_proj2,csp_pow2] = csp_test(data_both,filters2);
%%

%%
% Concatenating the data sets obtained using TIMELOCKANALYSIS - VARIANT III
data_both_mat = [tlpreprochits.trial; tlpreprocmisses.trial];

[filters3,d3] = csp_train(data_both_mat,design,2,'EED');      
[csp_proj3,csp_pow3] = csp_test(data_both_mat,filters3);
%%
