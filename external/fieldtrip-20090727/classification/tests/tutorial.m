% This Matlab script demonstrates how to use the neuroML toolbox
% together with FieldTrip in order to classify EEG/MEG data based on
% single trials. The dataset is already preprocessed and frequency analysed and
% consists of a cell array representing left motor execution (Lx), left motor
% imagery (Li), right motor execution (Rx), right motor imagery (Ri), and no motor
% response (No) classes.
%
% Assignments:
%
% 1. Study the script
% 2. Study the toolbox by calling: >> doc nmlt and checking the examples
%    folder.
% 3. Study the structure of the data file; this is the result of calling the
%    FieldTrip freqanalyse function on preprocessed data.
% 4. Select the classes you wish to compare and motivate your choice.
% 3. Try to find a subset of the (average over) channels, frequencies, and time
%    segments that gives good classification accuracy and explain this.
%    Hint: what prior knowledge do we have about this type of task?
% 4. Define your own classification procedure; how does this change the
%    results? Motivate your choice.
% 5. Choose another evaluation metric. Motivate your choice and explain the results.

% clear data and fix the random number generator (RNG)

fclose('all');
close all
clear all

format long

rand('twister',1); randn('state',1);

% load experimental data for a subject
data = load('~/data/christian/eegcharlotte/freqftrs_eeg_hjorth_test2.mat','dat'); data = data.dat;

% select conditions from {Rx, Ri, Lx, Li, No}
data = { data{1} data{2} data{3} data{4} data{5}}; 

% specify transformations of spectral data
% Here, we may preselect only a few channels, times, and/or frequencies,
% and we may average over channels, time, or frequencies.
cfg.channel     = {'F3' 'Fz' 'F4' 'FC5' 'FC1' 'FCz' 'FC2' 'FC6' 'C3' 'Cz' 'C4' 'CP5' 'CP1' 'CP2' 'CP6' 'P3' 'Pz' 'P4'}; % also see CHANNELSELECTION
cfg.frequency   = [1 150]; % a range in data{i}.freq
cfg.avgoverchan = 'no';
cfg.avgoverfreq = 'no';
cfg.avgovertime = 'no';

% perform FieldTrip transformations
% possible other transformations may be specified here
% (baseline correction etc)
[cfg,data] = prepare_timefreq_data(cfg,data{:});

% create design matrix; this specifies at least the class labels (1,2,3,...) for each example
design = [];
for c=1:length(data.biol)
    design = [design; c*ones(size(data.biol{c},1),1)];
end

% concatenate the data into one dataset
data = cat(1,data.biol{:});

% specify classification procedure
myproc = clfproc({ ...
    preprocessor('prefun',@(x)(log10(x))) ...
    standardizer() ...
    nb() ...
    });

% validation method
cv = crossvalidator('procedure',myproc,'cvfolds',10,'randomize',true,'verbose',true);

% compute crossvalidation results
cv = cv.validate(data,design);

% get average classification accuracy
cv.evaluate('metric','accuracy')







%
% Short introduction:
% 
% A typical machine learning task relies on the application of a number of 
% transformations to a dataset of interest in order to obtain some prediction.
% The main focus is on classification.
%
% Suppose we have an N x M dataset with N examples and M features for which 
% we wish to predict class labels C as captured in a design matrix (typically
% a vector containing the labels. Suppose furthermore that we wish to standardize 
% the data and perform a simple discriminant analysis using ten-fold cross-
% validation. This is as simple as:
% 
% >> cv = crossvalidator(clfproc({ standardizer() da()}));
% >> cv = cv.validate(data,design);
% 
% after which cv.post contains the posterior probabilities of the class labels
% and cv.design the associated design for the 10-fold crossvalidation.
% 
% In order to compute the classification accuracy, we call evaluate:
% 
% >> cv.evaluate('metric','accuracy')
% 
%  ans =
% 
%    0.857142857142857
%
%  Multiple datasets are entered using cell arrays of data/design. This is
%  useful when analyzing multiple datasets, for ensemble learning, and for
%  transfer learning.
%
%  In case of online use, we create a classification procedure 
%
%  >> myproc = clfproc({ standardizer() da()})
%
%  train it with some data
%
%  >> myproc = myproc.train(data,design)
%
%  and test it by entering a data vector
%
%  >> myproc.test(datavec)
%
%  Sometimes, it is useful to compare two algorithms with one another.
%  E.g., in order to compute whether or not two algorithms behave
%  differently, we can call
%
%  >> [reject,pvalue,level] = significance(v1.post,v1.design,v2.post,v2.design,varargin)
%
%  for two validators v1 and v2.
%
%  More extensive examples are given in the examples directory or by
%  typing 'doc method' for your method of interest.
%





%  other options
% featureselector('subset',1) ...
% ensemble('procedures',{ clfproc({standardizer() da()}) clfproc({standardizer() kernelmethod()}) },'combination','voting') ...
% pnn('spread',[0.01 0.1 1.0 2.0]) ...
% ensemble('classifiers',{da() l2svm()},'combination','voting');
% filterer('filter',@mutual_information,'verbose',true,'validator',crossvalidator('procedure',clfproc({da()}),'cvfolds',0.9)) ...
% pcanalyzer('proportion',0.8) ...
% one_against_one('procedures',da())
% ensemble('procedures',da())
% kernelmethod('method',@l2svm_cg,'C',10:10:100)









