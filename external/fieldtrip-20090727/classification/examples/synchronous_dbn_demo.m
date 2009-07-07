%% Synchronous DBNs for classification
%
% Synchronous DBNs can be used to model cue-based trial data. Hence, data
% with a finite horizon.
%
% Synchronous temporal data has dimensions N x (M x T) with N examples, M features, and 
% T time slices. In a synchronous experiment, each example is a finite horizon dynamic
% model of which the underlying state is assumed to be constant. In
% general, the data should have time as the last dimension
%
% The data consists of 7 different frequencies at 274 channels at time
% points [-0.5 0 0.5 1 1.5 2 2.5]. We can expect evoked response after the
% cue and alpha modulation after about 1 second.
%
% Copyright (C) 2008  Marcel van Gerven
%

%% Simple classification examples
function synchronous_dbn_demo()

%% 
% get data and design (class labels)
[data,design] = read_imaging_data();

% concatenate the data into one dataset
data = cat(1,data.biol{:});

%% 
% specify classification procedure

myproc = clfproc({ ...    
    preprocessor('prefun',@(x)(log10(x))) ...
    hmm('horizon',5,'coupled',true,'ar',false,'mixture',1,'verbose',true) ...    
    });

% validation method
cv = crossvalidator('procedure',myproc,'cvfolds',0.9,'randomize',true,'verbose',true);

%% 
% compute crossvalidation results
cv = cv.validate(data,design);

%% 
% get average classification accuracy
cv.evaluate('metric','accuracy')

%% 
% significance compared with a baseline classifier
cv.significance











