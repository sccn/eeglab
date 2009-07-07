%% Using the classification module together with FieldTrip data
% This example demonstrates how to use neuroimaging data obtained from FieldTrip
% together with the classification module. Here, we have more control over what
% will happen compared with statistics_crossvalidate. 
% In the example, we make use of covert attention data of one subject that has 
% already been frequency analyzed. 
%
% The data consists of 7 different frequencies at 274 channels at time
% points [-0.5 0 0.5 1 1.5 2 2.5]. We can expect evoked response after the
% cue and alpha modulation after about 1 second.
%
% Copyright (C) 2008  Marcel van Gerven
%

%% Simple classification examples
function classification_demo()

%% 
  % get data and design (class labels)
  [data,design] = read_imaging_data();
 
  % concatenate the data into one dataset
  data = cat(1,data.biol{:});

%% 
% specify classification procedure

myproc = clfproc({ standardizer() lr() });

%%
% validation method; randomize trial order and give verbose output
cv = crossvalidator('procedure',myproc,'cvfolds',0.8,'randomize',true,'verbose',true);

%% 
% compute crossvalidation results
cv = cv.validate(data,design);

%% 
% get average classification accuracy
cv.evaluate('metric','accuracy')

%% 
% significance compared with a baseline classifier
cv.significance

