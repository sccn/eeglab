%% this script can be used to test your classifier

fclose('all');
close all
clear all
format long

% fix the RNG in order to reproduce the experiment
rand('twister',1); randn('state',1);

%% 
% load data
datas = cell(1,2);
load ~/data/alphalat/subject11/freqdat8144;
datas{1} = data;
load ~/data/alphalat/subject11/freqdat8142;
datas{2} = data;
data = datas;
clear datas;

%% 
% specify transformations of spectral data.
% Here, we may preselect only a few channels, times, and/or frequencies,
% and we may average over channels, time, or frequencies.
% we call a FieldTrip private function (this should be replaced at some
% point)
cfg.channel     = {'MLO' 'MRO'};
cfg.frequency   = [8 14];
cfg.latency     = [0.5 2.5];
cfg.avgoverfreq = 'yes';
cfg.avgovertime = 'yes';
[cfg,data] = prepare_timefreq_data(cfg,data{:});

%%
% create design matrix
design = [];
for c=1:length(data.biol)
    design = [design; c*ones(size(data.biol{c},1),1)];
end

% transform design to angles [-pi...pi]
design(design==1) = pi; % left
design(design==2) = 0; % right

% concatenate the data into one dataset
data = cat(1,data.biol{:});

%% 
% specify regression procedure: insert your method here!
%myproc = { standardizer() mulreg('prefun',@(x)([sin(x) cos(x)]),'postfun',@(x)(atan2(x(:,1),x(:,2))))  };
%myproc = { standardizer() circreg(); };
%myproc = { standardizer() linreg(); };
myproc = { standardizer() gpregressor(); };

%%
% validation method; randomize trial order and give verbose output
cv = crossvalidator('balanced',false,'procedure',myproc,'cvfolds',0.9,'randomize',true,'verbose',true);

% debugging
%data = [-ones(100,1); ones(100,1)] + 0.1*randn(200,1);
%design = [pi*ones(100,1); 0*ones(100,1)];% + 0.1*randn(200,1);

%% 
% compute crossvalidation results
cv = cv.validate(data,design);

%% 
% get residual sum of squares
cv.evaluate('metric','circrss')

%% 
% significance (F-test and t-test not yet implemented for regressors)
%cv.significance

%%
% plot regression
bar([cv.post cv.design])

if size(cv.post,2)>1
  
  bar([cv.design cv.post(:,1)]);
  hold on;
  errorbar(0.15+(1:size(cv.post,1)),cv.post(:,1),cv.post(:,2),'ko')

end