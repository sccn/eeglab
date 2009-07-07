%% this script can be used to test your classifier

fclose('all');
close all
clear all
format long

% fix the RNG in order to reproduce the experiment
rand('twister',1); randn('state',1);

%% 
% load data
load ~/code/classification/toolboxes/bayesbrain/examples/freqli;
load ~/code/classification/toolboxes/bayesbrain/examples/freqri;
data = {freqLI freqRI}; clear freqLI; clear freqRI;

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

% concatenate the data into one dataset
data = cat(1,data.biol{:});

%% 
% specify classification procedure: insert your method here!

% define groups
sz = size(data);
groups = cell(1,prod(sz(2:end)));
for c=1:length(groups)
    groups{c} = [2*(c-1)+1 2*(c)];
end

% gamma / exponential does best

ii=2;
if ii==1

    % custom pdf
    truncgauss = @(x,mu,sigma)(normpdf(x,mu,sigma)./(normcdf(500,mu,sigma)-normcdf(0,mu,sigma)));
    
    % custom mle
    truncmle = @(x)(mle(x,'pdf',truncgauss,...
        'start', [nanmean(x) nanstd(x)], ...
        'lower', [0 0]));

    myproc = { ...
        gnb('conditional',truncgauss,'mle',truncmle)...
        };
    
elseif ii==2

    % custom pdf
    vonmises = @(x,mu,k)((1/(2*pi*besseli(0,k)))*exp(k*cos(mod(x,2*pi)-mu))); %gives the density at that point

    myproc = { ...
        gnb('conditional',vonmises,'mle',@mlemises)...
        };
    
elseif ii==3
    
    myproc = { gnb('conditional','exponential') };
end
    
%%
% validation method; randomize trial order and give verbose output
cv = crossvalidator('balanced',false,'procedure',myproc,'cvfolds',0.9,'randomize',true,'verbose',true);

%% 
% compute crossvalidation results

cv = cv.validate(data,design);

%% 
% get average classification accuracy

cv.evaluate('metric','accuracy')

%% 
% significance compared with a baseline classifier
cv.significance

