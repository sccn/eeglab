%% this script can be used to test your classifier

fclose('all');
close all
clear all
format long

% fix the RNG in order to reproduce the experiment
rand('twister',1); randn('state',1);

% debugging; generate data
debug = 0;
if debug
    data = 0.1*randn(1000,2);

    mycr = circreg();
    mycr.mu = 0;
    mycr.kappa = 1;
    mycr.beta = [1 10]';
    mycr.gamma = [1 10]';

    design = mycr.sample(data);
    
    myproc = { circreg('mode','mixed','verbose',true); };

else
    %%
    % load data for 4 directions
    datas = cell(1,4);
    load ~/data/alphalat/subject11/freqdat8141;
    datas{1} = data;
    load ~/data/alphalat/subject11/freqdat8142;
    datas{2} = data;
    load ~/data/alphalat/subject11/freqdat8143;
    datas{3} = data;
    load ~/data/alphalat/subject11/freqdat8144;
    datas{4} = data;
    data = datas; clear datas;

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
    design(design==1) = pi/2; % up
    design(design==2) = 0; % right
    design(design==3) = -pi/2; % down
    design(design==4) = pi; % left

    % concatenate the data into one dataset
    data = cat(1,data.biol{:});
    
    myproc = { standardizer() circreg('mode','mean','verbose',true); };
end

%%
% validation method; randomize trial order and give verbose output
cv = crossvalidator('balanced',false,'procedure',myproc,'cvfolds',10,'randomize',true,'verbose',true);

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
% plot regression: DIFFERENCE BETWEEN ANGLES
% met = mod(mod(cv.design,2*pi) - mod(cv.post,2*pi),2*pi);
% met = min(met,2*pi - met);
% bar(met)
