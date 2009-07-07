function [data,design] = read_imaging_data()
    % get data used in examples
    
    % clear data and fix the random number generator (RNG)

    fclose('all');
    close all
    clear all
    format long

    % fix the RNG in order to reproduce the experiment
    rand('twister',1); randn('state',1);

    % load data
    load ~/code/classification/toolboxes/bayesbrain/examples/freqli;
    load ~/code/classification/toolboxes/bayesbrain/examples/freqri;
    data = {freqLI freqRI}; clear freqLI; clear freqRI;

    % specify transformations of spectral data.
    % Here, we may preselect only a few channels, times, and/or frequencies,
    % and we may average over channels, time, or frequencies.

    cfg.channel     = {'MLO' 'MRO'};
    cfg.frequency   = [8 14];    
    cfg.latency     = [0.5 2.5];
    cfg.avgoverfreq = 'yes';

    % we call a FieldTrip private function (should be replaced at some point)
    [cfg,data] = prepare_timefreq_data(cfg,data{:});
    
    % create design matrix (class labels)
    design = [];
    for c=1:length(data.biol)
        design = [design; c*ones(size(data.biol{c},1),1)];
    end

   
end