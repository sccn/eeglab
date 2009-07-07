
 % channel_correl_analysis_exmp.m
 
 % This cript presents an example of correlation analysis on the power
 % spectral features extracted form the data - the emphasis is on how the
 % top ranked channels (using FR_RATIO criterion) correlate with the others
 
 
% ch_rank - FD_RATIO-based ranking that can be obtaiend by running the first part of chan_discrim_analysis_exmp.m (prerequisite variable)


% noptchan - number of optimal top rank channels to analyse
  noptchan = 4;
  chan_foilim = [8 12];  %it can be generated using channel_dscrim_analysis
  
% extract relevant features - here power spectral content from each channel (independent freq bands)  
cfgs = {};
if size(chan_foilim,1) == 1
    cfg.channel            = 'all';
    cfg.frequency          = chan_foilim;
    cfg.avgoverchan        = 'no';
    cfg.avgoverfreq        = 'yes';
    cfg.datarepresentation = 'concatenated';
    cfg.precision          = 'double';
    cfgs = {cfg};
elseif size(chan_foilim,1) > 1   %this size should actually be equal to the number of channels
    for i=1:size(chan_foilim,1)
        cfg.channel            = i;
        cfg.frequency          = chan_foilim(i,1:2);
        cfg.avgoverchan        = 'no';
        cfg.avgoverfreq        = 'yes';
        cfg.datarepresentation = 'concatenated';
        cfg.precision          = 'double';
        cfgs{i} = cfg;
    end
end

% power spectral feature extraction
[feature_data] = powspctrm_feature_extraction(stat_data,cfgs);
data = feature_data.biol;
design = feature_data.design(:,3);

% estimate different quantifiers of feature-feature and feature-class correlation
% (see function FEATURE_CORR.m)
[C_bpc,C_wfc,C_ff,C_ff1,C_ff2] = feature_corr(data,design);

cfg = [];
cfg.layout 	= 'CTF151.lay';
cfg.colorbar = 'yes';
cfg.hlmarker   = 'x';
cfg.hlcolor    = [0 0 0];
cfg.hlmarkersize = 8;

for i=1:noptchan
    
    cfg.highlight =  ch_rank(i);
    
    figure;
    % correlation between one of the top rank channels and the others
    % (feature correlation)
    subplot(3,1,1);
    topoplot(cfg,(C_ff(ch_rank(i),:)'));
    title(sprintf('Corr f-f         RANK %d   ch %d',i));
    
    % correlation between one of the top rank channels and the others
    % (class-specific feature correlation, i.e. data representing class 1 are only taken into account)
    subplot(3,1,2);
    topoplot(cfg,(C_ff1(ch_rank(i),:)'));
    title(sprintf('Corr f-f class1  RANK %d   ch %d',i));

    % correlation between one of the top rank channels and the others
    % (class-specific feature correlation, i.e. data representing class 2 are only taken into account)
    subplot(3,1,3);
    topoplot(cfg,(C_ff2(ch_rank(i),:)'));
    title(sprintf('Corr f-f class2  RANK %d   ch %d',i));
end

cfg.hlmarker   = 'o';
cfg.highlight = [ch_rank(1:noptchan)];
cfg.electrodes = 'labels';

figure;
topoplot(cfg,C_bpc');
title('Bi-serial point correlation between channels and conditions (classes)');

