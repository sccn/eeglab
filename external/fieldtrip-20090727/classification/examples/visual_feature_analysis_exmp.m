
%%

orig_data = {MISSfft HITfft};

channels            = [1:151];

collective_cfg = {};
for i = channels
    cfg.channel            = channels(i);
    cfg.avgoverchan        = 'no';
    cfg.avgoverfreq        = 'yes';
    cfg.datarepresentation = 'concatenated';
    cfg.precision          = 'double';
    cfg.frequency          = [8 12];% This could also be set as foilim(Criterion_mapping(i),1:2);
                                    % It would have to be run subsequently to chan_discrim_analysis_exemp.m to
                                    % obtain the most discriminative frequency bands per channel
                                    % (purely discrimination-driven approach so it is better to be careful w.r.t.
                                    % to physiological underpinning)
    collective_cfg{i} = cfg;
end

[feature_data] = powspctrm_feature_extraction(orig_data,collective_cfg);

%  Data
design = feature_data.design(:,3);
data = feature_data.biol;

% the first and the second moment in each class
mHits = mean(data(design==2,:),1);
mMisses = mean(data(design==1,:),1);
sHits = std(data(design==2,:),1);
sMisses = std(data(design==1,:),1);

%  Normalized Data
ndata = normalizemeanstd(data);

% the first and the second moment in each class
n_mHits = mean(ndata(design==2,:),1);
n_mMisses = mean(ndata(design==1,:),1);
n_sHits = std(ndata(design==2,:),1);
n_sMisses = std(ndata(design==1,:),1);
%%

%----------------------------------------------------
%%

% Power spectral features extracted from any channels can be visualised,
% here: ch_group is adopted from a preceding run of chan_discrim_analysis_exemp.m (it is sorted according to FD_RATIO) 
ch_group = ch_rank;
          
% usually the top ranked channels are visualised  (1, 2 or 3)        
ch =  ch_group([1 2]);

if numel(ch)==1

    figure;
    hist(ndataMisses(:,ch),30); title('NdataMISS');
    figure;
    hist(ndataHits(:,ch),30);  title('NdataHIT');
    
    figure;
    hist(dataMisses(:,ch),30); title('dataMISS');
    figure;
    hist(dataHits(:,ch),30);  title('dataHIT');

elseif numel(ch)==2

    figure;
    plot(ndataMisses(:,ch(1)),ndataMisses(:,ch(2)),'or');
    hold on;
    plot(ndataHits(:,ch(1)),ndataHits(:,ch(2)),'*g');   title('Ndata');

    figure;
    plot(dataMisses(:,ch(1)),dataMisses(:,ch(2)),'or');
    hold on;
    plot(dataHits(:,ch(1)),dataHits(:,ch(2)),'*g');   title('data');
%     xlim([0,50e-27]);
%     ylim([0,50e-27]);

    figure;
    plot(mMisses(ch(1)),mMisses(ch(2)),'or');   hold on;
    plot([mMisses(ch(1)) - sMisses(ch(1))  mMisses(ch(1)) + sMisses(ch(1))],[mMisses(ch(2)) mMisses(ch(2))]); hold on;
    plot([mMisses(ch(1))   mMisses(ch(1))],[mMisses(ch(2))- sMisses(ch(2)) mMisses(ch(2))+ sMisses(ch(2))]);

    hold on;
    plot(mHits(ch(1)),mHits(ch(2)),'*g');   hold on;
    plot([mHits(ch(1)) - sHits(ch(1))  mHits(ch(1)) + sHits(ch(1))],[mHits(ch(2)) mHits(ch(2))]); hold on;
    plot([mHits(ch(1))   mHits(ch(1))],[mHits(ch(2))- sHits(ch(2)) mHits(ch(2))+ sHits(ch(2))]);  title('meand and std dev for data');
%     xlim([0,10e-28]);
%     ylim([0,10e-28]);

elseif numel(ch)==3
    
    scaling = 0.1;
    
    figure;
    
    ellipsoid(mMisses(ch(1)),mMisses(ch(2)),mMisses(ch(3)),scaling.*sMisses(ch(1)),scaling.*sMisses(ch(2)),scaling.*sMisses(ch(3)));
    hold on;
    ellipsoid(mHits(ch(1)),mHits(ch(2)),mHits(ch(3)),scaling.*sHits(ch(1)),scaling.*sHits(ch(2)),scaling.*sHits(ch(3)));
%     xlim([0,10e-28]);
%     ylim([0,10e-28]);
%     zlim([0,10e-28]);
         
else
    disp('The dimesnionality above 3 is not visualised');
end