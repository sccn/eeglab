

% chan_discrim_analysis_exemp.m

stat_data = {MISSfft HITfft}; 

fc = []; diffc = []; nfc = []; ndiffc = [];

foilim = [9 10; 8 10; 8 12; 9 12; 8 11; 9 11];


% searching through the frequency bands (among those in rows of foilim array)
% to select the most discriminative information in each individual channel
for i=1:size(foilim,1)

    cfg.channel            = [1:151];
    cfg.frequency          = foilim(i,1:2);
    cfg.avgoverchan        = 'no';
    cfg.avgoverfreq        = 'yes';
    cfg.datarepresentation = 'concatenated';
    cfg.precision          = 'double';

    [feature_data] = powspctrm_feature_extraction(stat_data,{cfg});

    % FD_RATIO (Fischer's disriminant ratio and the absolute difference)
    % calculated for original data - channel-wise
    [fc(i,:) diffc(i,:)] = fd_ratio(feature_data.biol,feature_data.design(:,3));

end

% criterion for channel ranking - either fc or diffc
%*****************************************************
Criterion = fc';
%+++++++++++++++++++++++++++++++++++++++++++++++++++++
[Criterion Criterion_mapping ] = max(Criterion,[],2);
%*****************************************************

cfg = [];
cfg.layout 	= 'CTF151.lay';
cfg.colorbar = 'yes';
cfg.electrodes = 'labels';
figure;
topoplot(cfg,Criterion);  % criterion is plotted here 
%(the higher the value the most discriminative information is offered at a given channel - from FD_RATIO's perspective)

[sortCriterion ch_rank] = sort(Criterion',2,'descend');

tracc_rfd = [];
tracc_lda = [];
cv_rfd = [];
cv_lda = [];
ch_group = [];

% evaluation of the feature separability between two conditions (quantified
% with classification accuracy - LDA and RFD are applied) for different
% channel numbers: 1,2,3,...nchan (in each iteration the top ranked channels are selected)
for ch = 1:length(Criterion)
           
    ch_group = [ch_group ch_rank(ch)];
    
    cfg.channel            = ch_group;
    cfg.frequency          = foilim(Criterion_mapping(ch_rank(ch)),1:2); % channel-specific freq band
    cfg.avgoverchan        = 'no';
    cfg.avgoverfreq        = 'yes';
    cfg.datarepresentation = 'concatenated';
    cfg.precision          = 'double';
   
    [feature_data] = powspctrm_feature_extraction(stat_data,{cfg});

    % data and design matrix
    design = feature_data.design(:,3);
    data = feature_data.biol;

    
    
    %--------    RFD evaluation  ------------
    
    res = []; % for different values of regularisation C for RFD

    N_iter=10;  % multiple 5-fold crossvalidation (see classification_examp.m) for parameter selection 
    % and estimation of generalisation error (or classification accuracy)
    for C_svm = [1:2:10 20:20:100]
        myproc = clfproc({ ...
            standardizer() ...
            rfda('C',C_svm) ...
            });

        %*****   CROSSVALIDATION accuracy  **********
        cv = crossvalidator('procedure',myproc,'cvfolds',5,'randomize',true,'verbose',false);
        acc = zeros(N_iter,1);
        for i=1:N_iter
            cv = cv.validate(data,design);
            acc(i) = cv.evaluate('metric','accuracy');
        end

        disp(sprintf('Crossvalidation with RFD - ACC = %5.4f +- %5.4f',mean(acc),std(acc)));

        res = [res; mean(acc)];
    end
    
    cv_rfd = [cv_rfd res];

    %*****   TRAINING accuracy  **********
    cv = crossvalidator('procedure',myproc,'cvfolds',1,'randomize',true,'verbose',false);
    cv = cv.validate(data,design);
    tracc_rfd = [ tracc_rfd cv.evaluate('metric','accuracy')];
    disp(sprintf('RFD Training ACC = %5.4f',tracc_rfd(ch)));

    
    %--------   LDA evaluation ---------
    myproc = clfproc({ ...
        standardizer() ...
        da() ...
        });
    
    %*****   CROSSVALIDATION accuracy  **********
    cv = crossvalidator('procedure',myproc,'cvfolds',5,'randomize',true,'verbose',false);
    acc_lda = zeros(N_iter,1);
    for i=1:N_iter
        cv = cv.validate(data,design);
        acc_lda(i) = cv.evaluate('metric','accuracy');
    end

    disp(sprintf('Crossvalidation with LDA - ACC = %5.4f +- %5.4f',mean(acc),std(acc)));

    res_lda = mean(acc_lda);
    cv_lda = [cv_lda res_lda];

    %*****   TRAINING accuracy  **********
    cv = crossvalidator('procedure',myproc,'cvfolds',1,'randomize',true,'verbose',false);
    cv = cv.validate(data,design);
    tracc_lda = [tracc_lda  cv.evaluate('metric','accuracy')];
    disp(sprintf('LDA Training ACC = %5.4f',tracc_lda(ch)));

    % information about the progress
    disp('-----------------------------------------------');
    disp(sprintf('---      Channel config %3d out of %3d      ---',ch,length(Criterion)));
    disp('-----------------------------------------------');
    
end

figure;
plot(tracc_lda,'b--'); hold on;
plot(cv_lda,'b-'); hold on;

plot(tracc_rfd,'r--'); hold on;
plot(cv_rfd,'r-');
title('Classification Accuracy for training (dotted) and crossvalidation (solid) using LDA (blue) and RFD (red)');

