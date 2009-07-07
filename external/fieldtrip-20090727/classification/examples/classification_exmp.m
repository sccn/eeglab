% classification_exmp.m

% This script presents a classification approach to data analysis. It
% should be treated as an example.

addpath(genpath('h:\common\matlab\fieldtrip\classification'));

inp.isplanar = 'yes';
inp.whatData = 'JM';
inp.mtmtype = 'mtmfft';

% load the data set: here two data sets were analysed in the spectral
% domain using FREQANALYSIS
basicPATH = strcat('D:\new4Pawel\',inp.whatData,'\');
if strcmp(inp.isplanar,'yes')
    in_datafile_1 = 'HITfft_plan_mtmfft.mat';
    in_datafile_2 = 'MISSfft_plan_mtmfft.mat';
else
    in_datafile_1 = 'HITfft_mtmfft.mat';
    in_datafile_2 = 'MISSfft_mtmfft.mat';
end
load(cat(2,basicPATH,in_datafile_1),'HITfft');
load(cat(2,basicPATH,in_datafile_2),'MISSfft');


disp(' ');
disp(sprintf('%3d hits',size(HITfft.powspctrm,1)));
disp(sprintf('%3d misses',size(MISSfft.powspctrm,1)));

if strcmp(inp.whatData,'JM')
    disp('*****   REMOVAL   ****** ');
    rmHIT  = [5 7 14 15 19 21 22 26 28 30 32 35 39 40 47 65 75 83 85 89 90 94 96 98 108 116 117 119 133 136 139 141 147 160 161];
    rmMISS = [19 41 44 45 46 48 49 50 51 52 54 58 59 66 67 68 70 73 78 80 120 121 123 124 125 126 128 130 132 133 134 136 82 84 87 89 92 93 95 97 105 ...
        107 108 110 116 118 141 143 144 145 150 153 155 159 161 166 171 172 173 174 175 180 ];
    HITfft.powspctrm = HITfft.powspctrm(setdiff(1:size(HITfft.powspctrm,1),rmHIT),:,:);
    MISSfft.powspctrm = MISSfft.powspctrm(setdiff(1:size(MISSfft.powspctrm,1),rmMISS),:,:);
    HITfft.cumsumcnt = HITfft.cumsumcnt(setdiff(1:size(HITfft.cumsumcnt,1),rmHIT));
    MISSfft.cumsumcnt = MISSfft.cumsumcnt(setdiff(1:size(MISSfft.cumsumcnt,1),rmMISS));
    HITfft.cumtapcnt = HITfft.cumtapcnt(setdiff(1:size(HITfft.cumtapcnt,1),rmHIT));
    MISSfft.cumtapcnt = MISSfft.cumtapcnt(setdiff(1:size(MISSfft.cumtapcnt,1),rmMISS));

    disp(sprintf('%3d hits after removal',size(HITfft.powspctrm,1)));
    disp(sprintf('%3d misses after removal',size(MISSfft.powspctrm,1)));
        
end

if strcmp(inp.whatData,'JM')  
    % REMARK: these channels do not necessarily represent optimal channels for the given subject - just an example
    selected_channels{1} = {'MLO11','MLO12','MLO21','MLO22','MLO32','MLP11','MLP21','MLP22','MLP31','MLP32','MZO','MZP','MRO11','MRO12','MRO21','MRO22','MRO32'};
    selected_channels{2} = {'MRP11','MRP21','MRP22','MRP31','MRP32','MRT16','MRT26'};       
    selected_freq = [8 12; 9 10];    

else
    selected_channels = 'all';
    selected_freq = [8 12];
end

stat_data = {MISSfft HITfft}; %data to extract features from
cfgs = {};
for i = 1:size(selected_freq,1)
    cfg.channel            = selected_channels{i};
    cfg.frequency          = [8 12];
    cfg.avgoverchan        = 'no';
    cfg.avgoverfreq        = 'yes';
    cfg.datarepresentation = 'concatenated';
    cfg.precision          = 'double';
  
    cfgs{i} = cfg; 
end

[feature_data] = powspctrm_feature_extraction(stat_data,cfgs);

% preparation of the final feature data set for classification
design = feature_data.design(:,3);
data = feature_data.biol ;

%++++++++++++++++++++++++++++++++++++++++++++
% ----------   Classification  --------------
%++++++++++++++++++++++++++++++++++++++++++++

applyLDA = 0;
applySVM = 0;
applyLVQ = 1;

if applyLDA
    %-------------------------------------------
    %-------   LDA classification   ------------
    %-------------------------------------------
    myproc = clfproc({ ...
        %preprocessor('prefun',@(x)(log10(x))) ...
        standardizer() ...
        da() ...
        });

    % TRAINING error (or rather classification accuracy)
    tr = crossvalidator('procedure',myproc,'cvfolds',1,'randomize',false,'verbose',true);
    tr = tr.validate(data,design);
    tracc_lda = tr.evaluate('metric','accuracy');

    % 5-fold CROSSVALIDATION error (or rather classification accuracy)
    cv = crossvalidator('procedure',myproc,'cvfolds',5,'randomize',true,'verbose',true);
    M=10;
    acc = zeros(M,1);
    for i=1:M
        cv = cv.validate(data,design); % for mutliple n-fold crossvalidation (here 10x5-fold crossvalidation is performed)
        acc(i) = cv.evaluate('metric','accuracy');
    end
    disp('--------------   LDA results  -------------');
    disp(sprintf('Crossvalidation ACC with = %5.4f +- %5.4f',mean(acc),std(acc)));
    disp(sprintf('Training ACC with LDA    = %5.4f ',tracc_lda));
    disp('-------------------------------------------');

end


if applySVM
    %-------------------------------------------
    % ---------  SVM classification   ----------
    %-------------------------------------------

    res_svm = [];
    tracc_svm = [];

    % SVM's hyperparameter selection in 10x5-fold crossvalidation (see below)
    for C_svm = [1 3 5 10 20:20:100 200]
        myproc = clfproc({ ...
            %   preprocessor('prefun',@(x)(log10(x))) ...
            standardizer() ...
            svmmethod('method',@svm_km_l2,'C',C_svm,'ratio4estplatt',0.0,'kernel','linear','kerparam',1) ...
            });

        % 5-fold CROSSVALIDATION error (or rather classification accuracy)
        cv = crossvalidator('procedure',myproc,'cvfolds',5,'randomize',true,'verbose',true);

        M=10;
        acc = zeros(M,1);
        for i=1:M  % for mutliple n-fold crossvalidation (here 10x5-fold crossvalidation is performed)
            cv = cv.validate(data,design);
            acc(i) = cv.evaluate('metric','accuracy');
        end

        cvacc_svm = mean(acc);
        disp(sprintf('Crossvalidation ACC = %5.4f +- %5.4f',cvacc_svm,std(acc)));

        % contains information about classification accuracy obtained in
        % crossvalidation for various values of hyperparameter C (if a non-linear
        % kernel is used with SVM,e.g. gaussian, then additional kernel parameters
        % can also be selected in the same crossvalidation analysis)
        res_svm = [res_svm cvacc_svm];

        % TRAINING error (or rather classification accuracy)
        tr = crossvalidator('procedure',myproc,'cvfolds',1,'randomize',false,'verbose',true);
        tr = tr.validate(data,design);
        tracc_svm = [tracc_svm tr.evaluate('metric','accuracy')];
        disp(sprintf('Training ACC = %5.4f  for C = %3.2f',tracc_svm(end),C_svm));
    end

end


if applyLVQ
    %-------------------------------------------
    % ---------  LVQ classification   ----------
    %-------------------------------------------

    % Different variants of LVQ object are stored in a cell array:
    lvq_objs = {lvq('method','lvq1','numclasses',2,'initmethod','randinit1','initparam',20),...
        lvq('method','lvq3','numclasses',2,'initmethod','randinit2','initparam',20),...
        lvq('method','dslvq','numclasses',2,'initmethod','kmeansinit','initparam',[4 0]),...
        lvq('method','dslvq','numclasses',2,'initmethod','fcminit','initparam',[5 0]),...
        lvq('method','dslvq','numclasses',2,'initmethod','subclustinit','initparam',[0.5 1])};

    lvq_results = zeros(length(lvq_objs),2);   % two columns to store mean crossvalidation accuracy and training accuracy

    for variant=1:length(lvq_objs)   

        myproc = clfproc({ ...
            standardizer() ...
            lvq_objs{variant}});

        %5-fold CROSSVALIDATION error (or rather classification accuracy)
        cv = crossvalidator('procedure',myproc,'cvfolds',5,'randomize',true,'verbose',true);

        M=5; 
        acc = zeros(M,1);
        for i=1:M  % for mutliple n-fold crossvalidation (here 5x5-fold crossvalidation is performed)
            cv = cv.validate(data,design);
            acc(i) = cv.evaluate('metric','accuracy');
        end

        cvacc_lvq = mean(acc);
        disp(sprintf('Crossvalidation ACC = %5.4f +- %5.4f',cvacc_lvq,std(acc)));


        %TRAINING error (or rather classification accuracy)
        tr = crossvalidator('procedure',myproc,'cvfolds',1,'randomize',false,'verbose',true);
        tr = tr.validate(data,design);
        tracc_lvq = tr.evaluate('metric','accuracy');
        disp(sprintf('Training ACC = %5.4f  for variant %d',tracc_lvq,variant));

        lvq_results(variant,1) = cvacc_lvq;
        lvq_results(variant,2) = tracc_lvq;

    end

end
