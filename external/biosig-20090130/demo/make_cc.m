% Make classifier - Demo how to perform offline BCI analysis and for
% obtaining a classifier for online feedback. 
%
% Copyright (C) 2006, Alois Schloegl 
% 
% References: 
% [1] A. Schlögl, K. Lugger and G. Pfurtscheller (1997) 
% Using Adaptive Autoregressive Parameters for a Brain-Computer-Interface Experiment, Proceedings of the 19th Annual International Conference if the IEEE Engineering in Medicine and Biology Society ,vol 19, pp.1533-1535, 1997 http://hci.tugraz.at/schloegl/publications/schloegl1997b.pdf	
% 
% [2] A. Schlögl, C. Neuper G. Pfurtscheller (2002), 
% Estimating the mutual information of an EEG-based Brain-Computer-Interface. Biomedizinische Technik 47(1-2): 3-8, 2002 http://hci.tugraz.at/schloegl/publications/schloegl2002.pdf		
% 
% [3] A. Schlögl (2006) 
% GDF - A general dataformat for Biosignals (Version 2.0). - http://arxiv.org/abs/cs.DB/0608052 .
% 
% [4] A. Schlögl, C. Brunner, R. Scherer, A. Glatz; 
% BioSig - an open source software library for BCI research. Towards Brain-Computer Interfacing, (Eds.) G. Dornhege, J.R. Millan, T. Hinterberger, D.J. McFarland, K.-R.Müller; MIT press (accepted).
% 
% [5] A. Schlögl, J. Kronegg, J.E. Huggins, S. G. Mason;
% Evaluation criteria in BCI research. 
% (Eds.) G. Dornhege, J.R. Millan, T. Hinterberger, D.J. McFarland, K.-R.Müller; Towards Brain-Computer Interfacing, MIT press (accepted) 

%% load data 
[s,HDR] = sload('l1.gdf');

%% artifact correction, preprocessing 
% -none-

%%--%% feature extraction 
%[A,M,C] = hjorth(s,HDR.SampleRate); f = [A,M,C];
%[A,M,C] = barlow(s,HDR.SampleRate); f = [A,M,C];
%[A,M,C] = wackermann(s,HDR.SampleRate); f = [A,M,C];
f = bandpower(s,HDR.SampleRate,[10,12;16,24],1);
f(f<-10)=NaN; % flat lines can result in -infinity, ignore these segments 
% Adaptive Autoregressive parameters 
%a = [];
%for k=1:size(s,2),
%	X = tvaar(s(:,k),3,0.006);
%	X = tvaar(s(:,k),X);
%	a = [a,X.AAR];
%end;

%% identify and validate classifier 
seg  = reshape(1:8*HDR.SampleRate,HDR.SampleRate/5,8*5)'; % define segments, typical 1/10s .. 1/5s 
flag = seg(:,1)>(3*HDR.SampleRate);     % flag of possible segments
ng   = ceil([1:length(HDR.TRIG)]'/length(HDR.TRIG)*K);	% K_fold XV

CC   = findclassifier(f, HDR.TRIG, [HDR.Classlabel,ng], seg, flag, 'LDA');

MODE.TYPE = 'SVM';
MODE.hyperparameters.c_value = 1e4;
%CC   = findclassifier(f, HDR.TRIG, [HDR.Classlabel,ng], seg, flag, MODE);

MODE.TYPE = 'RBF';
MODE.hyperparameters.c_value = 1e4;
MODE.hyperparameters.gamma = 1000;
%CC   = findclassifier(f, HDR.TRIG, [HDR.Classlabel,ng], seg, flag, MODE);

%% display result 
CC.TSD.T = CC.TSD.T/HDR.SampleRate; % scale time axis 
plota(CC.TSD)

% CC.weights containts the weights of the classier, (one vector for each class)
% CC.weigths; % contains the classifier for online feedback, i.e. weights for the linear combiner. 
% CC.weigths(1,:) is the bias and CC.weights(2:end,:) are the weights for each feature 
% CC.TI contains time segment use for the final classifier 
% CC.TSD.{ACC00, KAP00, ...} contain time course of performance metric (cross-validated)

