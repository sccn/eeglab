
% DEMO2 demonstrates the use of the data set III from the BCI competition 2003 for 
%   The demo shows the offline analysis for obtaining a classifier and 
%   uses a jack-knife method (leave-one-trial out) for validation. 
%   AAR parameters are extracted 
%
%
% References: 
% [1] A. Schlögl, B. Kemp, T. Penzel, D. Kunz, S.-L. Himanen, A. Värri, G. Dorffner, G. Pfurtscheller.
%       Quality Control of polysomnographic Sleep Data by Histogram and Entropy Analysis.
%       Clin. Neurophysiol. 1999, Dec; 110(12): 2165 - 2170.
% [2] Alois Schlögl (2000)
%       The electroencephalogram and the adaptive autoregressive model: theory and applications
%       Shaker Verlag, Aachen, Germany, (ISBN3-8265-7640-3). 
% [3] Schlögl A., Neuper C. Pfurtscheller G.
%       Estimating the mutual information of an EEG-based Brain-Computer-Interface
%       Biomedizinische Technik 47(1-2): 3-8, 2002.
% [4] A. Schlögl, C. Keinrath, R. Scherer, G. Pfurtscheller,
%       Information transfer of an EEG-based Bran-computer interface.
%       Proceedings of the 1st International IEEE EMBS Conference on Neural Engineering, Capri, Italy, Mar 20-22, 2003. 
% [5] A. Schlögl, J. Kronegg, J.E. Huggins, S. G. Mason.
%       Evaluation criteria in BCI research.
%       (Eds.) G. Dornhege, J.R. Millan, T. Hinterberger, D.J. McFarland, K.-R.Müller,
%       Towards Brain-Computer Interfacing. MIT Press, p.327-342, 2007.
% [6] A. Schlögl, F.Y. Lee, H. Bischof, G. Pfurtscheller
%   	Characterization of Four-Class Motor Imagery EEG Data for the BCI-Competition 2005.
%   	Journal of neural engineering 2 (2005) 4, S. L14-L22
% [7] A. Schlögl, C. Brunner, R. Scherer, A. Glatz;
%   	BioSig - an open source software library for BCI research.
%   	(Eds.) G. Dornhege, J.R. Millan, T. Hinterberger, D.J. McFarland, K.-R. Müller;
%   	Towards Brain-Computer Interfacing, MIT Press, 2007, p.347-358. 
% [8] A. Schlögl, C. Brunner
%   	BioSig: A Free and Open Source Software Library for BCI Research.
%	Computer (2008, In Press)	

%	$Id: demo2.m,v 1.1 2009-01-30 06:04:39 arno Exp $
%	Copyright (C) 1999-2003,2006,2007,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'B0101T.gdf';
if ~exist(filename,'file')
	unix('wget http://hci.tugraz.at/~schloegl/bci/competition2008/B0101T.gdf');
end; 
	
% Step 1a: load data ============================%
fprintf(1,'Step 1: Prepare data.\n',filename);
fprintf(1,'\ta: Load data file %s.\n',filename);
[s,HDR]=sload(filename);

% Step 1b: extract trigger and classlabels (if not already available) ==%
fprintf(1,'\tb: Extract trigger and classlabels.\n');
%--------  extraction from trigger channel: 
% HDR.TRIG = gettrigger(s(:,triggerchannel)); 
%--------- extraction from event table 
ix = find((HDR.EVENT.TYP>hex2dec('300'))&(HDR.EVENT.TYP<hex2dec('30d'))); % 0x0300..0x03ff
[i,j,HDR.Classlabel] = unique(HDR.EVENT.TYP(ix));
HDR.TRIG = HDR.EVENT.POS(HDR.EVENT.TYP==hex2dec('300'));
%HDR.TRIG = HDR.EVENT.POS(ix);
t0 = HDR.EVENT.POS(ix);
t0 = (t0(1) - HDR.TRIG(1))/HDR.SampleRate; 	% time between trial start and cue; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Preprocessing and artifact processing ====================%
%   2a: Overflowdetection: eeghist.m, [Schlögl et al. 1999] 
fprintf(1,'Step 2: Preprocessing.\n');
fprintf(1,'\ta: Quality control with histogram analysis [Schloegl et al. 1999].\n');
%Q = eeg2hist(filename); 
Q = eeg2hist(filename,'manual'); % manual selection of threshold
clear H;
H.FileName  = Q.FileName;
H.THRESHOLD = Q.THRESHOLD; 
[HDR]=sopen(H,'r',0);[s,HDR]=sread(HDR);HDR=sclose(HDR);
%   2b: Muscle detection: detectmuscle.m
fprintf(1,'\tb: Detection of muscle artifacts.\n');

%   2c: resampling 
fprintf(1,'\tc: resampling.\n');
DIV = 1; 
s = rs(s,DIV,1);   % downsampling by a factor of DIV; 
%s = rs(s,256,100); % downsampling from 256 to 100 Hz. 
HDR.SampleRate = HDR.SampleRate/DIV; 
HDR.EVENT.POS = round(HDR.EVENT.POS/DIV); 
HDR.EVENT.DUR = round(HDR.EVENT.DUR/DIV); 
HDR.TRIG      = round(HDR.TRIG); 

%   2d: Correction of EOG artifacts: regress_eog.m, get_regress_eog.m   
% 		[Schlögl et al. 2007]
fprintf(1,'\td: reduce EOG artifacts.\n');
eogchan=identify_eog_channels(filename); 
	% eogchan can be matrix in order to convert 
      	%     monopolar EOG to bipolar EOG channels
eegchan=find(HDR.CHANTYP=='E'); % exclude any non-eeg channel. 
%R = regress_eog(s,eegchan,eogchan); 
%s = s*R.r0; 	% reduce EOG artifacts 

%   2e: spatial filters 
%       spatial filters can be used to focus on specific areas.
%       examples are bipolar (BIP), common average reference (CAR), 
%       Hjorth's Laplace derivation (LAP), Common spatiol patterns (CSP)
fprintf(1,'\te: spatial filters.\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Feature extraction ====================%
eegchan=find(HDR.CHANTYP=='E'); % select EEG channels 
fprintf(1,'Step 3: feature extraction.\n');

p = 9; 
MODE.MOP = [0,p,0];	% order of AAR model
MODE.UC  = 0.0085;		% update coefficient of AAR model %    3a: Time domain parameters 

%    3a: Time-dependent parameters - motivated by log(Hjorth)
fprintf(1,'\ta: TDP (log(Hjorth)).\n');
f1 = tdp(s(:,eegchan),p,MODE.UC);

%    3b: Adaptive Autoregressive Parameters
fprintf(1,'\tb: Adaptive Autoregressive parameters (AAR).\n');
f2 = [];
for ch = 1:length(eegchan),
	X = tvaar(s(:,eegchan(ch)),MODE.MOP,MODE.UC); 
     	X = tvaar(s(:,eegchan(ch)),X);		% AAR estimation
     	f2 = [f2,X.AAR,log(X.PE)];	
end; 

%    3c: bandpower
fprintf(1,'\tc: bandpower.\n');
bands = [10,12;16,24]; % define frequency bands 
win = 1; 	% length of smoothing window in seconds
s1 = s; s1(isnan(s1))=0; 
f3 = bandpower(s, HDR.SampleRate, bands, win);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Classification ==================== %
fprintf(1,'Step 4: classification.\n');
TrialLen = 9; % seconds
SegmentLen = 100; % samples
NoS = ceil(9*HDR.SampleRate/SegmentLen);
MODE.T   = reshape((1:NoS*SegmentLen),SegmentLen,NoS)';
% valid segments for building classifier must be after cue.
MODE.WIN = MODE.T(:,1) > 3*HDR.SampleRate+1;	% cue @ t=3s.
MODE.t0= t0;
MODE.t = [min(MODE.T(:)):max(MODE.T(:))]'/HDR.SampleRate;
MODE.Segments = MODE.T;
MODE.Fs = HDR.SampleRate;  

%%%% X-V based on Jackknife-procedure (leave-K-trials-out)
K = 1; 
NG = ceil([1:length(HDR.Classlabel)]'/K);

TYPE.TYPE = 'LDA';	% classifier type 
% TYPE.hyperparameters.gamma = 0; 	% and its hyperparameters
% other possible classifiers are: MDA, LD2, LD3, LD4, RDA, GDBC, SVM, RBF, NBC, aNBC, LDA/GSVD, MDA/GSVD, LDA/sparse, MDA/sparse 

% search best segment, cross-validation using jackknife procedure, 1-vs-rest classifier 
CC1 = findclassifier(f1, HDR.TRIG, [HDR.Classlabel,NG], MODE, [], TYPE);
CC2 = findclassifier(f2, HDR.TRIG, [HDR.Classlabel,NG], MODE, [], TYPE);
CC3 = findclassifier(f3, HDR.TRIG, [HDR.Classlabel,NG], MODE, [], TYPE);

% For online feedback, the weights of the linear classifier 
%   are available through 
fprintf(1,'\tb: weights of linear classifier.\n');
CC1.weights
CC2.weights
CC3.weights
% the first element represents the bias. 

% The chosen time segment used for computing the classifiers are: 
fprintf(1,'\tc: choosen time segment.\n');
MODE.T(CC1.TI,[1,end])/HDR.SampleRate, 
MODE.T(CC2.TI,[1,end])/HDR.SampleRate, 
MODE.T(CC3.TI,[1,end])/HDR.SampleRate, 
 
% Accordingly, the time-varying distance is available 
d1 = [ones(size(f1,1),1),f1]*CC1.weights;
d2 = [ones(size(f2,1),1),f2]*CC2.weights;
d3 = [ones(size(f3,1),1),f3]*CC3.weights;
% Note, if the same a1,a2,a3 were already used for classifier training, 
%   these results are subject to overfitting. 

% Alternatively, you can setup your own cross-validation procedure 
% using train_sc.m and test_sc.m.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: Show results  ==================== %
fprintf(1,'Step 5: classifier.\n');

% the various evaluation criteria are discussed in Schlogl, Kronegg et 2007. 
fprintf(1,'\t Fig 1: results from TDP+%s.\n',TYPE.TYPE);
figure(1);
plota(CC1)

fprintf(1,'\t Fig 2: results from AAR+%s results.\n',TYPE.TYPE);
figure(2);
plota(CC2)

fprintf(1,'\t Fig 3: results from BandPower+%s results.\n',TYPE.TYPE);
figure(3);
plota(CC3)


fprintf(1,'\t Fig 3+: various evaluation criteria [Schlögl et al. 2007] for comparing different features ');

LEG= {['TDP+',TYPE.TYPE],['AAR+',TYPE.TYPE],['BP+',TYPE.TYPE]};
M = length(unique(HDR.Classlabel));
FFIELD = {'ERR','r','I','SNR','AUC','ACC00','KAP00','I_wolpaw','I_Nykopp','STMI'}; 
TIT = {'Error rate','correlation coefficient','Mutual information','Signal-to-Noise ratio','Area-under-the-ROC-curve','Accuracy','Cohens kappa coefficient','Information transfer [Wolpaw]','Information transfer [Nykopp]','Steepness of mutual information'};
PhysDim = {'[1]','[1]','[bit]','[1]','[1]','[1]','[1]','[bit]','[bit]','[bit/s]'};
for k=1:length(FFIELD),
	figure(k+3); 
	ffield = FFIELD{k};
fprintf(1,'\t Fig %i: %s (CC.TSD.%s)\n',k+3,TIT{k},ffield);
	if strcmp(ffield,'STMI') 	% steepness of mutual information 
		t = CC1.T.t-CC1.T.t0;
		t(t<0.5) = NaN; 
		plot(CC1.T.t,[sum(CC1.TSD.I,2),sum(CC2.TSD.I,2),sum(CC3.TSD.I,2)]./t(:,[1,1,1])*(M-1)/M);
	elseif (size(getfield(CC1.TSD,ffield),2)>1)	
		plot(CC1.T.t,[sum(getfield(CC1.TSD,ffield),2),sum(getfield(CC2.TSD,ffield),2),sum(getfield(CC3.TSD,ffield),2)]*(M-1)/M);
	else
		plot(CC1.T.t,[getfield(CC1.TSD,ffield),getfield(CC2.TSD,ffield),getfield(CC3.TSD,ffield)]);
	end; 	
	legend(LEG);
	ylabel([ffield,' ',PhysDim{k}]); 
	title(TIT{k})
end; 

