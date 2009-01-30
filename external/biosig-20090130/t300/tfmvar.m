function [R] = tfmvar(s,TRIG,T,MOP,f,Fs,cl)
% TFMVAR Time-Frequency MVAR analysis
%   time-frequency analysis of 
%   multivariate stochastic processes. 
%
% [R] = tfmvar(s,TRIG,T,MOP,f,Fs, [CL,group])
%
% INPUT: 
%  s    signal data (one channel per column) 
%  TRIG trigger time points (in SAMPLES)
%  T    windows definition; each column defines one window)
%       T(1,:) and T(2,:) indicate start and end [in samples], respectivly  
%  MOP  model order of the MVAR model   
%  f    designated frequencies 
%  Fs   sampling rate. 
%  [CL,group]  is OPTIONAL
%	CL 	are the labels for different classes, conditions, states. 
%		CL must be a column vector having the same length than TRIG
%       group 	is useful for controlling the resampling
%		same numbers indicate that member belongs to the same group. 
%		E.g. if data from several subjects are concatanated, and the 
%		the trials of each subject have the same numbers, the standard error 
%		of the group-statistic is estimated. 
%		If group is empty [default], each trial gets a different number; 
%		Accordingly, a trial-based leave-on-out-method (LOOM) is used, 
%		for computing the standard error. 
%               
%
%
% OUTPUT: 
%     	M and SE contain the mean 
%	and the standard error of the mean  
%       of the following characteristic parameters. 
%       The size of the parameters is defined by the number of channels,
%       the number of windows the number of designated
%       frequencies [size(s,2), size(T,1), length(f)] respectively. 
%
% univariate:
%   S1		Autospectra
%   logS1   	log(abs(S1))
%   AR1         univariate autoregressive parameters 
%   C1          variance of predication error 
%
% multivariate: 
%   S		Auto- and Cross-spectra
%   h		transfer function 		
%   logS   	log(abs(S))
%   logh   	log(abs(h))
%   y1i         imaginary part of amplitude spectra 
%   h1i         imaginary part of transfer function
%   phaseS   	phase of S
%   phaseh   	phase of h
%   COH		coherence
%   coh		coherence neglecting the cross-correlation 
%		  due to the innovation process
%   pCOH 	partial coherence
%   PDC	 	partial directed coherence [2, 5]
%   DTF 	directed transfer function [3, 6]
%   dDTF 	modified DTF [8]
%   ffDTF 	modified DTF [8]
%   AR		MVAR parameters
%   C		covariance matrix of the innovation process	
%   DC		directed granger causality [2,3,5,6]
%   GGC		Geweke's Granger Causality (not quite the same as in [12,13]
%   Af		Frequency transform of A(z)
%
% [R] = tfmvar(s,TRIG,T,MOP,f,Fs)
%   R is a struct containing M and SE as well as a few more 
%     parameters for visualization
%
%  The standard error is calculated with a jackknife-method,
%  based on LEAVE-K-TRIALs-OUT. Therefore, the SE need to be 
%  rescaled, depending on the needs [10,11]. 
%     SE 
%	standard error of the mean from the bootstrap results 
%	This has usually no common meaning (pretty much useless). 
%     SE*(N-K)^(1/2) 
%	standard deviation of the means from the bootstrapping
%	It can be interpreted as the standard error of the total mean 
%	(across all trials).
%	This value becames smaller if the number of trials increase. 	
%     SE*(N-K) 	
%	average standard error of the mean (based on a single trial).
%	This value provides a realistic value for the confidence 
%	interval of the estimates and can be used to test the 
%	significance. 
%     SE*(N-K)*N^(1/2) 
%	[estimated] standard deviation of a single trial estimate
%	This value is important for a single-trial classification.  
%		
%
% see also: tsa/MVAR, tsa/MVFREQZ, PLOTA
%
% Reference(s):
% [1] Kay S. M., Marple S. L., Spectrum Analysis - A Modern Perspective, Proc. IEEE, 1981
% [2] Baccala L. A., Sameshima K., Partial Directed Coherence: A New Concept in Neural Structure Determination, Biological Cybernetics 84, 2001
% [3] Kaminski M., Blinowska K., Szelenberger W., Topographic Analysis of Coherence and Propagation of EEG Activity During Sleep and Wakefulness, Electroencephalography and Clinical Neurophysiology 102, 1997
% [4] Franaszczuk P. J., Bergey G. K., An Autoregressive Method for the Measurement of Synchronization of Interictal and Ictal EEG Signals, Biological Cybernetics 81, 1999
% [5] Sameshima K., Baccala L. A., Using Partial Directed Coherence to Describe Neuronal Ensemble Interactions, Journal of Neuroscience Methods 94, 1999
% [6] Kaminski M., Ding M., Truccolo W. A., Bressler S. L., Evaluating Causal Relations in Neural Systems: Granger Causality, Directed Transfer Function and Statistical Assessment of Significance, Biological Cybernetics 85, 2001
% [7] Liang H., Ding M., Bressler S. L., On the Tracking of Dynamic Functional Relations in Monkey Cerebral Cortex, Neurocomputing, 2000
% [8] Korzeniewska A., Manczak M., Kaminski M., Blinowska K. J., Kasicki S., Determination of Information Flow Direction Among Brain Structures By a Modified Directed Transfer Function (dDTF) Method, Journal of Neuroscience Methods 125, 2003
% [9] A. Schl\"ogl, Comparison of Multivariate Autoregressive Estimators. Signal processing, Elsevier B.V. (in press). 
%       available at http://dx.doi.org/doi:10.1016/j.sigpro.2005.11.007
% [10] http://www.physics.utah.edu/~detar/phycs6730/handouts/jackknife/jackknife/jackknife.html
% [11] http://www-stat.stanford.edu/~susan/courses/s208/node16.html
% [12] Geweke J., 1982. J.Am.Stat.Assoc., 77, 304-313.
% [13] Bressler S.L., Richter C.G., Chen Y., Ding M. (2007)
%	Cortical fuctional network organization from autoregressive modelling of loal field potential oscillations.
%	Statistics in Medicine, doi: 10.1002/sim.2935 

%	$Revision: 1.1 $
%	$Id: tfmvar.m,v 1.1 2009-01-30 06:04:44 arno Exp $
%	Copyright (C) 2004,2005,2006,2007 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.

FLAG.SEM = 1; 

R.datatype = 'TF-MVAR 2.0';	% C = PE(MOP+1)
%R.datatype = 'TF-MVAR';	% SE = SE/N,  C = PE(MOP+1)'

R.T = T; %(T(1,:)+T(2,:))/(2*Fs); 
R.SampleRate = Fs; 
R.F = f; 
R.MOP = MOP; 
R.nan_ratio = zeros(1,size(T,2));

if nargin<7,
	cl = ones(length(TRIG),1); 
end;	
tmp=cl;tmp(isnan(tmp))=0;
if any(rem(tmp,1) & ~isnan(cl)),
        fprintf(2,'Error %s: class information is not integer\n',mfilename);
        return;
end;
if length(TRIG)~=length(cl);
        fprintf(2,'Warning TFMVAR: number of Triggers do not match class information');
end;
if size(cl,1)~=length(TRIG);    
        fprintf(2,'Warning TFMVAR: Classlabels must be a column vector');
        if length(TRIG)==size(cl,2),
                cl = cl';
        end;
end;

TRIG = TRIG(:);
TRIG = TRIG(~any(isnan(cl),2));
cl   = cl(~any(isnan(cl),2),:);
if size(cl,2)>1,
        cl2 = cl(:,2);          % 2nd column contains the group definition, ( Leave-One (Group) - Out ) 
        cl  = cl(:,1); 
else
        cl2 = [1:length(cl)]';  % each trial is a group (important for cross-validation); Trial-based LOOM  
end;
CL = unique(cl);


% univariate
[N,m] = size(s);

sz = [m,length(f),size(T,2)];
R.M.S1 = zeros(sz);
R.M.logS1 = zeros(sz);

sz(2)   = MOP;
R.M.AR1 = zeros(sz);
sz(2)   = MOP+1;
R.M.C1  = zeros(sz);

if FLAG.SEM,
        R.SE1 = R.M;
end;

sz    = [m,length(f),length(TRIG)];
S1    = zeros(sz);
sz(2) = MOP;
A1    = zeros(sz);
sz(2) = MOP+1;
C1    = zeros(sz);

% multivariate 
[N,m] = size(s);

sz = [m,m,length(f),size(T,2)];
R.M.S = zeros(sz);
R.M.h = zeros(sz);
%R.M.iY = zeros(sz);
R.M.ih = zeros(sz);
R.M.dT = zeros(sz);
R.M.phaseS = zeros(sz);
R.M.phaseh = zeros(sz);
R.M.logS   = zeros(sz);
R.M.logh   = zeros(sz);
R.M.COH    = zeros(sz);
R.M.coh    = zeros(sz);
R.M.PDC    = zeros(sz);
R.M.PDCF   = zeros(sz);
R.M.GPDC   = zeros(sz);
R.M.DTF    = zeros(sz);
R.M.pCOH   = zeros(sz);
R.M.dDTF   = zeros(sz);
R.M.ffDTF  = zeros(sz);
R.M.pCOH2  = zeros(sz);
R.M.GGC  = zeros(sz);
R.M.Af  = zeros(sz);
sz(3)  = MOP;
R.M.A  = zeros(sz);     
sz(3)  = 1;
R.M.DC = zeros(sz);
R.M.C  = zeros(sz);
sz(2)  = m*MOP;
R.M.AR = zeros(sz);

R.SE  = R.M;

r=R; clear R; 
for k0 = 1:length(CL); 
	R{k0} = r;
	ix0 = find(cl==CL(k0));	
	trig = TRIG(ix0);
       	CL2 = unique(cl2(ix0));
	R{k0}.N = length(CL2); 
     	
tic;
for k1 = 1:size(T,2),
        [S0,sz0] = trigg(s,trig,T(1,k1),T(2,k1),1+MOP);
        R{k0}.nan_ratio(k1) = mean(isnan(S0(:)));
        % univariate
        [AR1,RC1,PE1] = durlev(acovf(S0,MOP));
        
        for k=1:m;
                [h1,f] = freqz(sqrt(PE1(k,MOP+1)/(Fs*2*pi)),ar2poly(AR1(k,:)),f,Fs);
                H1(k,:)= h1(:)'; %F(:,k)=f(:);
        end;
        
        %[S(:,:,:,k2),  h(:,:,:,k2), PDC(:,:,:,k2), COH(:,:,:,k2), DTF(:,:,:,k2), DC(:,:,1,k2), pCOH(:,:,:,k2), dDTF(:,:,:,k2), ffDTF(:,:,:,k2), pCOH2(:,:,:,k2),coh(:,:,:,k2)] = mvfreqz(X.B,X.A,X.C,f,Fs);
        
        R{k0}.M.S1(:,:,k1) = H1; 
        R{k0}.M.logS1(:,:,k1) = log(abs(H1)); 
        R{k0}.M.AR1(:,:,k1) = AR1; 
        R{k0}.M.C1(:,:,k1)  = PE1;
        
        % multivariate 
        [A,RC,PE] = mvar(S0',MOP,5);
        X.A = [eye(size(S0,1)),-A];
        X.B = [eye(size(S0,1))];
        X.C = PE(:,MOP*size(S0,1)+(1:size(S0,1)));
        X.datatype = 'MVAR';
        
        [S, h, PDC, COH, DTF, DC, pCOH, dDTF, ffDTF, pCOH2, PDCF, coh, GGC, Af, GPDC] = mvfreqz(X.B,X.A,X.C,f,Fs);

	R{k0}.M.dT(:,:,:,k1)     = angle(S)./reshape(repmat(2*pi*f(:),[m*m,1]),[m,m,length(f)]);
        R{k0}.M.phaseS(:,:,:,k1) = angle(S);
        R{k0}.M.phaseh(:,:,:,k1) = angle(h);
        R{k0}.M.S(:,:,:,k1)      = S;   % power spectra
        R{k0}.M.h(:,:,:,k1)      = h;   % complex amplitude spectra (auto- and cross-spectra)
        R{k0}.M.ih(:,:,:,k1)     = imag(h);   % complex amplitude spectra (auto- and cross-spectra)
%        R{k0}.M.iY(:,:,:,k1)     = h.*repmat(sqrtm(X.C/(2*pi*Fs)),[1,1,length(f)]);   % complex amplitude spectra (auto- and cross-spectra)
        R{k0}.M.logS(:,:,:,k1)   = log(abs(S));
        R{k0}.M.logh(:,:,:,k1)   = log(abs(h));
        R{k0}.M.PDC(:,:,:,k1)    = PDC;
        R{k0}.M.PDCF(:,:,:,k1)   = PDCF;
        R{k0}.M.GPDC(:,:,:,k1)   = GPDC;
        R{k0}.M.COH(:,:,:,k1)    = abs(COH);
        R{k0}.M.coh(:,:,:,k1)    = abs(coh);
        R{k0}.M.iCOH(:,:,:,k1)   = imag(COH);
        R{k0}.M.icoh(:,:,:,k1)   = imag(coh);
        R{k0}.M.DTF(:,:,:,k1)    = DTF;
        R{k0}.M.pCOH(:,:,:,k1)   = pCOH;
        R{k0}.M.dDTF(:,:,:,k1)   = dDTF;
        R{k0}.M.ffDTF(:,:,:,k1)  = ffDTF;
        R{k0}.M.pCOH2(:,:,:,k1)  = pCOH2;
        R{k0}.M.GGC(:,:,:,k1)    = GGC;
        R{k0}.M.Af(:,:,:,k1)     = Af;
        R{k0}.M.AR(:,:,1,k1)     = A;
        R{k0}.M.DC(:,:,1,k1)     = DC;
        R{k0}.M.C(:,:,1,k1)      = X.C;
        
        if (length(CL2)>1),  % jackknife resampling. 
        
		sz    = [m,m,length(f),length(CL2)];
		S     = zeros(sz);
		h     = zeros(sz);
		COH   = zeros(sz);
		coh   = zeros(sz);
		PDC   = zeros(sz);
		PDCF  = zeros(sz);
		DTF   = zeros(sz);
		pCOH  = zeros(sz);
		dDTF  = zeros(sz);
		ffDTF = zeros(sz);
		pCOH2 = zeros(sz);		
		sz(3) = 1;
		C     = zeros(sz);
		DC    = zeros(sz);
		sz(2) = m*MOP;
		AR    = zeros(sz);


                for k2 = 1:length(CL2),
                        sel = (cl2(ix0)~=CL2(k2)); 
                        [S0,sz0] = trigg(s, trig(sel), T(1,k1), T(2,k1), 1+MOP);
                        %fprintf(2,'\nExtract epochs done\nNumber of Trials: %i : %i\n',length(t0),length(t1));
                        
                        % univariate
                        [AR1,RC1,PE1] = durlev(acovf(S0,MOP));
                        
                        for k=1:m;
                                [h1,f] = freqz(sqrt(PE1(k,MOP+1)/(Fs*2*pi)),ar2poly(AR1(k,:)),f,Fs);
                                H1(k,:)= h1(:)'; %F(:,k)=f(:);
                        end;
                        
                        %[S(:,:,:,k2),  h(:,:,:,k2), PDC(:,:,:,k2), COH(:,:,:,k2), DTF(:,:,:,k2), DC(:,:,1,k2), pCOH(:,:,:,k2), dDTF(:,:,:,k2), ffDTF(:,:,:,k2), pCOH2(:,:,:,k2),coh(:,:,:,k2)] = mvfreqz(X.B,X.A,X.C,f,Fs);
                        
                        S1(:,:,k2) = H1; 
                        A1(:,:,k2) = AR1; 
                        C1(:,:,k2) = PE1;
                        
                        % multivariate 
                        [A,RC,PE] = mvar(S0',MOP,5);
                        X.A = [eye(size(S0,1)),-A];
                        X.B = [eye(size(S0,1))];
                        X.C = PE(:,MOP*size(S0,1)+(1:size(S0,1)));
                        X.datatype = 'MVAR';
                        
                        [S(:,:,:,k2),  h(:,:,:,k2), PDC(:,:,:,k2), COH(:,:,:,k2), DTF(:,:,:,k2), DC(:,:,1,k2), pCOH(:,:,:,k2), dDTF(:,:,:,k2), ffDTF(:,:,:,k2), pCOH2(:,:,:,k2),PDCF(:,:,:,k2), coh(:,:,:,k2), GGC(:,:,:,k2), Af(:,:,:,k2), GPDC(:,:,:,k2)] = mvfreqz(X.B,X.A,X.C,f,Fs);

%                        Y(:,:,:,k1)  = h(:,:,:,k2).*repmat(sqrtm(X.C/(2*pi*Fs)),[1,1,length(f)]);   % complex amplitude spectra (auto- and cross-spectra)

                        AR(:,:,1,k2) = A;
                        C(:,:,1,k2)  = X.C;
                        
                        fprintf(2,'%.0f\t%.0f\t%.1f\n',[k1,k2,toc]);tic;
                end;
                
                % univariate
                [R{k0}.SE.S1(:,:,k1),   R{k0}.M.S1(:,:,k1)  ] = sem(abs(S1),3);
                [R{k0}.SE.logS1(:,:,k1),R{k0}.M.logS1(:,:,k1)] = sem(log(abs(S1)),3);
                [R{k0}.SE.AR1(:,:,k1),  R{k0}.M.AR1(:,:,k1) ] = sem(A1,3);
                [R{k0}.SE.C1(:,:,k1),   R{k0}.M.C1(:,:,k1)  ] = sem(C1,3);

                % multivariate 
                [R{k0}.SE.dT(:,:,:,k1),     R{k0}.M.dT(:,:,:,k1)]     = sem(angle(S)./permute(reshape(repmat(2*pi*f(:)',[m*m*length(CL2),1]),[m,m,length(CL2),length(f)]),[1,2,4,3]),4);
                [R{k0}.SE.phaseS(:,:,:,k1), R{k0}.M.phaseS(:,:,:,k1)] = sem(angle(S),4);
                [R{k0}.SE.phaseh(:,:,:,k1), R{k0}.M.phaseh(:,:,:,k1)] = sem(angle(h),4);
%                [R{k0}.SE.iY(:,:,:,k1),     R{k0}.M.iY(:,:,:,k1)]     = sem(imag(Y),4);
                [R{k0}.SE.ih(:,:,:,k1),     R{k0}.M.ih(:,:,:,k1)]     = sem(imag(h),4);
                [R{k0}.SE.dT1(:,:,:,k1),    R{k0}.M.dT1(:,:,:,k1)]    = sem(angle(coh)./permute(reshape(repmat(2*pi*f(:)',[m*m*length(CL2),1]),[m,m,length(CL2),length(f)]),[1,2,4,3]),4);
                [R{k0}.SE.phaseS1(:,:,:,k1),R{k0}.M.phaseS1(:,:,:,k1)]= sem(angle(coh),4);
                [R{k0}.SE.phaseS2(:,:,:,k1),R{k0}.M.phaseS2(:,:,:,k1)]= sem(angle(COH),4);  % for testing: should be same than phaseS
                [R{k0}.SE.S(:,:,:,k1),      R{k0}.M.S(:,:,:,k1)     ] = sem(S,4);
                [R{k0}.SE.iSpectrum(:,:,:,k1),      R{k0}.M.iSpectrum(:,:,:,k1)     ] = sem(imag(S),4);
                [R{k0}.SE.h(:,:,:,k1),      R{k0}.M.h(:,:,:,k1)     ] = sem(h,4);
                [R{k0}.SE.logS(:,:,:,k1),   R{k0}.M.logS(:,:,:,k1)  ] = sem(log(abs(S)),4);
                [R{k0}.SE.logh(:,:,:,k1),   R{k0}.M.logh(:,:,:,k1)  ] = sem(log(abs(h)),4);
                [R{k0}.SE.PDC(:,:,:,k1),    R{k0}.M.PDC(:,:,:,k1)   ] = sem(PDC,4);
                [R{k0}.SE.PDCF(:,:,:,k1),   R{k0}.M.PDCF(:,:,:,k1)  ] = sem(PDCF,4);
                [R{k0}.SE.GPDC(:,:,:,k1),   R{k0}.M.GPDC(:,:,:,k1)  ] = sem(GPDC,4);
                [R{k0}.SE.COH(:,:,:,k1),    R{k0}.M.COH(:,:,:,k1)   ] = sem(abs(COH),4);
                [R{k0}.SE.coh(:,:,:,k1),    R{k0}.M.coh(:,:,:,k1)   ] = sem(abs(coh),4);
                [R{k0}.SE.iCOH(:,:,:,k1),   R{k0}.M.iCOH(:,:,:,k1)  ] = sem(imag(COH),4);
                [R{k0}.SE.icoh(:,:,:,k1),   R{k0}.M.icoh(:,:,:,k1)  ] = sem(imag(coh),4);
                [R{k0}.SE.DTF(:,:,:,k1),    R{k0}.M.DTF(:,:,:,k1)   ] = sem(DTF,4);
                [R{k0}.SE.pCOH(:,:,:,k1),   R{k0}.M.pCOH(:,:,:,k1)  ] = sem(pCOH,4);
                [R{k0}.SE.dDTF(:,:,:,k1),   R{k0}.M.dDTF(:,:,:,k1)  ] = sem(dDTF,4);
                [R{k0}.SE.ffDTF(:,:,:,k1),  R{k0}.M.ffDTF(:,:,:,k1) ] = sem(ffDTF,4);
                [R{k0}.SE.pCOH2(:,:,:,k1),  R{k0}.M.pCOH2(:,:,:,k1) ] = sem(pCOH2,4);
                [R{k0}.SE.GGC(:,:,:,k1),    R{k0}.M.GGC(:,:,:,k1)   ] = sem(GGC,4);
                [R{k0}.SE.Af(:,:,:,k1),     R{k0}.M.Af(:,:,:,k1)    ] = sem(abs(Af),4);
                [R{k0}.SE.AR(:,:,1,k1),     R{k0}.M.AR(:,:,1,k1)    ] = sem(AR,4);
                [R{k0}.SE.DC(:,:,1,k1),     R{k0}.M.DC(:,:,1,k1)    ] = sem(DC,4);
                [R{k0}.SE.C(:,:,1,k1),      R{k0}.M.C(:,:,1,k1)     ] = sem(C,4);
                
        end; 
        %% In order to obtain the standard error of the mean, 
        %% this SE must be multiplied by N-1, with N = length(TRIG)
        
end; 	%for k1


%%%%% scaling of standard error 
GF = fieldnames(R{k0}.SE);
for k2=1:length(GF),
      	R{k0}.SE = setfield(R{k0}.SE,GF{k2},getfield(R{k0}.SE,GF{k2})*sqrt((R{k0}.N-1)*(length(CL2)-1)));
end;

if any(R{k0}.nan_ratio),
        fprintf(2,'Warning TFMVAR: up to %5.2f%% are missing values (NaN).\n',max(R{k0}.nan_ratio)*100);
        fprintf(2,'This may cause underestimating the standard error.\n');
end;	


end; 	%for k0



