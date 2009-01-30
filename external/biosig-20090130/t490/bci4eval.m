function [o] = bci4eval(tsd,TRIG,cl,pre,post,Fs)
% BCI4eval evaluates a BCI-result for two and more classes
%
%   Two classes are evaluated like in [1,2]:
%   - It returns the classification error, the signal to noise ratio, 
%   the mutual information, as well as mean, standard error, 
%   within-class accuracy and standard deviation for both classes. 
%   - time course of these resulting parameters are supported
%
%   More than two classes are evaluated with 
%   - Kappa coefficient including standard deviation 
%   - Accuracy
%   
%   Missing values can be encoded as NaN.
%
% X = bci4eval(tsd,trig,cl,pre,post,Fs)
% INPUT:
%       tsd     continous output 
%               for 2 classes, tsd must have size Nx1 
%               size NxM for M-classes, for each row the largest value 
%               determines the assigned class 
%       trig    trigger time points
%       cl      classlabels
%       pre     offset of trial start 
%       post    offset of trial end 
%       Fs      sampling rate;
%
% OUTPUT: 
%       X is a struct with various results  
%       2-classes:
%               X.MEAN1, XMEAN2: mean of both classes      
%               X.ERR           error rate 
%               X.p_value       significance level of paired t-test 
%               X.SNR           signal-to-noise ratio
%               X.I             mutual information
%               X.AUC           area-under-the-(ROC) curve
%       N(>2)-classes:
%               X.KAP00         Cohen's kappa coefficient
%               X.Ksd00         standard error of kappa coefficient 
%               X.ACC00         accuracy 
%
%               X.MEAN0         average output of non-active class                
%               X.MEAN1         average output of active class
%               X.SNR           signal-to-noise ratio for each class
%               X.I             mutual information for each class
%               X.AUC           area-under-the-(ROC) curve for each class
%               X.r             correlation coefficient (parametric) 
%               X.rankcorrelation    rank correlation (non-parametric)
%               X.I_Nykopp      Nykopp's mutual information
%               X.I_Wolpaw      Wolpaws mutual information 
%
%
% see also: SUMSKIPNAN, PLOTA, BCI3EVAL
%
% REFERENCES: 
%  [1] Schlögl A., Neuper C. Pfurtscheller G.
%	Estimating the mutual information of an EEG-based Brain-Computer-Interface
%	Biomedizinische Technik 47(1-2): 3-8, 2002.
%  [2] A. Schlögl, C. Keinrath, R. Scherer, G. Pfurtscheller,
%	Information transfer of an EEG-based Bran-computer interface.
%	Proceedings of the 1st International IEEE EMBS Conference on Neural Engineering, pp.641-644, Mar 20-22, 2003. 
%  [3]  A. Schlögl, Evaluation of the dataset III of the BCI-competition 2003. 
%	http://ida.first.fraunhofer.de/projects/bci/competition/results/TR_BCI2003_III.pdf
% [4] Schlögl A, Kronegg J, Huggins JE, Mason SG;
%	Evaluation criteria in BCI research.
%	(Eds.) G. Dornhege, J.R. Millan, T. Hinterberger, D.J. McFarland, K.-R.Müller;
%	Towards Brain-Computer Interfacing, MIT Press, p327-342, 2007


%    $Id: bci4eval.m,v 1.1 2009-01-30 06:04:50 arno Exp $
%    Copyright (C) 2003,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    This is part of the BIOSIG-toolbox http://biosig.sf.net/

% BioSig is free software; you can redistribute it and/or
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


if nargin<6
        Fs = 1;
end;

DIM = 2; 
CL = unique(cl(~isnan(cl)));M = length(CL); 
M1 = size(tsd,2);
%CL= 0:M; M,M1,CL,

if any([1,M]==size(tsd,2))
        [x,sz] = trigg(tsd,TRIG,pre,post);
        D = reshape(x,sz);
        D = squeeze(D);
else
        if size(tsd,1)==length(cl)
                D = tsd';
        elseif size(tsd,2)==length(cl)
                D = tsd;
        else
                error('BCI4EVAL: size of data and size of Classlabels does not fit');
        end;
        sz = [1,size(D)];
        pre=1;
        post=size(D,1);
end;
if size(D,(M==size(tsd,2))+2) ~= length(cl),
        size(D,(M==size(tsd,2))+2),
        [size(D), size(cl),size(CL),size(tsd)],
        error('BCI4EVAL: length of Trigger and Length of Classlabels must fit')
end;

% Time axis
o.T = [pre:post]'/Fs;
if (M==2) & (sz(1)==1),
        for k = 1:length(CL),
                X{k} = squeeze(D(:,cl==CL(k),:));
        end;

        % classification error/accuracy 
        o.ERR = (1-mean(sign([-X{1},X{2}]),DIM))/2;
        o.ACC00 = 1-o.ERR; 
        
        % within-class accuracy
        o.BCG1 = (1 + mean(sign(-X{1}), DIM))/2;
        o.BCG2 = (1 + mean(sign( X{2}), DIM))/2;
        
        %%%%% 2nd order statistics
        [i1.SUM, o.N1, i1.SSQ] = sumskipnan(X{1},DIM);       
        [i2.SUM, o.N2, i2.SSQ] = sumskipnan(X{2},DIM);       
        
        o.MEAN1 = i1.SUM./o.N1;	% mean
        v1      = i1.SSQ-i1.SUM.*o.MEAN1;	% n*var
        o.SD1   = sqrt(v1./o.N1); % standard deviation 
        %o.SE1 = sqrt(v1)./o.N1; % standard error of the mean 
        
        o.MEAN2 = i2.SUM./o.N2;
        v2    = i2.SSQ-i2.SUM.*o.MEAN2;
        o.SD2 = sqrt(v2./o.N2);
        %o.SE2 = sqrt(v2)./o.N2;
    
        %%%%% Signal-to-Noise Ratio 

        % intra-class variability
        vd = var([-X{1},X{2}],[],DIM);        

        % paired t-test
        [se,m]=sem([-X{1},X{2}],DIM);
        o.p_value = 2 * tcdf(-abs(m./se),o.N1+o.N2-2);
        
        o.SNR = 1/4*(o.MEAN2-o.MEAN1).^2./vd; 
        
        %%%%% Mutual Information 
        o.I   = 1/2*log2(o.SNR+1);
        
        %%%%% ROC, AUC
        if DIM==2, D=D'; end; 
        [s,ix]=sort(D,1);
        TNR = cumsum(cl(ix)==CL(1),1);
        FNR = cumsum(cl(ix)==CL(2),1);
        TNR = TNR./repmat(TNR(end,:),size(ix,1),1);
        FNR = FNR./repmat(FNR(end,:),size(ix,1),1);
        AUC = sum(diff(FNR,[],1) .* (TNR(1:end-1,:)+TNR(2:end,:)))/2;
        o.AUC = AUC'; 

        o.r = corrcoef(D,double(cl(:)));
        o.N = double(~isnan(D)')*double(~isnan(cl(:)));

        % o.rankcorrelation = corrcoef(D,double(cl(:)),'rank');         %
        % is SLOW 
        
        for k=1:size(D,2),
                [kap,sd,H,z,OA,SA] = kappa(cl(:),D(:,k)>0);
                tpr = H(1,1)/sum(H(1,:));
                fpr = H(2,2)/sum(H(2,:));
                if tpr>fpr,
                        Aprime(1,k) = 1/2 +(+tpr-fpr)*(1+tpr-fpr)/(4*tpr*(1-fpr));
                else
                        Aprime(1,k) = 1/2 -(-tpr+fpr)*(1-tpr+fpr)/(4*fpr*(1-tpr));
                end
                dprime(1,k) = norminv(tpr)-norminv(fpr); 
                %ACC1,k) = sum(diag(H))/sum(H(:)); %see above
        end;
        o.dprime = dprime';
        o.Aprime = Aprime';
        
        %%%%% OUTPUT
        o.datatype = 'TSD_BCI7';  % useful for PLOTA
        
end

if (M==sz(1)),
        for k=1:sz(1),
                x = bci4eval(tsd(:,k),TRIG,cl==CL(k),pre,post,Fs);
                if k==1; 
                        o.ERR      = x.ERR;
                        o.MEAN1    = x.MEAN1;
                        o.MEAN2    = x.MEAN2;
                        o.SD1      = x.SD1;
                        o.SD2      = x.SD2;
                        o.SNR      = x.SNR;
                        o.I        = x.I;
                        o.AUC      = x.AUC;
                        o.Aprime   = x.Aprime;
                        o.dprime   = x.dprime;
                        o.r        = x.r; 
                        o.N        = x.N; 
                        %o.rankcorrelation = x.rankcorrelation; % is slow
                        
                else
                        o.ERR(:,k)      = x.ERR;
                        o.MEAN1(:,k)    = x.MEAN1;
                        o.MEAN2(:,k)    = x.MEAN2;
                        o.SD1(:,k)      = x.SD1;
                        o.SD2(:,k)      = x.SD2;
                        o.SNR(:,k)      = x.SNR;
                        o.I(:,k)        = x.I;
                        o.AUC(:,k)      = x.AUC;
                        o.Aprime(:,k)   = x.Aprime;
                        o.dprime(:,k)   = x.dprime;
                        o.r(:,k)        = x.r; 
                        %o.rankcorrelation(:,k) = x.rankcorrelation; 
                end;
        end;
        o.CL = CL'; 

elseif sz(1)==1,
	% two-class problem with single column tsd 
        [x,sz] = trigg([-tsd,tsd],TRIG,pre,post);
        D = reshape(x,sz);
        D = squeeze(D);
else
        error('number of classes and number of traces do not fit');
end;

        
        [m,IX] = max(D,[],1);
        IX(isnan(m)) = NaN;
        IX = squeeze(IX);
        CMX = repmat(zeros,[size(IX,1),length(CL)*[1,1]]);
        for k = 1:length(CL),
        for j = 1:length(CL),
                CMX(:,k,j)=sum(IX(:,CL(k)==cl)==j,2);
        end;
        end;
        o.CMX = CMX; 
        
        o.KAP00 = zeros(size(CMX,1),1);
        o.Ksd00 = zeros(size(CMX,1),1);
        o.ACC00 = zeros(size(CMX,1),1);
        o.R     = zeros(size(CMX,1),1);

        for k   = 1:size(CMX,1),
                [o.KAP00(k),o.Ksd00(k),h,z,o.ACC00(k),sa,o.I_Nykopp(k,1)] = kappa(squeeze(CMX(k,:,:)));            
        end;
        o.datatype = 'TSD_BCI9';  % useful for PLOTA
        [tmp,o.tix]= max([o.KAP00,o.R,sum(o.I,2),wolpaw_entropy(o.ACC00,M)]); 
        o.optCMX   = squeeze(CMX(o.tix(1),:,:));%,length(CL)*[1,1]);
        o.I_wolpaw = wolpaw_entropy(o.ACC00,M);
        
