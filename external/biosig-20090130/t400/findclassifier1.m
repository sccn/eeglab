function [CC,Q,tsd,md]=findclassifier1(D,TRIG,cl,T,t0,SWITCH)
% FINDCLASSIFIER1
% 
% [CC,Q,TSD,MD]=findclassifier1(D,TRIG,Class,class_times,t_ref);
%
% D 	data, each row is one time point
% TRIG	trigger time points
% Class class information
% class_times	classification times, combinations of times must be in one row 
% t_ref	reference time for Class 0 (optional)
%
% CC 	contains LDA and MD classifiers
% Q  	is a list of classification quality for each time of 'class_times'
% TSD 	returns the LDA classification 
% MD	returns the MD  classification 
%
% [CC,Q,TSD,MD]=findclassifier(AR,find(trig>0.5)-257,~mod(1:80,2),reshape(1:14*128,16,14*8)');
%
%
% Reference(s): 
% [1] Schlögl A., Neuper C. Pfurtscheller G.
% 	Estimating the mutual information of an EEG-based Brain-Computer-Interface
%  	Biomedizinische Technik 47(1-2): 3-8, 2002.
% [2] A. Schlögl, C. Keinrath, R. Scherer, G. Pfurtscheller,
%	Information transfer of an EEG-based Bran-computer interface.
%	Proceedings of the 1st International IEEE EMBS Conference on Neural Engineering, Capri, Italy, Mar 20-22, 2003 


%   Copyright (C) 1999-2004 by Alois Schloegl <a.schloegl@ieee.org>	
%	$Id: findclassifier1.m,v 1.1 2009-01-30 06:04:48 arno Exp $


% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
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


warning('this function is obsolete and replaced by FINDCLASSIFIER');


tsd=[];md=[];

if nargin<6,
        SWITCH=0;
end;
if nargin>4,
        if isempty(t0),
                t0=logical(ones(size(T,1),1));
        end;
end;
tmp=cl;tmp(isnan(tmp))=0;
if any(rem(tmp,1) & ~isnan(cl)),
        fprintf(2,'Error %s: class information is not integer\n',mfilename);
        return;
end;
if length(TRIG)~=length(cl);
        fprintf(2,'number of Triggers do not match class information');
end;

CL = unique(cl(~isnan(cl)));
CL = sort(CL);
TRIG = TRIG(:);
if ~all(D(:,1)==1)
        %        D1=[ones(size(D,1)-1,1),diff(D)];
        D =[ones(size(D,1),1),D];
        %else
        %	D1=[ones(size(D,1)-1,1),diff(D(:,2:end))];
end;

% add sufficient NaNs at the beginning and the end 
tmp = min(TRIG)+min(min(T))-1;
if tmp<0,
        TRIG = TRIG - tmp;
        D = [repmat(nan,[-tmp,size(D,2)]);D];
end;        
tmp = max(TRIG)+max(max(T))-size(D,1);
if tmp>0,
        D = [D;repmat(nan,[tmp,size(D,2)])];
end;        

% estimate classification result for all time segments in T - without crossvalidation 
CMX = zeros([size(T,1),length(CL)*[1,1]]);
for k = 1:size(T,1),
        cmx = zeros(length(CL));
        for l = 1:length(CL), 
                t = perm(TRIG(cl==CL(l)),T(k,:));
                %t = t(t<=size(D,1));
                [C{k,l},tmp] = covm(D(t(:),:),'M');
        end;
        %[Q(k),d{k}] = qcmahal({C0r,C{k,:}});
        [CC.QC(k),d{k}] = qcmahal({C{k,:}});
        lnQ(k) = mean(log(d{k}(~eye(length(d{k})))));
        for l = 1:length(CL),
                t = perm(TRIG(cl==CL(l)),T(k,:));
                %t = t(t<=size(D,1));
                [tmp] = mdbc({C{k,:}},D(t(:),:));
                [tmp,ix] = min(tmp,[],2);
                tmp = isnan(tmp);
                ix(tmp) = NaN; %NC(1)+1;	% not classified; any value but not 1:length(MD)
                ix(~tmp) = CL(ix(~tmp));
                tmp = histo3([ix;CL(:)]);
                cmx(tmp.X,l) = tmp.H-1;            
        end;
        CMX(k,:,:)  = cmx;
        CC.KAPPA(k) = kappa(cmx);
        CC.ACC(k)   = sum(diag(cmx))/sum(cmx(:));
end;	
% identify best classification time 
if nargin>4,
        tmp = CC.QC;
        tmp(~t0) = 0;
        [maxQ,CC.TI] = max(tmp); %d{K},
else
        [maxQ,CC.TI] = max(CC.QC); %d{K},
end;
%CC.TI = K;

% build MD classifier 
CC.MD  = {C{CC.TI,:}};
CC.IR  = mdbc({C{CC.TI,:}});
CC.D   = d{CC.TI};
CC.Q   = CC.QC(CC.TI);
CC.CMX = squeeze(CMX(CC.TI,:,:));

m1=decovm(CC.MD{1});
m2=decovm(CC.MD{2});
tmp=mdbc(CC.MD,[1,m1;1,m2]);
CC.scale=[1,1]*max(abs(tmp(:)));  % element 1 

[maxQ,CC.lnTI] = max(lnQ); %d{K},
CC.DistMXln = d{CC.lnTI};
CC.MDln = {C{CC.lnTI,:}};

% alternative classifier using two different time segments.
if 1,
        [Q,d] = qcmahal({C{:}}');
        CC.T2.D = d;
        [ix,iy] = find(d==max(d(:)));
        ix=mod(ix(1)-1,size(C,1))+1;
        iy=mod(iy(1)-1,size(C,1))+1;
        CC.T2.TI = [ix(1),iy(1)];
        CC.T2.MD = {C{ix(1),1},C{iy(1),2}};
end;

% build LDA classifier
if length(CL)==2,
        % LDA 
        C0 = zeros(size(C{CC.TI,1}));
        for l=1:length(CL);
                [M{l},sd,COV,xc,N,R2] = decovm(C{CC.TI,l});
                C0 = C0 + C{CC.TI,l};
        end;
        [M0,sd,COV0,xc,N,R2] = decovm(C0);
        w     = COV0\(M{1}'-M{2}');
        w0    = M0*w;
        %CC.LDA.b = w0;
        %CC.LDA.w = -w;
        CC.lda = [w0; -w];
        
        % MD 
        tsd = ldbc(CC.MD,D);
end;
md = mdbc(CC.MD,D);
lld= llbc(CC.MD,D);

% bias correction not used anymore.
CC.BIAS.LDA=0;%mean(tsd);
CC.BIAS.MDA=0;%mean(diff(-md,[],2));
CC.BIAS.GRB=0;%mean(diff(-exp(-md/2),[],2));

%[tmp,IX] = min(md,[],2);

%% cross-validation with jackknife (leave one trial out)
nc  = max(max(T))-min(min(T))+1;

JKD   = repmat(nan,[nc,length(CL),length(TRIG)]);
JKD1  = repmat(nan,[nc,length(TRIG)]);
JKD2  = repmat(nan,[nc,length(TRIG)]);
JKLD  = repmat(nan,[nc,length(TRIG)]);
for l = find(~isnan(cl(:)'));1:length(cl);
        c = find(cl(l)==CL);
        t = TRIG(l)+T(CC.TI,:);
        %t = t(t<=size(D,1));
        [tmp,tmpn] = covm(D(t(:),:),'M');
        
        cc 	= CC.MD;
        cc{c}   = CC.MD{c}-tmp;
        
        %t = TRIG(l)+(1:nc);
        %t = t(t<=size(D,1));
        t = TRIG(l)+(min(min(T)):max(max(T)));
        
        [d,ix] = llbc(cc,D(t,:));
        if length(CL)==2,
                JKD3(:,l)=d(:,1);
                JKD4(:,l)=d(:,2);
        end;
        
        d = mdbc(cc,D(t,:));
        JKD(:,:,l) = d;
        [tmp,MDIX(:,l)] = min(d,[],2);
        
        if length(CL)==2,
                JKD1(:,l) = d(:,1);
                JKD2(:,l) = d(:,2);
                
                LDA(:,l) = ldbc(cc);
                JKLD(:,l) = D(t,:)*LDA(:,l);
        end;	
end;
[CC.ldaC0,NN] = covm(LDA','D0');
%CC.ldaC0=CC.ldaC0./NN*min(0,sum(~isnan(CL))-1); 
% since NN==min(0,sum(~isnan(CL))-1), no need to rescale

% Concordance matrix with cross-validation 
CC.mmx= zeros([size(MDIX,1),length(CL)^2]);
CC.I0 = zeros([size(MDIX,1),length(CL)]);
CC.I  = zeros([size(MDIX,1),1]);
tmp   = zeros([size(MDIX,1),length(CL)]);
for k = 1:length(CL),
        jkd = squeeze(JKD(:,k,:));
        o = bci3eval(jkd(:,cl~=k),jkd(:,cl==k),2);
        
        CC.TSD{k}  = o;
        CC.I0(:,k) = log2(2*var(jkd,[],2)./(var(jkd(:,cl==k),[],2) + var(jkd(:,cl~=k),[],2)))/2;
        
        [sum0,n0,ssq0]=sumskipnan(jkd(:,cl==k),2);
        [sum1,n1,ssq1]=sumskipnan(jkd(:,cl~=k),2);
        s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
        s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
        s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
        SNR = 2*s./(s0+s1); % this is SNR+1 
        CC.I1(:,k) = log2(SNR)/2;
        
        for l = 1:length(CL),
                tmp(:,l) = sum(MDIX(:,cl==CL(k))==CL(l),2);    
                if CL(k) == CL(l),
                        acc = tmp(:,l);
                end;
        end;
        CC.mmx(:,(1-length(CL):0)+k*length(CL)) = tmp;
        CC.acc(:,k) = acc./sum(tmp,2);
end;
CC.CMX00 = reshape(sum(CC.mmx(T(CC.TI,:),:),1),[1,1]*length(CL))/size(T,2);
CC.I = sum(CC.I0,2);
CC.ACC00 = sum(CC.mmx(:,1:length(CL)+1:end),2)/sum(~isnan(cl));	
CC.KAP00 = zeros(size(MDIX,1),1);
for k = 1:size(MDIX,1),
        CC.KAP00(k) = kappa(reshape(CC.mmx(k,:),[1,1]*length(CL)));
end;

if length(CL) > 2, 
        return; 
end; 


if bitand(SWITCH,1),
        CC.LDA.ERR00 = (mean(sign(JKLD),2)+1)/2;
        CC.MDA.ERR00 = (mean(sign(JKD1-JKD2),2)+1)/2;
        CC.GRB.ERR00 = (mean(sign(exp(-JKD2/2)-exp(-JKD1/2)),2)+1)/2;
end;


d=JKLD;
tmp1 = d(1-min(T(:))+T(CC.TI,:),cl==CL(1));
[sum0,n0,ssq0] = sumskipnan(tmp1(:));       
tmp2 = d(1-min(T(:))+T(CC.TI,:),cl==CL(2));
[sum1,n1,ssq1] = sumskipnan(tmp2(:));       
CC.LDA.AUC      = auc(tmp1,tmp2);
CC.LDA.ERR(1,1) = mean(sign([tmp1(:)]))/2+1/2;
CC.LDA.ERR(1,2) = mean(sign([tmp2(:)]))/2+1/2;
CC.LDA.ERR(2,1) = mean(sign([mean(tmp1,1)']))/2+1/2;
CC.LDA.ERR(2,2) = mean(sign([mean(tmp2,1)']))/2+1/2;
s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
SNR = 2*s./(s0+s1); % this is SNR+1 
CC.LDA.I   = log2(SNR)/2;
CC.LDA.SNR = SNR - 1;
if 0,
        clear tmp1 tmp2; 
        tmp1 = stat2(d(:,cl==CL(1)),2);       
        tmp2 = stat2(d(:,cl==CL(2)),2);       
        CC.LDA.TSD=stat2res(tmp1,tmp2);
        CC.LDA.TSD.ERR=1/2-mean(sign([-d(:,cl==CL(1)),d(:,cl==CL(2))]),2)/2;
elseif bitand(SWITCH,1),
        CC.LDA.TSD=bci3eval(d(:,cl==CL(1)),d(:,cl==CL(2)),2);
end;

d = JKD1 - JKD2;
tmp1 = d(1-min(T(:))+T(CC.TI,:),cl==CL(1));
[sum0,n0,ssq0] = sumskipnan(tmp1(:));       
tmp2 = d(1-min(T(:))+T(CC.TI,:),cl==CL(2));
[sum1,n1,ssq1] = sumskipnan(tmp2(:));       
CC.MDA.AUC      = auc(tmp1,tmp2);
CC.MDA.ERR(1,1) = mean(sign([tmp1(:)]))/2+1/2;
CC.MDA.ERR(1,2) = mean(sign([tmp2(:)]))/2+1/2;
CC.MDA.ERR(2,1) = mean(sign([mean(tmp1,1)']))/2+1/2;
CC.MDA.ERR(2,2) = mean(sign([mean(tmp2,1)']))/2+1/2;
s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
SNR = 2*s./(s0+s1); % this is SNR+1 
CC.MDA.I   = log2(SNR)/2;
CC.MDA.SNR = SNR - 1;
if 0,
        clear tmp1 tmp2;
        tmp1 = stat2(d(:,cl==CL(1)),2);       
        tmp2 = stat2(d(:,cl==CL(2)),2);       
        CC.MDA.TSD=stat2res(tmp1,tmp2);
        CC.MDA.TSD.ERR=1/2-mean(sign([-d(:,cl==CL(1)),d(:,cl==CL(2))]),2)/2;
elseif bitand(SWITCH,1),
        CC.MDA.TSD=bci3eval(d(:,cl==CL(1)),d(:,cl==CL(2)),2);
end;

d = sqrt(JKD1) - sqrt(JKD2);
tmp1 = d(1-min(T(:))+T(CC.TI,:),cl==CL(1));
[sum0,n0,ssq0] = sumskipnan(tmp1(:));       
tmp2 = d(1-min(T(:))+T(CC.TI,:),cl==CL(2));
[sum1,n1,ssq1] = sumskipnan(tmp2(:));       
CC.MD2.AUC      = auc(tmp1,tmp2);
CC.MD2.ERR(1,1) = mean(sign([tmp1(:)]))/2+1/2;
CC.MD2.ERR(1,2) = mean(sign([tmp2(:)]))/2+1/2;
CC.MD2.ERR(2,1) = mean(sign([mean(tmp1,1)']))/2+1/2;
CC.MD2.ERR(2,2) = mean(sign([mean(tmp2,1)']))/2+1/2;
s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
SNR = 2*s./(s0+s1); % this is SNR+1 
CC.MD2.I   = log2(SNR)/2;
CC.MD2.SNR = SNR - 1;
if 0,
        clear tmp1 tmp2;
        tmp1 = stat2(d(:,cl==CL(1)),2);       
        tmp2 = stat2(d(:,cl==CL(2)),2);       
        CC.MD2.TSD=stat2res(tmp1,tmp2);
        CC.MD2.TSD.ERR=1/2-mean(sign([-d(:,cl==CL(1)),d(:,cl==CL(2))]),2)/2;
elseif bitand(SWITCH,1),
        CC.MD2.TSD=bci3eval(d(:,cl==CL(1)),d(:,cl==CL(2)),2);
end;

%%%
if any(isnan(cl)),
        t = perm(TRIG(isnan(cl)),T(CC.TI,:));
        t = t(t<=size(D,1));
        tmp= rs(D(t(:),:),size(T,2),1);
        [CC.OUT.LDA] = ldbc(CC.MD,tmp);
        CC.OUT.LDAcl = CL((CC.OUT.LDA>0)+1);
        [CC.OUT.MDA] = mdbc(CC.MD,tmp);
        [tmp,ix] = min(CC.OUT.MDA,[],2);
        tmp = isnan(tmp);
        ix(tmp) = NaN;   % invalid output, not classified
        ix(~tmp) = CL(ix(~tmp));
        CC.OUT.MDAcl = ix;
end;

return;

d = JKD3 - JKD4;
tmp1 = d(1-min(T(:))+T(CC.TI,:),cl==CL(1));
[sum0,n0,ssq0] = sumskipnan(tmp1(:));       
tmp2 = d(1-min(T(:))+T(CC.TI,:),cl==CL(2));
[sum1,n1,ssq1]  = sumskipnan(tmp2(:));       
CC.MLL.AUC      = auc(tmp1,tmp2);
CC.MLL.ERR(1,1) = mean(sign([tmp1(:)]))/2+1/2;
CC.MLL.ERR(1,2) = mean(sign([tmp2(:)]))/2+1/2;
CC.MLL.ERR(2,1) = mean(sign([mean(tmp1,1)']))/2+1/2;
CC.MLL.ERR(2,2) = mean(sign([mean(tmp2,1)']))/2+1/2;
s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
SNR = 2*s./(s0+s1); % this is SNR+1 
CC.MLL.I   = log2(SNR)/2;
CC.MLL.SNR = SNR - 1;
if 0,
	clear tmp1 tmp2;
        tmp1 = stat2(d(:,cl==CL(1)),2);       
        tmp2 = stat2(d(:,cl==CL(2)),2);       
        CC.MLL.TSD=stat2res(tmp1,tmp2);
        CC.MLL.TSD.ERR=mean(sign([-d(:,cl==CL(1)),d(:,cl==CL(2))]),2)/2+1/2;
elseif bitand(SWITCH,1),
        CC.MLL.TSD=bci3eval(d(:,cl==CL(1)),d(:,cl==CL(2)),2);
end;

d = exp(-JKD1/2)-exp(-JKD2/2);
tmp1 = d(1-min(T(:))+T(CC.TI,:),cl==CL(1));
[sum0,n0,ssq0] = sumskipnan(tmp1(:));       
tmp2 = d(1-min(T(:))+T(CC.TI,:),cl==CL(2));
[sum1,n1,ssq1] = sumskipnan(tmp2(:));       
CC.GRB.AUC      = auc(tmp1,tmp2);
CC.GRB.ERR(1,1) = mean(sign([tmp1(:)]))/2+1/2;
CC.GRB.ERR(1,2) = mean(sign([tmp2(:)]))/2+1/2;
CC.GRB.ERR(2,1) = mean(sign([mean(tmp1,1)']))/2+1/2;
CC.GRB.ERR(2,2) = mean(sign([mean(tmp2,1)']))/2+1/2;
s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
SNR = 2*s./(s0+s1); % this is SNR+1 
CC.GRB.I   = log2(SNR)/2;
CC.GRB.SNR = SNR - 1;
if 0,
        clear tmp1 tmp2; 
        tmp1 = stat2(d(:,cl==CL(1)),2);       
        tmp2 = stat2(d(:,cl==CL(2)),2);       
        CC.GRB.TSD=stat2res(tmp1,tmp2);
        CC.GRB.TSD.ERR=1/2-mean(sign([-d(:,cl==CL(1)),d(:,cl==CL(2))]),2)/2;
elseif bitand(SWITCH,1),
        CC.GRB.TSD=bci3eval(d(:,cl==CL(1)),d(:,cl==CL(2)),2);
end;

