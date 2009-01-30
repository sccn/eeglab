function [CC,Q,tsd,md]=findclassifier2(D,TRIG,cl,T,t0,SWITCH)
% FINDCLASSIFIER2
%       is very similar to FINDCLASSIFIER1 but uses a different 
%       criterion for selecting the time segment.
% [CC,Q,TSD,MD]=findclassifier2(D,TRIG,Class,class_times,t_ref);
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
% [CC,Q,TSD,MD]=findclassifier2(AR,find(trig>0.5)-257,~mod(1:80,2),reshape(1:14*128,16,14*8)');
%
%
% Reference(s): 
% [1] Schlögl A., Neuper C. Pfurtscheller G.
% 	Estimating the mutual information of an EEG-based Brain-Computer-Interface
%  	Biomedizinische Technik 47(1-2): 3-8, 2002.
% [2] A. Schlögl, C. Keinrath, R. Scherer, G. Pfurtscheller,
%	Information transfer of an EEG-based Bran-computer interface.
%	Proceedings of the 1st International IEEE EMBS Conference on Neural Engineering, Capri, Italy, Mar 20-22, 2003 


%   $Id: findclassifier2.m,v 1.1 2009-01-30 06:04:49 arno Exp $
%   Copyright (C) 1999-2005 by Alois Schloegl <a.schloegl@ieee.org>	
%   This is part of the BIOSIG-toolbox http://biosig.sf.net/


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

cl = cl(:);
CL = unique(cl(~isnan(cl)));

[CL,iCL] = sort(CL);
TRIG = TRIG(:);
if ~all(D(:,1)==1)
        %        D1=[ones(size(D,1)-1,1),diff(D)];
        D = [ones(size(D,1),1),D];
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
        %[CC.QC3(k,1:2),d3{k}] = qcloglik({C{k,:}});
        CC.QC2(k) = sum(sqrt(d{k}(:)))/(size(d{k},1)*(size(d{k},1)-1));
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
                
                [tmp] = gdbc({C{k,:}},D(t(:),:));
                [tmp,ix] = min(tmp,[],2);
                tmp = isnan(tmp);
                ix(tmp) = NaN; %NC(1)+1;	% not classified; any value but not 1:length(MD)
                ix(~tmp) = CL(ix(~tmp));
                tmp = histo3([ix;CL(:)]);
                cmx2(iCL(tmp.X),l) = tmp.H-1;            
        end;
        CMX(k,:,:)   = cmx;
        CC.KAPPA(k)  = kappa(cmx);
        CC.ACC(k)    = sum(diag(cmx))/sum(cmx(:));
        CC.KAPPA2(k) = kappa(cmx2);
        CC.ACC2(k)   = sum(diag(cmx2))/sum(cmx2(:));

        t = perm(TRIG(~isnan(cl)),T(k,:));
        tmp = repmat(cl(~isnan(cl))',size(T,2),1);
        [tmp,acc_g(k,:),kap_g(k,:),H] = gdbc({C{k,:}},D(t(:),:),tmp(:));
end;	
CC.GDBC.acc = acc_g;
CC.GDBC.kap = kap_g;
CC.GDBC.H = H;
CC.GDBC.p = tmp;

% identify best classification time 
if nargin>4,
        tmp = CC.QC;
        tmp(~t0) = 0;
        [maxQ,CC.TI] = max(tmp); %d{K},
else
        [maxQ,CC.TI] = max(CC.QC); %d{K},
end;
CC.qc = [CC.QC',CC.KAPPA',CC.ACC',CC.QC2'];
%plot([CC.QC',CC.KAPPA',CC.ACC',CC.QC2',CC.QC3',t0]);drawnow
[maxQ(1),TI(1)] = max(CC.QC'.*t0); %d{K},
[maxQ(2),TI(2)] = max(CC.KAPPA'.*t0); %d{K},
[maxQ(3),TI(3)] = max(CC.ACC'.*t0); %d{K},
[maxQ(4),TI(4)] = max(CC.QC2'.*t0); %d{K},
%[maxQ(5),TI(5)] = max(CC.QC3(:,1).*t0); %d{K},
%[maxQ(6),TI(6)] = max(CC.QC3(:,2).*t0); %d{K},
CC.TIS = TI;
CC.TI  = TI(2);
%if any(TI(1)~=TI),
%        fprintf(1,'FINDCLASSIFIER2: \nTIS:  %i\t%i\t%i\t%i',TI);
%        fprintf(1,'\n     %5.3f\t%5.3f\t%5.3f\t%5.3f\n',maxQ);
%end;

% build MD classifier 
CC.D   = d{CC.TI};
CC.Q   = CC.QC(CC.TI);
CC.MD  = {C{CC.TI,:}};

if SWITCH<0, return; end; 

CC.IR  = mdbc({C{CC.TI,:}});
CC.CMX = squeeze(CMX(CC.TI,:,:));
CC.CMX = CMX;
CC.tsc = T(CC.TI,[1,end]);

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

% define datatype
CC.datatype = 'Classifier';
CC.Classes = CL;


%% cross-validation with jackknife (leave one trial out)
nc  = max(max(T))-min(min(T))+1;

JKD   = repmat(nan,[nc,length(CL),length(TRIG)]);
JKLL  = JKD;
JKD1  = repmat(nan,[nc,length(TRIG)]);
JKD2  = repmat(nan,[nc,length(TRIG)]);
JKD3  = repmat(nan,[nc,length(TRIG)]);
JKD4  = repmat(nan,[nc,length(TRIG)]);
JKD5  = repmat(nan,[nc,length(TRIG)]);
JKD6  = repmat(nan,[nc,length(TRIG)]);
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
        
        d = llbc(cc,D(t,:));
        JKLL(:,:,l) = d;
        [tmp,LLIX(:,l)] = max(d,[],2);
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
	end;
	                
        d = gdbc(cc,D(t,:));
        JKGD(:,:,l) = d;
        [tmp,GDIX(:,l)] = max(d,[],2);
        
        d = ldbc2(cc,D(t,:));
        LDA2(:,:,l) = d;
        [tmp,LD2IX(:,l)] = max(d,[],2);
        
        d = ldbc3(cc,D(t,:));
        LDA3(:,:,l) = d;
        size(LDA3),
        [tmp,LD3IX(:,l)] = max(d,[],2);
        
        d = ldbc4(cc,D(t,:));
        LDA4(:,:,l) = d;
        [tmp,LD4IX(:,l)] = max(d,[],2);
        
        if length(CL)==2,
                JKD5(:,l) = d(:,1);
                JKD6(:,l) = d(:,2);
                
                LDA(:,l) = ldbc(cc);
                JKLD(:,l) = D(t,:)*LDA(:,l);
        end;	
end;
if length(CL)==2,
        [CC.ldaC0,NN] = covm(LDA','D0');
        %CC.ldaC0=CC.ldaC0./NN*min(0,sum(~isnan(CL))-1); 
        % since NN==min(0,sum(~isnan(CL))-1), no need to rescale
end;	

% Concordance matrix with cross-validation 
CC.lmx= zeros([size(LLIX,1),length(CL)^2]);
CC.LLH.I0 = zeros([size(LLIX,1),length(CL)]);
CC.LLH.I  = zeros([size(LLIX,1),1]);

CC.mmx= zeros([size(MDIX,1),length(CL)^2]);
CC.MD2.I0 = zeros([size(MDIX,1),length(CL)]);
CC.MD2.I  = zeros([size(MDIX,1),1]);
tmp   = zeros([size(MDIX,1),length(CL)]);
for k = 1:length(CL),

        jkd = squeeze(LDA4(:,k,:));
        o = bci3eval(jkd(:,(cl~=k)&~isnan(cl)),jkd(:,cl==k),2);
        CC.LD4.TSD{k}  = o;
        CC.LD4.I0(:,k) = log2(2*var(jkd,[],2)./(var(jkd(:,cl==k),[],2) + var(jkd(:,(cl~=k)&~isnan(cl)),[],2)))/2;
        [sum0,n0,ssq0]=sumskipnan(jkd(:,cl==k),2);
        [sum1,n1,ssq1]=sumskipnan(jkd(:,(cl~=k)&~isnan(cl)),2);
        s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
        s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
        s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
        SNR = 2*s./(s0+s1); % this is SNR+1 
        CC.LD4.I1(:,k) = log2(SNR)/2;

        jkd = squeeze(LDA3(:,k,:));
        o = bci3eval(jkd(:,(cl~=k)&~isnan(cl)),jkd(:,cl==k),2);
        CC.LD3.TSD{k}  = o;
        CC.LD3.I0(:,k) = log2(2*var(jkd,[],2)./(var(jkd(:,cl==k),[],2) + var(jkd(:,(cl~=k)&~isnan(cl)),[],2)))/2;
        [sum0,n0,ssq0]=sumskipnan(jkd(:,cl==k),2);
        [sum1,n1,ssq1]=sumskipnan(jkd(:,(cl~=k)&~isnan(cl)),2);
        s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
        s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
        s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
        SNR = 2*s./(s0+s1); % this is SNR+1 
        CC.LD3.I1(:,k) = log2(SNR)/2;

        jkd = squeeze(LDA2(:,k,:));
        o = bci3eval(jkd(:,(cl~=k)&~isnan(cl)),jkd(:,cl==k),2);
        CC.LD2.TSD{k}  = o;
        CC.LD2.I0(:,k) = log2(2*var(jkd,[],2)./(var(jkd(:,cl==k),[],2) + var(jkd(:,(cl~=k)&~isnan(cl)),[],2)))/2;
        [sum0,n0,ssq0]=sumskipnan(jkd(:,cl==k),2);
        [sum1,n1,ssq1]=sumskipnan(jkd(:,(cl~=k)&~isnan(cl)),2);
        s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
        s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
        s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
        SNR = 2*s./(s0+s1); % this is SNR+1 
        CC.LD2.I1(:,k) = log2(SNR)/2;
        
        jkd = squeeze(JKGD(:,k,:));
        o = bci3eval(jkd(:,(cl~=k)&~isnan(cl)),jkd(:,cl==k),2);
        CC.MD3.TSD{k}  = o;
        CC.MD3.I0(:,k) = log2(2*var(jkd,[],2)./(var(jkd(:,cl==k),[],2) + var(jkd(:,(cl~=k)&~isnan(cl)),[],2)))/2;
        [sum0,n0,ssq0]=sumskipnan(jkd(:,cl==k),2);
        [sum1,n1,ssq1]=sumskipnan(jkd(:,(cl~=k)&~isnan(cl)),2);
        s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
        s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
        s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
        SNR = 2*s./(s0+s1); % this is SNR+1 
        CC.MD3.I1(:,k) = log2(SNR)/2;
        
        jkd = squeeze(JKD(:,k,:));
        o = bci3eval(jkd(:,(cl~=k)&~isnan(cl)),jkd(:,cl==k),2);
        CC.MD2.TSD{k}  = o;
        CC.MD2.I0(:,k) = log2(2*var(jkd,[],2)./(var(jkd(:,cl==k),[],2) + var(jkd(:,(cl~=k)&~isnan(cl)),[],2)))/2;
        [sum0,n0,ssq0]=sumskipnan(jkd(:,cl==k),2);
        [sum1,n1,ssq1]=sumskipnan(jkd(:,(cl~=k)&~isnan(cl)),2);
        s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
        s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
        s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
        SNR = 2*s./(s0+s1); % this is SNR+1 
        CC.MD2.I1(:,k) = log2(SNR)/2;
        
        jkd = squeeze(JKD(:,k,:).^2);
        o = bci3eval(jkd(:,(cl~=k)&~isnan(cl)),jkd(:,cl==k),2);
        CC.MDA.TSD{k}  = o;
        CC.MDA.I0(:,k) = log2(2*var(jkd,[],2)./(var(jkd(:,cl==k),[],2) + var(jkd(:,(cl~=k)&~isnan(cl)),[],2)))/2;
        [sum0,n0,ssq0]=sumskipnan(jkd(:,cl==k),2);
        [sum1,n1,ssq1]=sumskipnan(jkd(:,(cl~=k)&~isnan(cl)),2);
        s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
        s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
        s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
        SNR = 2*s./(s0+s1); % this is SNR+1 
        CC.MDA.I1(:,k) = log2(SNR)/2;
        
        jkd = exp(-squeeze(JKD(:,k,:))/2);
        o = bci3eval(jkd(:,(cl~=k)&~isnan(cl)),jkd(:,cl==k),2);
        CC.GRB.TSD{k}  = o;
        CC.GRB.I0(:,k) = log2(2*var(jkd,[],2)./(var(jkd(:,cl==k),[],2) + var(jkd(:,(cl~=k)&~isnan(cl)),[],2)))/2;
        [sum0,n0,ssq0] = sumskipnan(jkd(:,cl==k),2);
        [sum1,n1,ssq1] = sumskipnan(jkd(:,(cl~=k)&~isnan(cl)),2);
        s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
        s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
        s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
        SNR = 2*s./(s0+s1); % this is SNR+1 
        CC.GRB.I1(:,k) = log2(SNR)/2;
        
        jkd = exp(-squeeze(JKLL(:,k,:))/2);
        o = bci3eval(jkd(:,(cl~=k)&~isnan(cl)),jkd(:,cl==k),2);
        CC.LLH.TSD{k}  = o;
        CC.LLH.I0(:,k) = log2(2*var(jkd,[],2)./(var(jkd(:,cl==k),[],2) + var(jkd(:,(cl~=k)&~isnan(cl)),[],2)))/2;
        [sum0,n0,ssq0] = sumskipnan(jkd(:,cl==k),2);
        [sum1,n1,ssq1] = sumskipnan(jkd(:,(cl~=k)&~isnan(cl)),2);
        s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
        s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
        s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
        SNR = 2*s./(s0+s1); % this is SNR+1 
        CC.LLH.I1(:,k) = log2(SNR)/2;
        
        for l = 1:length(CL),
                tmp(:,l) = sum(MDIX(:,cl==CL(k))==CL(l),2);    
                if CL(k) == CL(l),
                        acc = tmp(:,l);
                end;
        end;
        CC.MDA.mmx(:,(1-length(CL):0)+k*length(CL)) = tmp;
        CC.MDA.acc(:,k) = acc./sum(tmp,2);

        for l = 1:length(CL),
                tmp(:,l) = sum(GDIX(:,cl==CL(k))==CL(l),2);    
                if CL(k) == CL(l),
                        acc = tmp(:,l);
                end;
        end;
        CC.MD3.mmx(:,(1-length(CL):0)+k*length(CL)) = tmp;
        CC.MD3.acc(:,k) = acc./sum(tmp,2);

        for l = 1:length(CL),
                tmp(:,l) = sum(LD2IX(:,cl==CL(k))==CL(l),2);    
                if CL(k) == CL(l),
                        acc = tmp(:,l);
                end;
        end;
        CC.LD2.mmx(:,(1-length(CL):0)+k*length(CL)) = tmp;
        CC.LD2.acc(:,k) = acc./sum(tmp,2);

        for l = 1:length(CL),
                tmp(:,l) = sum(LD3IX(:,cl==CL(k))==CL(l),2);    
                if CL(k) == CL(l),
                        acc = tmp(:,l);
                end;
        end;
        CC.LD3.mmx(:,(1-length(CL):0)+k*length(CL)) = tmp;
        CC.LD3.acc(:,k) = acc./sum(tmp,2);

        for l = 1:length(CL),
                tmp(:,l) = sum(LD4IX(:,cl==CL(k))==CL(l),2);    
                if CL(k) == CL(l),
                        acc = tmp(:,l);
                end;
        end;
        CC.LD4.mmx(:,(1-length(CL):0)+k*length(CL)) = tmp;
        CC.LD4.acc(:,k) = acc./sum(tmp,2);

        for l = 1:length(CL),
                tmp(:,l) = sum(LLIX(:,cl==CL(k))==CL(l),2);    
                if CL(k) == CL(l),
                        acc = tmp(:,l);
                end;
        end;
        CC.LLH.mmx(:,(1-length(CL):0)+k*length(CL)) = tmp;
        CC.LLH.acc(:,k) = acc./sum(tmp,2);
end;
CC.MDA.I = sum(CC.MDA.I0,2);
CC.MD2.I = sum(CC.MD2.I0,2);
CC.MD3.I = sum(CC.MD3.I0,2);
CC.LD2.I = sum(CC.LD2.I0,2);
CC.GRB.I = sum(CC.GRB.I0,2);
CC.LLH.I = sum(CC.LLH.I0,2);

CC.MDA.CMX00 = reshape(sum(CC.MDA.mmx(T(CC.TI,:),:),1),[1,1]*length(CL))/size(T,2);
CC.MDA.ACC00 = sum(CC.MDA.mmx(:,1:length(CL)+1:end),2)/sum(~isnan(cl));	
CC.MDA.KAP00 = zeros(size(MDIX,1),1);

CC.MD3.CMX00 = reshape(sum(CC.MD3.mmx(T(CC.TI,:),:),1),[1,1]*length(CL))/size(T,2);
CC.MD3.ACC00 = sum(CC.MD3.mmx(:,1:length(CL)+1:end),2)/sum(~isnan(cl));	
CC.MD3.KAP00 = zeros(size(GDIX,1),1);

CC.LD2.CMX00 = reshape(sum(CC.LD2.mmx(T(CC.TI,:),:),1),[1,1]*length(CL))/size(T,2);
CC.LD2.ACC00 = sum(CC.LD2.mmx(:,1:length(CL)+1:end),2)/sum(~isnan(cl));	
CC.LD2.KAP00 = zeros(size(LD2IX,1),1);

CC.LD3.CMX00 = reshape(sum(CC.LD3.mmx(T(CC.TI,:),:),1),[1,1]*length(CL))/size(T,2);
CC.LD3.ACC00 = sum(CC.LD3.mmx(:,1:length(CL)+1:end),2)/sum(~isnan(cl));	
CC.LD3.KAP00 = zeros(size(LD3IX,1),1);

CC.LD4.CMX00 = reshape(sum(CC.LD4.mmx(T(CC.TI,:),:),1),[1,1]*length(CL))/size(T,2);
CC.LD4.ACC00 = sum(CC.LD4.mmx(:,1:length(CL)+1:end),2)/sum(~isnan(cl));	
CC.LD4.KAP00 = zeros(size(LD4IX,1),1);

CC.LLH.CMX00 = reshape(sum(CC.LLH.mmx(T(CC.TI,:),:),1),[1,1]*length(CL))/size(T,2);
CC.LLH.ACC00 = sum(CC.LLH.mmx(:,1:length(CL)+1:end),2)/sum(~isnan(cl));	
CC.LLH.KAP00 = zeros(size(LLIX,1),1);
for k = 1:size(MDIX,1),
        CC.MDA.KAP00(k) = kappa(reshape(CC.MDA.mmx(k,:),[1,1]*length(CL)));
        CC.MD3.KAP00(k) = kappa(reshape(CC.MD3.mmx(k,:),[1,1]*length(CL)));
        CC.LD2.KAP00(k) = kappa(reshape(CC.LD2.mmx(k,:),[1,1]*length(CL)));
        CC.LD3.KAP00(k) = kappa(reshape(CC.LD3.mmx(k,:),[1,1]*length(CL)));
        CC.LD4.KAP00(k) = kappa(reshape(CC.LD4.mmx(k,:),[1,1]*length(CL)));
        CC.LLH.KAP00(k) = kappa(reshape(CC.LLH.mmx(k,:),[1,1]*length(CL)));
end;

if length(CL) > 2, 
        return; 
end; 


if bitand(SWITCH,1),
        CC.LLH.ERR00 = (mean(sign(JKLL),2)+1)/2;
        CC.LDA.ERR00 = (mean(sign(JKLD),2)+1)/2;
        CC.MD3.ERR00 = (mean(sign(JKGD),2)+1)/2;
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

d = JKD6 - JKD5;
tmp1 = d(1-min(T(:))+T(CC.TI,:),cl==CL(1));
[sum0,n0,ssq0] = sumskipnan(tmp1(:));       
tmp2 = d(1-min(T(:))+T(CC.TI,:),cl==CL(2));
[sum1,n1,ssq1] = sumskipnan(tmp2(:));       
CC.MD3.AUC      = auc(tmp1,tmp2);
CC.MD3.ERR(1,1) = mean(sign([tmp1(:)]))/2+1/2;
CC.MD3.ERR(1,2) = mean(sign([tmp2(:)]))/2+1/2;
CC.MD3.ERR(2,1) = mean(sign([mean(tmp1,1)']))/2+1/2;
CC.MD3.ERR(2,2) = mean(sign([mean(tmp2,1)']))/2+1/2;
s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
SNR = 2*s./(s0+s1); % this is SNR+1 
CC.MD3.I   = log2(SNR)/2;
CC.MD3.SNR = SNR - 1;
if 0,
        clear tmp1 tmp2; 
        tmp1 = stat2(d(:,cl==CL(1)),2);       
        tmp2 = stat2(d(:,cl==CL(2)),2);       
        CC.MD3.TSD=stat2res(tmp1,tmp2);
        CC.MD3.TSD.ERR=1/2-mean(sign([-d(:,cl==CL(1)),d(:,cl==CL(2))]),2)/2;
elseif bitand(SWITCH,1),
        CC.MD3.TSD=bci3eval(d(:,cl==CL(1)),d(:,cl==CL(2)),2);
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
        ix(tmp)  = NaN;   % invalid output, not classified
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

