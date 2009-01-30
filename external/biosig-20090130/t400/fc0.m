function [CC,Q,tsd,md,cc]=fc0(D,TRIG,T,arg4)
% FC finds a classifier for asnychroneous data
% 
% [CC,Q,TSD,MD]=fc0(D,TRIG,class_times [,SWITCH]);
%
% D 	data, each row is one time point
% TRIG	trigger time points
% class_times each row defines a segment used for Classification
% SWITCH 0 [default] minimizes result
%	2 only have of the possible t1-t2 combinations are used
%
% CC 	contains LDA and MD classifiers
% Q  	is a list of classification quality for each time of 'class_times'
% TSD 	returns the LDA classification 
% MD	returns the MD  classification 
%
%
% [CC,Q,TSD,MD]=fc(AR,find(trig>0.5)-257,reshape(1:14*128,16,14*8)');
%
% see also: COVM.M, QCMAHAL, MDBC, LDBC, 

%
%	$Revision: 1.1 $
%	$Id: fc0.m,v 1.1 2009-01-30 06:04:48 arno Exp $
%	Copyright (c) 1999-2002, 2004 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/
%

% CHANGELOG
%	10.12.2001	changed Covariance-Matrix from unbiased to biased
%	20.12.2001	ROC and AUC included 	
% 	04.01.2002	zeros class added 
%			md changed to diff of log of MD
%	18.01.2002	Gaussian Radial basis functions
%	12.02.2002	included in BCI7-paradigm
%	13.02.2002	CC.D included
%       30.04.2002	JACKKNIFE included
%	19.07.2002	SWITCH included 

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

if nargin>3;
        SWITCH = arg4;
end;

TRIG=TRIG(:);
if ~all(D(:,1)==1)
%        D1=[ones(size(D,1)-1,1),diff(D)];
        D =[ones(size(D,1),1),D];
%else
%	D1=[ones(size(D,1)-1,1),diff(D(:,2:end))];
end;

dT = (min(T(:)):max(T(:)));
nc      = length(dT);
%%% pad sufficient NaN's
off=min(TRIG(1)+min(T(:))-1,0);
D=[repmat(nan,-off,size(D,2));D];
TRIG=TRIG-off;
off=max(max(T(:))+TRIG(length(TRIG))-length(D),0);
D=[D;repmat(nan,off,size(D,2))];


[CC0,NN0] = covm(D,'E');
for k = 1:size(T,1),	
        t = perm(TRIG,T(k,:));
        t = t(:);
        [cc{k},nn{k}] = covm(D(t,:),'M');
        %CC{k}=cc{k}./nn{k};
end;
[Q , d] = qcmahal(cc);
if bitand(SWITCH,1)
        CC.d=d;
end;
if bitand(SWITCH,2),
        d=d+d';
        d=d(:,1:size(d,2)/2);
end;

[ix,iy] = find(d == max(d(:)));

ix = ix(1); iy = iy(1);
if ix>iy,
	tmp=ix;ix=iy;iy=tmp;
end;
CC.TI   = [ix,iy];
CC.MD   = {cc{ix},cc{iy}};
[CC.Q,CC.D]= qcmahal(CC.MD);
Q=CC.Q;


%CC.IR   = mdbc(CC.MD);

JKD1=zeros(length(dT),length(TRIG));
JKD2=zeros(length(dT),length(TRIG));
%CC.SGN = zeros(length(TRIG),nc);
%% jackknife
for l=1:length(TRIG),
        T1 = TRIG(l) + T(CC.TI(1),:);
        T2 = TRIG(l) + T(CC.TI(2),:);
        
        tmp = D(T1,:);
        [tmp1,tmp2]=covm(tmp,'M');
        cc1 = CC.MD{1}-tmp1;
        
        tmp = D(T2,:);
        [tmp1,tmp2]=covm(tmp,'M');
        cc2 = CC.MD{2}-tmp1;

        t  = TRIG(l) + dT;
        d = mdbc({cc1,cc2},D(t,:));
        JKD1(:,l)=d(:,1);
        JKD2(:,l)=d(:,2);
        
        d = llbc({cc1,cc2},D(t,:));
        JKD3(:,l)=d(:,1);
        JKD4(:,l)=d(:,2);
        
        JKLD(:,l) = ldbc({cc1,cc2}, D(t,:));
end;
if bitand(SWITCH,1),
        CC.LDA.ERR00 = (mean(sign(JKLD),2)+1)/2;
        CC.MDA.ERR00 = (mean(sign(JKD1-JKD2),2)+1)/2;
        CC.GRB.ERR00 = (mean(sign(exp(-JKD2/2)-exp(-JKD1/2)),2)+1)/2;
end;

d = JKD3 - JKD4;
tmp1 = d(1-min(T(:))+T(CC.TI(1),:),:);
[sum0,n0,ssq0] = sumskipnan(tmp1(:));       
tmp2 = d(1-min(T(:))+T(CC.TI(2),:),:);
[sum1,n1,ssq1] = sumskipnan(tmp2(:));       
CC.MLL.AUC      = auc(tmp1,tmp2);
CC.MLL.ERR(1,:) = mean(sign([tmp1(:),tmp2(:)]))/2+1/2;
CC.MLL.ERR(2,:) = mean(sign([mean(tmp1,1)',mean(tmp2,1)']))/2+1/2;
s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
SNR = 2*s./(s0+s1); % this is SNR+1 
CC.MLL.I   = log2(SNR)/2;
CC.MLL.SNR = SNR - 1;
clear tmp; tmp.datatype='STAT2';
[tmp.SUM,tmp.N,tmp.SSQ] = sumskipnan(d,2);       
if bitand(SWITCH,1),
        CC.MLL.TSD=tmp;
end;

d = JKD1 - JKD2;
tmp1 = d(1-min(T(:))+T(CC.TI(1),:),:);
[sum0,n0,ssq0] = sumskipnan(tmp1(:));       
tmp2 = d(1-min(T(:))+T(CC.TI(2),:),:);
[sum1,n1,ssq1] = sumskipnan(tmp2(:));       
CC.MDA.AUC      = auc(tmp1,tmp2);
CC.MDA.ERR(1,:) = mean(sign([tmp1(:),tmp2(:)]))/2+1/2;
CC.MDA.ERR(2,:) = mean(sign([mean(tmp1,1)',mean(tmp2,1)']))/2+1/2;
s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
SNR = 2*s./(s0+s1); % this is SNR+1 
CC.MDA.I   = log2(SNR)/2;
CC.MDA.SNR = SNR - 1;
clear tmp; tmp.datatype='STAT2';
[tmp.SUM,tmp.N,tmp.SSQ] = sumskipnan(d,2);       
if bitand(SWITCH,1),
        CC.MDA.TSD=tmp;
end;

d = exp(-JKD1/2)-exp(-JKD2/2);
tmp1 = d(1-min(T(:))+T(CC.TI(1),:),:);
[sum0,n0,ssq0] = sumskipnan(tmp1(:));       
tmp2 = d(1-min(T(:))+T(CC.TI(2),:),:);
[sum1,n1,ssq1] = sumskipnan(tmp2(:));       
CC.GRB.AUC      = auc(tmp1,tmp2);
CC.GRB.ERR(1,:) = mean(sign([tmp1(:),tmp2(:)]))/2+1/2;
CC.GRB.ERR(2,:) = mean(sign([mean(tmp1,1)',mean(tmp2,1)']))/2+1/2;
s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
SNR = 2*s./(s0+s1); % this is SNR+1 
CC.GRB.I   = log2(SNR)/2;
CC.GRB.SNR = SNR - 1;
clear tmp; tmp.datatype='STAT2';
[tmp.SUM,tmp.N,tmp.SSQ] = sumskipnan(d,2);       
if bitand(SWITCH,1),
        CC.GRB.TSD=tmp;
end;

d=JKLD;
tmp1 = d(1-min(T(:))+T(CC.TI(1),:),:);
[sum0,n0,ssq0] = sumskipnan(tmp1(:));       
tmp2 = d(1-min(T(:))+T(CC.TI(2),:),:);
[sum1,n1,ssq1] = sumskipnan(tmp2(:));       
CC.LDA.AUC      = auc(tmp1,tmp2);
CC.LDA.ERR(1,:) = mean(sign([tmp1(:),tmp2(:)]))/2+1/2;
CC.LDA.ERR(2,:) = mean(sign([mean(tmp1,1)',mean(tmp2,1)']))/2+1/2;
s0  = (ssq0-sum0.*sum0./n0)./(n0-1);
s1  = (ssq1-sum1.*sum1./n1)./(n1-1);
s   = (ssq0+ssq1-(sum0+sum1).*(sum0+sum1)./(n0+n1))./(n0+n1-1);
SNR = 2*s./(s0+s1); % this is SNR+1 
CC.LDA.I   = log2(SNR)/2;
CC.LDA.SNR = SNR - 1;
clear tmp; tmp.datatype='STAT2';
[tmp.SUM,tmp.N,tmp.SSQ] = sumskipnan(d,2);       
if bitand(SWITCH,1),
        CC.LDA.TSD=tmp;
end;




return

%%%% reference intervall
if nargin>4, 
	T = perm(TRIG,t(k,:));
        T = T(T<=size(D,1));
        tmp = D(T(:),:);
        C0r = tmp'*tmp;	
end;

for k = 1:size(t,1),
        for l = 1:length(CL), 
                T = perm(TRIG(cl==CL(l)),t(k,:));
                T = T(T<=size(D,1));
                tmp = D(T(:),:);
		C{k,l} = tmp'*tmp;   	
        end;
        %[Q(k),d{k}] = qcmahal({C0r,C{k,:}});
        [Q(k),d{k}] = qcmahal({C{k,:}});
        lnQ(k) = mean(log(d{k}(~eye(length(d{k})))));
end;
[maxQ,CC.TI] = max(Q); %d{K},
%CC.TI = K;
CC.MD = {C{CC.TI,:}};
CC.IR = mdbc({C{CC.TI,:}});
CC.D = d{CC.TI};
Q = Q(CC.TI);

[maxQ,CC.lnTI] = max(lnQ); %d{K},
CC.DistMXln = d{CC.lnTI};
CC.MDln = {C{CC.lnTI,:}};

% LDA 
C0 = zeros(size(C{CC.TI,1}));
for l=1:length(CL);
        [M{l},sd,COV,xc,N,R2] = decovm(C{CC.TI,l});
	C0 = C0 + C{CC.TI,l};
end;
[M0,sd,COV0,xc,N,R2] = decovm(C0);
w     = COV0\(M{1}'-M{2}');
w0    = M0*w;
CC.LDA.b = w0;
CC.LDA.w = -w;
CC.lda = [w0; -w];

% MD 
md = mdbc(CC.MD,D);


[tmp,IX] = min(md,[],2);
for k = 1:size(t,1),
	H0{k} = zeros(length(CL));
	for l = 1:size(t,2),
		T = TRIG+t(k,l);
	        T = T(T<=size(D,1));
		[tmp,sd,H] = kappa(cl(:)-min(cl)+1,IX(T),length(CL)); 	
		H0{k} = H0{k} + H;
	end;
	[CC.KAPPA(k),se,tmp,z,CC.OA(k),CC.SA(k,1:2)] = kappa(H0{k});
end;
CC.H0=H0;

if length(CL)>2, return; end; 

CC.LDA=eval_offline(D*CC.lda,TRIG,cl);
CC.MDA=eval_offline(md(:,1)-md(:,2),TRIG,cl);
CC.GRB=eval_offline(exp(-md(:,2)/2)-exp(-md(:,1)/2),TRIG,cl);
CC.GR1=eval_offline(exp(-md(:,2))-exp(-md(:,1)),TRIG,cl);

nc  = ceil(max(diff(TRIG))*1.5);
ERR00md = zeros(1,nc);
ERR00gr = zeros(1,nc);
ERR00ld = zeros(1,nc);
N00     = zeros(1,nc);

%% jackknife
for l=1:length(cl),
        c = find(cl(l)==CL);
        T = TRIG(l)+t(CC.TI,:);
        T = T(T<=size(D,1));
        tmp = D(T(:),:);
        cc{c}   = CC.MD{c}-tmp'*tmp;
        cc{3-c} = CC.MD{3-c};
        T = TRIG(l)+(1:nc);
        T = T(T<=size(D,1));
        [tmp] = mdbc(cc,D(T,:));
        ERR00md(1:length(T))= ERR00md(1:length(T))+sign(1.5-c)*sign(tmp(:,2)-tmp(:,1))';                
        ERR00gr(1:length(T))= ERR00gr(1:length(T))+sign(1.5-c)*sign(exp(-tmp(:,1)/2)-exp(-tmp(:,2)/2))';                
        tmp = ldbc(cc,D(T,:));
        ERR00ld(1:length(T))= ERR00ld(1:length(T))+sign(c-1.5)*sign(tmp)';                
        N00(1:length(T))  = N00(1:length(T))+1;                
end;

CC.ERR00md   = (1-ERR00md./N00)/2;
CC.ERR00gr   = (1-ERR00gr./N00)/2;
CC.ERR00ld   = (1-ERR00ld./N00)/2;
CC.N00     = N00;

