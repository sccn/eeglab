function [cc02,cc32]=make_cc7(fn,eegchan,trigchan,Fs)
% Builds a classifier based on AAR parameters from BCI recordings. 
%  The result includes also initial values for the next run. 
%
% CC=make_cc(fn,eegchan,trigchan,Fs)
%  
%  e.g. 
%  CC=make_cc('x21fb*')
%  CC=make_cc({'x21fb1.mat','x21fb2.mat'})
%  CC=make_cc(fn,[1,3],4,128)
%
% default: 
% 	eegchan=[1,3];
% 	trigchan=4;
% 	Fs=128;
%
% References:
% [1] Schlögl A., Neuper C. Pfurtscheller G. (2002)
%   Estimating the mutual information of an EEG-based Brain-Computer-Interface.
%   Biomedizinische Technik 47(1-2): 3-8, 2002
% [2] A. Schlögl, C. Keinrath, R. Scherer, G. Pfurtscheller, (2003)
%   Information transfer of an EEG-based Bran-computer interface.
%   Proceedings of the 1st International IEEE EMBS Conference on Neural Engineering, Capri, Italy, Mar 20-22, 2003. 
% [3] A. Schlögl, D. Flotzinger and G. Pfurtscheller (1997)
%   Adaptive Autoregressive Modeling used for Single-trial EEG Classification
%   Biomedizinische Technik 42: (1997), 162-167. 
% [4] A. Schlögl, C. Neuper and G. Pfurtscheller (1997)
%   Subject-specific EEG pattern during motor imagery
%   Proceedings of the 19th Annual International Conference if the IEEE Engineering in Medicine and Biology Society , vol 19, pp.1530-1532, 1997.
% [5] A. Schlögl, K. Lugger and G. Pfurtscheller (1997)
%   Using Adaptive Autoregressive Parameters for a Brain-Computer-InterfaceExperiment,
%   Proceedings of the 19th Annual International Conference if the IEEE Engineering in Medicine and Biology Society ,vol 19 , pp.1533-1535, 1997.
%

%	$Revision: 1.1 $
%	$Id: make_cc_tc2004.m,v 1.1 2009-01-30 06:04:39 arno Exp $
%	Copyright (C) 1999-2004 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

cname = computer;

if nargin<2,
	eegchan  = NaN;
end;
if nargin<3,
	trigchan = 4;
        if ~isempty(findstr(fn,'Bci7c/'))
                trigchan = 8;
        end;
end;
if nargin<4,
	Fs = 128;
end;
 
M0  = 7;
MOP = 3;
uc  = 30:5:80;

[S,H] = sload(fn);

%S = regress_eog(S,1:3,4:7); disp('!!! EOG removal with regression  !!!');

Fs    = H.SampleRate;
%ix = (H.EVENT.TYP>hex2dec('300')) & (H.EVENT.TYP<hex2dec('30d')); 
%H.Classlabel = mod(H.EVENT.TYP(ix),256);
TRIG = [];
if isfield(H,'EVENT'),  % GDF, BCI8 data
        TRIG = H.EVENT.POS(H.EVENT.TYP==hex2dec('300'));
        if isempty(TRIG),       % BKR, BCI7x data
                TRIG = H.EVENT.POS((H.EVENT.TYP>hex2dec('300'))&(H.EVENT.TYP<hex2dec('305')));
        end;
end;
if isempty(TRIG),       % BKR, BCI7x data
        TRIG = gettrigger(S(:,trigchan));
        if isfield(H,'TriggerOffset');
                TRIG = TRIG - H.TriggerOffset*Fs/1000;
                LEN  = H.BCI.Paradigm.TrialDuration*Fs/1000;
        elseif ~isempty(findstr(lower(fn),'bci7a')),    
                fprintf(1,'MAKE_CC7: assume trigger occures at t=2s \n');
                TRIG = TRIG - 3*Fs;
                LEN  = 9*Fs;
        else
                fprintf(1,'MAKE_CC7: assume trigger occures at t=2s \n');
                TRIG = TRIG - 2*Fs;
                LEN  = 9*Fs;
        end;
end;
cl = H.Classlabel(:);

if isnan(eegchan)
        eegchan = [1,3];
end;

if length(TRIG)~=length(cl);
        fprintf(2,'Attention: \n'),
end;
if 0, isfield(H,'ArtifactSelection'),
        for k = find(H.ArtifactSelection)';
                S(TRIG(k)+(1:9*Fs),eegchan) = NaN;
        end;
        
        TRIG = TRIG(~H.ArtifactSelection);
        cl   = cl  (~H.ArtifactSelection);
end;

if length(TRIG)~=length(cl);
        fprintf(2,'number of Triggers (%i) does not fit size of class information (%i)',length(TRIG),length(cl));
	return;        
end;

if ~any(size(eegchan)==1)
	S = S(:,1:size(eegchan,1))*eegchan;
	eegchan=1:size(eegchan,2); 
end;


if Fs==250,
        Fs = Fs/2;
        S  = rs(S,2,1);        
        TRIG = round(TRIG/2);
end;

% Muscle Detection 
% [INI,s,E,arti] = detectmuscle2(S(:,1:3));

% find initial values 
randn('seed',0);
[a0,A0] = getar0(S(:,eegchan),1:15,1000,Fs/2);
%save o3_init_A0 a0 A0
%load o3_init_A0 a0 A0
tmp = ceil(9*Fs/16)*16;
T  = reshape((1:tmp),16,tmp/16)';
t0 = zeros(tmp/16,1);
if ~isfield(H,'BCI')
        t0(33:60) = 1;   % look for interval [4s,7.5s]
else
        if isfield(H.BCI.Paradigm,'FallingPeriod')          % Basket Paradigm 
                T1 = H.BCI.Paradigm.FallingPeriod;
        elseif isfield(H.BCI.Paradigm,'FeedbackTiming')
                T1 = H.BCI.Paradigm.FeedbackTiming;
        end;
        %T1 = T1/1000*8;        % wrong !!!!
        T1 = T1/1000*Fs;        % correct !!!
        t0((T(:,1)>T1(1)) & (T(:,1)<T1(2))) = 1;
end;

t0 = logical(t0);

for p = 3; 6;1:M0;
for k = 7; 1:length(uc);
	UC0 = 2^(-uc(k)/8);
        ar = zeros(size(S,1),p*length(eegchan));
        e  = zeros(size(S,1),  length(eegchan));
        for ch = 1:length(eegchan),
                [ar(:,(1-p:0)+p*ch),e(:,ch),REV(ch)] = aar(S(:,eegchan(ch)), [2,3], p, UC0, a0{p},A0{p});
        end;

        [cc,Q,tsd,md] = findclassifier2(ar,TRIG, cl,T,t0,3);
        cc.T = (min(T(:)):max(T(:)))'/Fs;
        cc.MDA.TSD.T = (min(T(:)):max(T(:)))'/Fs;
        cc.MD2.TSD.T = (min(T(:)):max(T(:)))'/Fs;
        cc.MD3.TSD.T = (min(T(:)):max(T(:)))'/Fs;
        cc.LDA.TSD.T = (min(T(:)):max(T(:)))'/Fs;
        cc.GRB.TSD.T = (min(T(:)):max(T(:)))'/Fs;
        
        tmp   = ar(~any(isnan(e),2),:);
        cc.C0 = covm(tmp,'E');
        cc.W0 = covm(diff(tmp),'M');
        cc.V  = meansq(e);
        for ch = 1:length(eegchan),
                li        = (1-p:0)+ch*p;
                cc.a0{ch} = cc.C0(1,li+1);
                cc.A0{ch} = cc.C0(li+1,li+1) - cc.a0{ch}'*cc.a0{ch};
	        cc.W{ch}  = cc.W0(li,li);
        end;
        
        cc.Method = 'aar';
        cc.Fs  = Fs;
        cc.UC  = UC0;
        cc.p   = p;
        cc.MOP = p;
        cc.REV = REV;
        %cc.a0  = a0{p};
        %cc.A0  = A0{p};
        cc.EEGCHAN = eegchan;
        %cc.TRIGCHAN = trigchan;
        %CC{k+(p-1)*length(uc)}=cc;
        cc02 = cc; 

end;
end;

return

for p = 3;6;1:M0;
for k = 7;6;1:length(uc);
	UC0 = 2^(-uc(k)/8);
        ar = zeros(size(S,1),p*length(eegchan));
        e  = zeros(size(S,1),  length(eegchan));
        for ch = 1:length(eegchan),
                [ar(:,(1-p:0)+p*ch),e(:,ch),REV(ch)] = aar(S(:,eegchan(ch)), [0,0], p, UC0, cc0.a0{ch},cc0.A0{ch},cc0.W{ch},cc0.V(ch));
        end;
        
        [cc,Q,tsd,md] = findclassifier2(ar,TRIG, cl,T,t0,3);
        cc.T = (min(T(:)):max(T(:)))'/Fs;
        cc.MDA.TSD.T = (min(T(:)):max(T(:)))'/Fs;
        cc.MD2.TSD.T = (min(T(:)):max(T(:)))'/Fs;
        cc.MD3.TSD.T = (min(T(:)):max(T(:)))'/Fs;
        cc.LDA.TSD.T = (min(T(:)):max(T(:)))'/Fs;
        cc.GRB.TSD.T = (min(T(:)):max(T(:)))'/Fs;
        
        tmp   = ar(~any(isnan(e),2),:);
        cc.C0 = covm(tmp,'E');
        cc.W0 = covm(diff(tmp),'M');
        cc.V  = meansq(e);
        for ch = 1:length(eegchan),
                li        = (1-p:0)+ch*p;
                cc.a0{ch} = cc.C0(1,li+1);
                cc.A0{ch} = cc.C0(li+1,li+1) - cc.a0{ch}'*cc.a0{ch};
	        cc.W{ch}  = cc.W0(li,li);
        end;
        
        cc.Method = 'aar';
        cc.Fs  = Fs;
        cc.UC  = UC0;
        cc.p   = p;
        cc.MOP = p;
        cc.REV = REV;
        %cc.a0  = a0{p};
        %cc.A0  = A0{p};
        cc.EEGCHAN = eegchan;
        %cc.TRIGCHAN = trigchan;
        %CC{k+(p-1)*length(uc)}=cc;
        cc32 = cc;
end;
end;
