function QRS=nqrsdetect(S,fs);
% nqrsdetect - detection of QRS-complexes
%
%   QRS=nqrsdetect(S,fs);
%
% INPUT
%   S       ecg signal data
%   fs      sample rate
%
% OUTPUT
%   QRS     fiducial points of qrs complexes
%
%
% see also: QRSDETECT
%
% Reference(s):
% [1]: V. Afonso, W. Tompkins, T. Nguyen, and S. Luo, "ECG beat detection using filter banks,"
% 	IEEE Trans. Biomed. Eng., vol. 46, no. 2, pp. 192--202, Feb. 1999
%
% [2]: A.V. Oppenheim, R.W. Schafer, and J.R. Buck,  Discrete-Time Signal
% 	Processing, second edition, Prentice Hall, 1999, chapter 4.7.3

% Copyright (C) 2006 by Rupert Ortner
%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful, ...
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
%% USA


S=S(:);
S=full(S);
N=round(fs);   %Filter order
%---------------------------------------
%Replaces filter bank in [1]
Bw=5.6;     %filter bandwidth
Bwn=1/(fs/2)*Bw;    
M=round((fs/2)/Bw); %downsampling rate

Wn0=Bwn;    %bandwidth of the first filter
Wn1=[Bwn 2*Bwn];    %bandwidth of the second filter
Wn2=[2*Bwn 3*Bwn];
Wn3=[3*Bwn 4*Bwn];
Wn4=[4*Bwn 5*Bwn];

h0=fir1(N,Wn0); %impulse response of the first filter
h1=fir1(N,Wn1,'bandpass');
h2=fir1(N,Wn2,'bandpass');
h3=fir1(N,Wn3,'bandpass');
h4=fir1(N,Wn4,'bandpass');

%Polyphase implementation of the filters
y=cell(1,5);    
y{1}=polyphase_imp(S,h0,M); %W0 (see [1]) filtered and downsampled signal
y{2}=polyphase_imp(S,h1,M); %W1 
y{3}=polyphase_imp(S,h2,M); %W2
y{4}=polyphase_imp(S,h3,M); %W3
y{5}=polyphase_imp(S,h4,M); %W4
%----------------------------------------------

cut=ceil(N/M);  %Cutting off of initial transient because of the filtering
y1=[zeros(cut,1);y{1}(cut:length(y{1}))];  
y2=[zeros(cut,1);y{2}(cut:length(y{2}))];
y3=[zeros(cut,1);y{3}(cut:length(y{3}))];
y4=[zeros(cut,1);y{4}(cut:length(y{4}))];
y5=[zeros(cut,1);y{5}(cut:length(y{5}))];
%----------------------------------------

P1=sum([abs(y2) abs(y3) abs(y4)],2); %see [1] equation (13)
P2=sum([abs(y2) abs(y3) abs(y4) abs(y5)],2);
P4=sum([abs(y3) abs(y4) abs(y5)],2);

FL1=MWI(P1); %Feature 1 according to Level 1 in [1]
FL2=MWI(P2); %Feature 2 according to Level 2
FL4=MWI(P4); %Feature 4 according to Level 4
%--------------------------------------
%Level 1 [1]
d=sign(diff(FL1));
d1=[0;d];
d2=[d;0];
f1=find(d1==1);
f2=find(d2==-1);
EventsL1=intersect(f1,f2); %Detected events
%-------------------------------------------------------
%Level 2 [1]
meanL1=sum(FL2(EventsL1),1)/length(EventsL1);
NL=meanL1-meanL1*0.1;   %Start Noise Level
SL=meanL1+meanL1*0.1;   %Start Signal Level
threshold1=0.08;    %Threshold detection block 1
threshold2=0.7;     %Threshold detection block 2
[SignalL21,Noise1,DS1,Class1]=detectionblock(FL2,EventsL1,NL,SL,threshold1);
[SignalL22,Noise2,DS2,Class2]=detectionblock(FL2,EventsL1,NL,SL,threshold2);
%---------------------------------------------------
%Level 3 [1]
ClassL3=[]; 
for i=1:length(EventsL1)
    C1=Class1(i);
    C2=Class2(i);
    if C1==1
        if C2==1
            ClassL3=[ClassL3 1];   %Classification as Signal
        else
            delta1=(DS1(i)-threshold1)/(1-threshold1);
            delta2=(threshold2-DS2(i))/threshold2;
            if delta1>delta2
                ClassL3=[ClassL3 1]; %Classification as Signal
            else
                ClassL3=[ClassL3 0];  %Classification as Noise
            end
        end
    else
        if C2==1;
            ClassL3=[ClassL3 1]; %Classification as Signal
        else
            ClassL3=[ClassL3 0];  %Classification as Noise
        end
    end
end
SignalL3=EventsL1(find(ClassL3));   %Signal Level 3
NoiseL3=EventsL1(find(ClassL3==0)); %Noise Level 3
%--------------------------------------------
%Level 4 [1]
threshold=0.3;
VSL=(sum(FL4(SignalL3),1))/length(SignalL3); 
VNL=(sum(FL4(NoiseL3),1))/length(NoiseL3);   
SL=(sum(FL4(SignalL3),1))/length(SignalL3);   %Initial Signal Level
NL=(sum(FL4(NoiseL3),1))/length(NoiseL3);      %Initial Noise Level
SignalL4=[];
NoiseL4=[];
DsL4=[];   %Detection strength Level 4
for i=1:length(EventsL1)
    Pkt=EventsL1(i);    
    if ClassL3(i)==1;   %Classification after Level 3 as Signal
       SignalL4=[SignalL4,EventsL1(i)]; 
       SL=history(SL,FL4(Pkt));       
       Ds=(FL4(Pkt)-NL)/(SL-NL);       %Detection strength
       if Ds<0
           Ds=0;
       elseif Ds>1
           Ds=1;
       end
       DsL4=[DsL4 Ds];
    else        %Classification after Level 3 as Noise
       Ds=(FL4(Pkt)-NL)/(SL-NL);   
       if Ds<0
           Ds=0;
       elseif Ds>1
           Ds=1;
       end
       DsL4=[DsL4 Ds];
       if Ds>threshold          %new classification as Signal 
           SignalL4=[SignalL4,EventsL1(i)];  
           SL=history(SL,FL4(Pkt));     
       else                      %new classification as Noise
           NoiseL4=[NoiseL4,EventsL1(i)];
           NL=history(NL,FL4(Pkt));      
       end
   end
end
%------------------------------------------------
%Level 5  
%if the time between two RR complexes is too long => go back and check the
%events again with lower threshold
SignalL5=SignalL4;
NoiseL5=NoiseL4;
periods=diff(SignalL4);
M1=100;
a=1;
b=1/(M1)*ones(M1,1);
meanperiod=filter(b,a,periods); %mean of the RR intervals
SL=sum(FL4(SignalL4))/length(SignalL4);
NL=sum(FL4(NoiseL4))/length(NoiseL4);
threshold=0.2;
for i=1:length(periods)
    if periods(i)>meanperiod*1.5     %if RR-interval is to long
        intervall=SignalL4(i):SignalL4(i+1);
        critical=intersect(intervall,NoiseL4);   
        for j=1:length(critical)
            Ds=(FL4(critical(j))-NL)/(SL-NL); 
            if Ds>threshold         %Classification as Signal
                SignalL5=union(SignalL5,critical(j));   
                NoiseL5=setxor(NoiseL5,critical(j));
            end
        end
    end
end
%---------------------------------------------------
%Umrechnung auf Originalsignal (nicht downgesamplet)
Signaln=conversion(S,FL2,SignalL5,M,N,fs);
%----------------------------------------------------
%Level 6 If interval of two RR-complexes <0.24 => go back and delete one of them
height=FL2(SignalL5);   
Signal=Signaln;
temp=round(0.1*fs);
difference=diff(Signaln);  %Difference between two signal points
k=find(difference<temp);
for i=1:length(k)
    pkt1=SignalL5(k(i));
    pkt2=SignalL5(k(i)+1);
    verg=[height(k(i)),height(k(i)+1)]; 
    [x,j]=max(verg);    
    if j==1
        Signal=setxor(Signal,Signaln(k(i)+1)); %Deleting first Event
    else
        Signal=setxor(Signal,Signaln(k(i))); %Deleting second Event
    end
end
QRS=Signal;


%-------------------------------------------------------------------
%-------------------------------------------------------------------
%-------------------------------------------------------------------
%subfunctions

function y=MWI(S)

% MWI - Moving window integrator, computes the mean of two samples
%   y=MWI(S)
%
% INPUT
%   S       Signal
%
% OUTPUT
%   y       output signal
a=[0;S];
b=[S;0];
c=[a,b];
y=sum(c,2)/2;
y=y(1:length(y)-1);
%------------------------------------------------
function y=polyphase_imp(S,h,M)

% polyphase_imp - polyphase implementation of decimation filters [2]
%   y=polyphase_imp(S,h,M)
%
% INPUT
%   S       ecg signal data
%   h       filter coefficients
%   M       downsampling rate
%
% OUTPUT
%   y       filtered signal
%

%Determining polyphase components ek
e=cell(M,1);
l=1;
m=mod(length(h),M);
while m>0
    for n=1:ceil(length(h)/M)
        el(n)=h(M*(n-1)+l);
    end
    e{l}=el;
    l=l+1;
    m=m-1;
end
clear el;
for i=l:M
    for n=1:floor(length(h)/M)
        el(n)=h(M*(n-1)+i);      
    end
    e{i}=el;
end
%Filtering
max=ceil((length(S)+M)/M); 
Sdelay=S;
for i=1:M
    Sd=downsample(Sdelay,M);
    a=filter(e{i},1,Sd);     
    if length(a)<max
        a=[a;zeros(max-length(a),1)]; 
    end
    w(:,i)=a;
    Sdelay=[zeros(i,1);S];
end
y=sum(w,2);
%----------------------------------------------------------
function [Signal,Noise,VDs,Class]=detectionblock(mwi,Events,NL,SL,threshold)

% detectionblock - computation of one detection block 
%
%   [Signal,Noise,VDs,Class]=detectionblock(mwi,Events,NL,SL,threshold)
%
% INPUT
%   mwi         Output of the MWI
%   Events      Events of Level 1 (see [1])
%   NL          Initial Noise Level
%   SL          Initial Signal Level
%   threshold   Detection threshold (between [0,1])
%
% OUTPUT
%   Signal      Events which are computed as Signal
%   Noise       Events which are computed as Noise
%   VDs         Detection strength of the Events
%   Class       Classification: 0=noise, 1=signal

Signal=[];
Noise=[];
VDs=[];
Class=[];
sumsignal=SL;
sumnoise=NL;
for i=1:length(Events)
    P=Events(i);
    Ds=(mwi(P)-NL)/(SL-NL); %Detection strength
    if Ds<0
        Ds=0;
    elseif Ds>1
        Ds=1;
    end
    VDs=[VDs Ds];
    if Ds>threshold     %Classification as Signal
        Signal=[Signal P];
        Class=[Class;1];
        sumsignal=sumsignal+mwi(P);
        SL=sumsignal/(length(Signal)+1);    %Updating the Signal Level
    else        %Classification as Noise
        Noise=[Noise P];
        Class=[Class;0];
        sumnoise=sumnoise+mwi(P);
        NL=sumnoise/(length(Noise)+1);  %Updating the Noise Level
    end
end
%------------------------------------------------------------
function [pnew]=conversion(S,FL2,pold,M,N,fs)

% conversion - sets the fiducial points of the downsampled Signal on the
% samplepoints of the original Signal
% 
%   [pnew]=conversion(S,FL2,pold,M,N,fs)
%
% INPUT
%   S           Original ECG Signal
%   FL2         Feature of Level 2 [1]
%   pold        old fiducial points
%   M           M downsampling rate
%   N           filter order
%   fs          sample rate
%
% OUTPUT
%   pnew        new fiducial points
%

Signaln=pold;    
P=M;
Q=1;
FL2res=resample(FL2,P,Q);       %Resampling
nans1=isnan(S);
nans=find(nans1==1);
S(nans)=mean(S);    %Replaces NaNs in Signal
for i=1:length(Signaln)
    Signaln1(i)=Signaln(i)+(M-1)*(Signaln(i)-1);    
end
%------------------- Sets the fiducial points on the maximum of FL2
Signaln2=Signaln1;  
Signaln2=Signaln2';     
int=2*M;    %Window length for the new fiducial point
range=1:length(FL2res);
for i=1:length(Signaln2)
    start=Signaln2(i)-int/2;
    if start<1
        start=1;
    end
    stop=Signaln2(i)+int/2;
    if stop>length(FL2res)
        stop=length(FL2res);
    end
    intervall=start:stop;       %interval
    FL2int=FL2res(intervall);
    pkt=find(FL2int==max(FL2int));  %Setting point on maximum of FL2
    if length(pkt)==0   % if pkt=[];
        pkt=Signaln2(i)-start;
    else
        pkt=pkt(1); 
    end
    delay=N/2+M;
    Signaln3(i)=pkt+Signaln2(i)-int/2-delay;    %fiducial points according to FL2
end
%Sets the fiducial points on the maximum or minimum
%of the signal
Bw=5.6;   
Bwn=1/(fs/2)*Bw;
Wn=[Bwn 5*Bwn];
N1=32;
b=fir1(N1,Wn,'bandpass');
Sf=filtfilt(b,1,S);     %Filtered Signal with bandwidth 5.6-28 Hz
beg=round(1.5*M);
fin=1*M;
for i=1:length(Signaln3)
    start=Signaln3(i)-beg;
    if start<1
        start=1;
    end
    stop=Signaln3(i)+fin;
    if stop>length(Sf)
        stop=length(Sf);
    end
    intervall=start:stop;   %Window for the new fiducial point
    Sfint=abs(detrend(Sf(intervall),0));
    pkt=find(Sfint==max(Sfint));    %Setting point on maximum of Sfint
    if length(pkt)==0   %if pkt=[];
        pkt=Signaln3(i)-start;
    else
        pkt=pkt(1); 
    end
    pkt=pkt(1);
    Signaln4(i)=pkt+Signaln3(i)-beg-1;
end
Signal=Signaln4';   %New fiducial points according to the original signal

cutbeginning=find(Signal<N);    %Cutting out the first points because of initial transient of the filter in polyphase_imp
fpointsb=Signal(cutbeginning);
cutend=find(Signal>length(S)-N); %Cutting out the last points
fpointse=Signal(cutend);
pnew=setxor(Signal,[fpointsb;fpointse]);
%-------------------------------------------
function yn=history(ynm1,xn)

% history - computes y[n]=(1-lambda)*x[n]+lambda*y[n-1]
%
%   yn=history(ynm1,xn)

lambda=0.8; %forgetting factor
yn=(1-lambda)*xn+lambda*ynm1;

