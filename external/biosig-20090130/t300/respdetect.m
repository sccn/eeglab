function r=respdetect(S,fs);

% respdetect - detection of respiration
%
%   r=respdetect(S,fs);
%
% INPUT
%   S       signal data
%   fs      sample rate
%
% OUTPUT
%   r       fiducial points of the respiration

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



S=full(S);
nans1=isnan(S);
nans=find(nans1==1);
S(nans)=mean(S);    %Replaces NaNs

Wls=120;     %window length (seconds)
wl=Wls*fs;   %window length (sample points)
%----------------------------------filtering
[n,wn,bet,ty]=kaiserord([1 2],[1 0],[0.1 0.1],fs);
b=fir1(n,wn,ty,kaiser(n+1,bet));
Sf=filtfilt(b,1,S);
Sf=detrend(Sf,0);
%------------------------------------
begin=Sf(1:wl);     %start window of length Wls
beginplus=Sf(1:wl+wl);  
Sigsq=begin.^2;
q=0.5;
qu=quantile(Sigsq,q);
pkt1t=gettrigger(Sigsq,qu);
he=find(begin(pkt1t)>0);
pkt1t=pkt1t(he);    %only the positive ones

begindiff=diff(beginplus);
nulldiff=gettrigger(-begindiff,0); %local maxima
maxm=[];    
for j=1:length(pkt1t);   
    akt=pkt1t(j);
    gri=find(nulldiff>akt);
    gr=nulldiff(gri);
    maxm=[maxm,gr(1)];  %local maxima of the window (above threshold qu)
end
%----------------------------------------------
%rest of the signal, after the initial window
Sigsq=Sf;
Sigdiff=diff(Sf);
Sattelp=gettrigger(-Sigdiff,0);
Sattelp=[Sattelp;length(Sf)]; 
vecmax=maxm;    %maxima
ende=1;
while ende==1
    vecmaxn=length(vecmax);
    maxpkt=vecmax(vecmaxn);  %last point of vecmax
    numb=15;    %number of elements over one averages
    if vecmaxn<=numb
        numb=vecmaxn-1;
    end
    med=mean(Sigsq(vecmax(vecmaxn-numb:vecmaxn))); %average of the elements
    startw=maxpkt;  %beginning of the window
    endw=maxpkt+wl; %end of the window
    if endw>length(Sigsq)   
        endw=length(Sigsq);
    end
    Sigp=Sigsq(startw:endw);    %window
    pkt1t=gettrigger(Sigp,med*q);
    if length(pkt1t)>=1;
        pkttriggneu=pkt1t(1)+maxpkt;
        gri=find(Sattelp>=pkttriggneu);
        pktneu=Sattelp(gri(1));     %first point of gri is the new point of vecmax
        vecmax=[vecmax,pktneu];     %local maxima of the signal 
        ende=1;
    else
        ende=0; %end while
    end
end
Sfdet=detrend(Sf,0);
pktn=gettrigger(-Sfdet,0);    %zero-crossings
pktn=[pktn;length(Sf)];     
null=[];    
for j=1:length(vecmax);   %first zero-crossing after vecmax
    akt=vecmax(j);
    gri=find(pktn>=akt);
    gr=pktn(gri);
    null=[null,gr(1)];
end
r=unique(null);


