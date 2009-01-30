function [QRStime_corr,dr,dr_corr,U,ANNOT] = qrscorr(QRSindex,Fs,Nmax);
% Identification and correction of detection errors on QRS-beat sequences
%
% This algorithm identifies anomalies like FP, FN or ectopic beats.
%  The FP and FN are corrected by deleting or inserting one or multiple beats.
%
% INPUT:
%  QRSindex:     Sample values of the detected QRS-Complexes(minimum (8+(Nmax-1)Samples)
%  Fs:           Sample Frequency
%  Nmax:         Maximum number of insertions, deletions, movements at one position (default = 20)
%
% OUTPUT:
%  QRStime_corr: Time values of the corrected QRS-Complexes
%  dr:           Time values of the derivative of the original instantaneous
%                 heart rate:
%     dr = 2 * | (t(k-1) - 2t(k) + t(k+1)) / ((t(k-1) - t(k))*(t(k-1) - t(k+1))*(t(k) - t(k+1))) |
%  dr_corr:      Time values of the derivative of the corrected instantaneous
%                 heart rate:
%  U:            Calculated threshold for dr (Umax = 0.5 [s-2])
%  ANNOT:        Structure that annotates the corrected and uncorrected beats
%  ANNOT:         ANNOT.ins: Time values of all inserted beats
%                 ANNOT.del: Time values of all deleted beats
%                 ANNOT.mov: Time values of all moved beats
%                 ANNOT.NOR: Time values of all uncorrected beats
%
%
% [QRStime_corr,ANNOT] = qrscorr(QRSindex,Fs,Nmax);
%
% Example:
%  [QRStime_corr,ANNOT] = qrscorr(QRSindex,256,15);
%
%
% Reference:
% [1] J. Mateo, P. Laguna, Analysis of Heart Rate Variability in Presence
%      of Ectopic Beats Using the Heart Timing Signal
%     IEEE Transactions on biomedical engineering,
%      Vol.50, No.3, March 2003

% Filename: qrscorr.m
% This file is part of the BIOSIG-toolbox http://biosig.sf.net/
% Last modified: 2003/09/16
% Copyright (c) 2003 by Johannes Peham
%	$Revision: 1.1 $
%	$Id: qrscorr.m,v 1.1 2009-01-30 06:04:44 arno Exp $
%
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


%Input arguments
   if nargin==2, Nmax=20;end; %Nmax default = 20

%number of detected QRS-Complexes must be at least:
%2(preload) + 3(afterload) + (Nmax-1)(maximal no. of corrected beats at one position) + 3(real data)
   if length(QRSindex)<(2+3+(Nmax-1)+3), error('input QRS-indices must be at least (8 + (Nmax-1))!'); end;
   if (Nmax<1)|(length(Nmax)<1), error('Nmax must be at least 1');end;

%Converting into row vector
   QRSindex=QRSindex(:)';
   Nmax=Nmax(1);

%Initialisation of the indices
   t=QRSindex./Fs; %Converting Samples to time
   dr=drcalc(t);

%Calculating the mean RR-intervall
   RRmean=mean(diff(t));

% Calculating the threshold U
   U = min(std(dr)*4.3, 0.5);  %Upper limit of U = 0.5


%**************************************************************%
% Classification of the anomaly
%**************************************************************%
%Initialisation
   Tcorr=t;
   dr_corr=dr;
   i=3; %starting at t(3)
   ANNOT.del=[];
   ANNOT.ins=[];
   ANNOT.mov=[];
   ANNOT.NOR=[];

   while i<=(length(Tcorr)-3-(Nmax-1)),

      %N is counting the multiple corrections
      N=0;%initialising
  
      if dr_corr(i) >= U %Derivative smaller then the treshold?
         [Tcorr,ANNOT,N]=classification(Tcorr,ANNOT,i,U,RRmean,Nmax);
         dr_corr = drcalc(Tcorr);
      else
         Tcorr(1:2)=NaN; Tcorr(end-2-(Nmax-1):end)=NaN;
         del=(length(find(ANNOT.del == Tcorr(i))));%was the beat deleted?
         ins=(length(find(ANNOT.ins == Tcorr(i))));%was the beat inserted?
         mov=(length(find(ANNOT.mov == Tcorr(i))));%was the beat moved?
         if (del==0)&(ins==0)&(mov==0),
            ANNOT.NOR(end+1)=Tcorr(i);%annotate normal beat only if the beat was NOT corrected in any way
         end;
      end;
      i=i+1+N; %jump over multiple anomal beats
   end;

%Calculation of the derivative of the heart rate
   dr_corr=drcalc(Tcorr);

%Cutting the data
   dr=dr(2:end-1);
   Tcorr=Tcorr(3:end-3-(Nmax-1));
   dr_corr=dr_corr(4:end-4-(Nmax-1));

%Output variables
   QRStime_corr = Tcorr;

return;  
%-----------end main---------------------------------------------%


%-----------subroutines------------------------------------------%

function [Tcorr,ANNOT,N]=classification(Tcorr,ANNOT,k,U,RRmean,Nmax)
%Classifies the kth beat of the QRS series in ANNOT and corrects the time
% series Tcorr
%
%[Tcorr,ANNOT]=classification(Tcorr,ANNOT,k,U,Nmax)
%
% Tcorr: time values of the corrected QRSindices
% ANNOT: structure that annotates the corrected and uncorrected beats
% ANNOT: ANNOT.ins: Time values of all inserted beats
%        ANNOT.del: Time values of all deleted beats
%        ANNOT.mov: Time values of all moved beats
%        ANNOT.NOR: Time values of all uncorrected beats
% N: number of additional corrected beats (if one was corrected N=0)
% k: index of the beat that should be corrected
% U: threshold for the derivative of the heart rate
% Nmax: maximal number of insertions, deletions, movements

%Initialisation
T=Tcorr; %making a copy of Tcorr
if ((Nmax>(length(Tcorr)-8))|(Nmax<=0)), error('Nmax = maximal del,ins,mov exceeds data length or equals zero');end;
N=0; %number of additional insertions, deletions


while (N < Nmax)

    
%3) Inserting intermediate between tk-1 tk
    Tcorr((k+N):length(Tcorr)+1+N)=Tcorr((k-1):length(Tcorr));
    for i=1:N+1
        Tcorr(k-1+i) = Tcorr(k-1) + i*(Tcorr(k+1+N)-Tcorr(k-1))/(2+N);
    end;
    dr_corr=drcalc(Tcorr,k-1,k+N);

    if dr_corr < U,
        ANNOT.ins(end+1:end+1+N)=Tcorr(k:k+N);%annotate insertion
       %Correction is only valid for Tcorr(3) to Tcorr(end-3)
        Tcorr(1:2)=NaN; Tcorr(end-2-(Nmax-1):end)=NaN;
        return;
    else Tcorr=T;
    end;

%4) Inserting intermediate between tk tk+1
    Tcorr((k+1+N):length(Tcorr)+1+N)=Tcorr(k:length(Tcorr));
    for i=1:N+1
        Tcorr(k+i) = Tcorr(k) + i*(Tcorr(k+2+N)-Tcorr(k))/(2+N);
    end;
    dr_corr=drcalc(Tcorr,k,k+1+N);

    if dr_corr < U, 
        ANNOT.ins(end+1:end+1+N)=Tcorr(k+1:k+1+N);%annotate insertion
        ANNOT.NOR(end+1)=Tcorr(k);%annotate previous beat as normal
       %Correction is only valid for Tcorr(3) to Tcorr(end-3)
        Tcorr(1:2)=NaN; Tcorr(end-2-(Nmax-1):end)=NaN;
        return;
    else Tcorr=T;
    end;

    
%1) Deleting tk
    if length(Tcorr) <= 8, return;end; %minimum length of the indices is 8!
    Tcorr((k):(length(Tcorr)-1-N))=Tcorr((k+1+N):length(Tcorr-N));
        
    dr_corr=drcalc(Tcorr,k-1,k);
    
    if (dr_corr < U)&((Tcorr(k) - Tcorr(k-1))<2*RRmean), 
        Tcorr((end-N):end)=[];
        ANNOT.del(end+1:end+1+N)=T(k:k+N);%annotate deletion
        ANNOT.NOR(end+1)=Tcorr(k);%annotate current position as normal
        %Correction is only valid for Tcorr(3) to Tcorr(end-3)
        Tcorr(1:2)=NaN; Tcorr(end-2-(Nmax-1):end)=NaN;
        N=0;
        return;
    else Tcorr=T;
    end;

%2) Deleting tk+1
    if length(Tcorr) <= 8, return;end; %minimum length of the indices is 8!
    Tcorr((k+1):(length(Tcorr)-1-N))=Tcorr((k+2+N):length(Tcorr-N));
        
    dr_corr=drcalc(Tcorr,k,k+1);
    
    if (dr_corr < U)&((Tcorr(k+1) - Tcorr(k))<2*RRmean), 
        Tcorr((end-N):end)=[];
        ANNOT.del(end+1:end+1+N)=T(k+1:k+1+N);%annotate deletion
        ANNOT.NOR(end+1)=T(k);%annotate previous beat as normal
       %Correction is only valid for Tcorr(3) to Tcorr(end-3)
        Tcorr(1:2)=NaN; Tcorr(end-2-(Nmax-1):end)=NaN;
        N=0;
        return;
    else Tcorr=T;
    end;

    
%6) Moving tk+1 in the middle of tk tk+2

    for i=1:N+1
        Tcorr(k+i) = Tcorr(k) + i*((Tcorr(k+2+N) - Tcorr(k))/(2+N));
    end;

    dr_corr=drcalc(Tcorr,k,k+N);

    if dr_corr < U, 
        Tcorr=T; %if movement yield to an improvement the position k
                 % was an ectopic beat and has to be specially corrected (-->ECTBcorr)
        ANNOT.mov(end+1:end+1+N)=Tcorr(k+1:k+1+N);%annotate ectopic beat
        ANNOT.NOR(end+1)=Tcorr(k);%annotate previous beat as normal
        N=N+1;%next position is after the last moved k+(N+1)
        %Correction is only valid for Tcorr(3) to Tcorr(end-3)
        Tcorr(1:2)=NaN; Tcorr(end-2-(Nmax-1):end)=NaN;
        return;%annotate Insertion
    else Tcorr=T;
    end;

%5) Moving tk in the middle of tk-1 tk+1

    for i=1:N+1
        Tcorr(k+i-1) = Tcorr(k-1) + i*((Tcorr(k+1+N) - Tcorr(k-1))/(2+N));
    end;

    dr_corr=drcalc(Tcorr,k-1,k-1+N);

    if dr_corr < U, 
        Tcorr=T; %if movement yield to an improvement the position k or k:k+N
                 % was an ectopic beat and has to be specially corrected (-->ECTBcorr)
        if (length(find(ANNOT.mov == Tcorr(k))) == 0) %anootate only, if it was not already annotated in 6)!
            ANNOT.mov(end+1:end+1+N)=Tcorr(k:k+N);%annotate ectopic beat
        end;
        %Correction is only valid for Tcorr(3) to Tcorr(end-3)
        Tcorr(1:2)=NaN; Tcorr(end-2-(Nmax-1):end)=NaN;
        return;
    else Tcorr=T;
    end;



N=N+1;
end;

N=0;%no multiple correction
ANNOT.NOR(end+1)=Tcorr(k);
%Correction is only valid for Tcorr(3) to Tcorr(end-3)
Tcorr(1:2)=NaN; Tcorr(end-2-(Nmax-1):end)=NaN;
return

%**************************************************************%
%Calculation of the derivative of the instantaneous heart rate
%**************************************************************%
function dr = drcalc(t,j,jo);
%dr is the the derivative of the instantaneous heart rate referring to [1]
%dr = drcalc(t) calculation of all positions
%               result=[NaN dr(2) dr(3) ... dr(N-1) NaN]
%dr = drcalc(t,k) calculation of position k result is 1x1
%dr = drcalc(t,k1,k2) calculation from position k1 to k2 with the length
%                     from k1 to k2

if nargin == 3
   if ((j>jo)) error('k1>k2');
   elseif ((jo>=length(t))|(j<=1)) error('index k1 k2 exceeds dimensions'); 
   end;
   for i=j:jo
       if ( (t(i-1)==t(i)) | (t(i-1)==t(i+1)) | (t(i)==t(i+1))), dr(i-j+1)=0;
       else
           dr(i-j+1)=2*abs((t(i-1)-2*t(i)+t(i+1))/((t(i-1)-t(i))*(t(i-1)-t(i+1))*(t(i)-t(i+1))));
       end;
   end;
elseif nargin == 2
   if ((j>=length(t))|(j<=1)) error('index k exceeds dimensions');end;
    dr=2*abs((t(j-1)-2*t(j)+t(j+1))/((t(j-1)-t(j))*(t(j-1)-t(j+1))*(t(j)-t(j+1))));
elseif nargin ==1
   for j=2:(length(t)-1)
      dr(j)=2*abs((t(j-1)-2*t(j)+t(j+1))/((t(j-1)-t(j))*(t(j-1)-t(j+1))*(t(j)-t(j+1))));
   end;
   dr(1)=NaN;
   dr(length(t))=NaN;
else error('too few/many input parameters');
end;
return;
