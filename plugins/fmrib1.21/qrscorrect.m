%qrscorrect() -  Correct for false positive and negative QRS peaks
%   and align peaks for maximum correlation.
%
% Usage:
%   >>Peaks=qrscorrect(Peaks,Ecg,fs)
%
% Inputs:
%   Peaks: Indeces of located QRS peaks
%   Ecg: vector of ecg data
%   sampling frequency.
%
% Ouput:
%   Peaks: Corrected peaks
%
% Note: progress bar designed to work as part of fmrib_qrsdetect.
%
% Author: Rami Niazy, FMRIB Centre, University of Oxford.
%
% 
%
% Copyright (c) 2004 Rami Niazy, FMRIB Centre, University of Oxford.
%


% Copyright (C) 2004 Rami K. Niazy, FMRIB Centre, University of Oxford
%                    rami@fmrib.ox.ac.uk
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%

%   Dec 23, 2004
%   use binary prcorr2
%   instead of corrcoef

%   Dec 17, 2004
%   use corrcoef instead of
%   corr2

%   Dec 15, 2004
%   Median of all peaks used
%   for correction instead of
%   section median.

%   Nov 4, 2004
%   Changed order of correction
%   Added additional alignment

%   Nov 3, 2004
%   Added check and 
%   adjustment of edge peaks clearance

%   Oct 26, 2004
%   Fixed bug for non-int fs case

function Peaks=qrscorrect(Peaks,Ecg,fs)

% init
%------
P=zeros(1,length(Ecg));
P(Peaks)=1;
dP_std=std(diff(Peaks));
dP_med=median(diff(Peaks));
pad=round(5*fs);
secL=round(10*fs);
sections=floor(length(P)/secL)+1;

if length(P((sections-1)*secL+1:end))<pad
    sections=sections-1;
end

barth=2;
barth_step=barth;
if isempty(javachk('desktop'))
    clearcom='clc';
else
    clearcom='home';
end

% Correct for FP
%---------------
for s=1:sections
    
    %wait bar
    if s==1
        barth=5;
        barth_step=barth;
        Flag25=0;
        Flag50=0;
        Flag75=0;
        fprintf('\nStage 2 of 5: Correcting for False Positive detection.\n');
       
    end
    %-----------------calc------------------------------------------
    
    
    
    f_flag=0;
    if s==1
        sec=P(1:secL+pad);
    elseif s==sections
        sec=P((s-1)*secL+1-pad:end);
    else
        sec=P((s-1)*secL+1-pad:s*secL+pad);
    end
    
    sec_P=find(sec==1);
    sec_dP=diff(sec_P);
%     dP_med=median(sec_dP);
   
    
    false_p=1;
    while ~isempty(false_p)
        false_p=find(sec_dP<(dP_med-3*dP_std));%-3*dP_std
        if ~isempty(false_p)
            f_flag=1;
            sec(sec_P(false_p(1)+1))=0;
            sec_P=find(sec==1);
            sec_dP=diff(sec_P);
            false_p=find(sec_dP<(dP_med-3*dP_std));
        end
    end
    
    if f_flag==1
        if s==1
            P(1:secL+pad)=sec;
        elseif s==sections
            P((s-1)*secL+1-pad:end)=sec;
        else
            P((s-1)*secL+1-pad:s*secL+pad)=sec;
        end
    end
    
    %---------------------------------------------------------------
    %update wait bar
    percentdone=floor(s*100/sections);
    if floor(percentdone)>=barth
        if percentdone>=25 & Flag25==0
            fprintf('25%% ')
            Flag25=1;
        elseif percentdone>=50 & Flag50==0
            fprintf('50%% ')
            Flag50=1;
        elseif percentdone>=75 & Flag75==0
            fprintf('75%% ')
            Flag75=1;
        elseif percentdone==100
            fprintf('100%%\n')
        else
            fprintf('.')
        end

        while barth<=percentdone
            barth=barth+barth_step;
        end
        if barth>100
            barth=100;
        end
    end 
end


%Check Edge Peaks
%----------------
Peaks=find(P==1);
pc1=2;
pc2=1;

while (pc1>1 | pc2>0)
	dP=diff(Peaks);
	hwinl=round(mean(dP)/2);
	searchw=round(0.33*mean(dP));
    
	pc1=1;
	if length(Peaks)>5
        for p=1:5
            try 
                tmp=Ecg(Peaks(p)-hwinl-searchw:Peaks(p)+hwinl);
                break;
            catch
                pc1=pc1+1;
            end
        end
	end
	
	if pc1>1
        Peaks(1:pc1-1)=[];
	end
	
	pc2=0;
	if length(Peaks)>5
        for p=length(Peaks):-1:length(Peaks-4)
            try 
                tmp=Ecg(Peaks(p)-hwinl:Peaks(p)+hwinl+searchw);
                break;
            catch
                pc2=pc2+1;
            end
        end
	end
	
	if pc2>0
        Peaks(end-pc2:end)=[];
	end
end

%Finding corrected peaks
%-----------------------
searchwL=length(searchw);
peaksL=length(Peaks);
qrs=zeros(hwinl*2+1,peaksL);
qrs2=zeros(hwinl*2+1,peaksL);
barth=2;

% Find average Heart Beat for use in alignment
%---------------------------------------------
for p=1:peaksL
    qrs(:,p)=Ecg(Peaks(p)-hwinl:Peaks(p)+hwinl);
end

m_qrs=mean(qrs,2);
% figure;plot(qrs)

% Align Peaks (1)
%----------------
for p=1:peaksL
    
    %wait bar
    if p==1
        barth=5;
        barth_step=barth;
        Flag25=0;
        Flag50=0;
        Flag75=0;
        fprintf('\nStage 3 of 5: Aligning QRS Peaks  (1)\n');
    end
    
    %-----------------calc------------------------------------------
    
    ppn=1;
    for B=Peaks(p)-searchw:Peaks(p)+searchw
        C(ppn)=prcorr2(Ecg(B-hwinl:B+hwinl),m_qrs);;
        ppn=ppn+1;
	end
	[CV,CP]=max(C);
	Beta=CP-(searchw+1);
	Peaks(p)=Peaks(p)+Beta;
    
    %---------------------------------------------------------------
    %update wait bar
    percentdone=floor(p*100/(peaksL));
    if floor(percentdone)>=barth
        if percentdone>=25 & Flag25==0
            fprintf('25%% ')
            Flag25=1;
        elseif percentdone>=50 & Flag50==0
            fprintf('50%% ')
            Flag50=1;
        elseif percentdone>=75 & Flag75==0
            fprintf('75%% ')
            Flag75=1;
        elseif percentdone==100
            fprintf('100%%\n')
        else
            fprintf('.')
        end

        while barth<=percentdone
            barth=barth+barth_step;
        end
        if barth>100
            barth=100;
        end
    end    
end
barth=2;
P=zeros(1,length(Ecg));
P(Peaks)=1;

% Correct for FN
%---------------
for s=1:sections
    
    %wait bar
    if s==1
        barth=5;
        barth_step=barth;
        Flag25=0;
        Flag50=0;
        Flag75=0;
        fprintf('\nStage 4 of 5: Correcting for False Negative detection.\n');
    end
    
    %-----------------calc------------------------------------------
    
    
    
    f_flag=0;
    if s==1
        sec=P(1:secL+pad);
    elseif s==sections
        sec=P((s-1)*secL+1-pad:end);
    else
        sec=P((s-1)*secL+1-pad:s*secL+pad);
    end
    
    sec_P=find(sec==1);
    sec_dP=diff(sec_P);
    dP_med=median(sec_dP);
       
    false_n=1;
    while ~isempty(false_n)
        false_n=find(sec_dP>(1.5*dP_med));
        if ~isempty(false_n)
            f_flag=1;
            sec(sec_P(false_n(1))+round(dP_med))=1;%sec_dP(false_n(1)-1)
            sec_P=find(sec==1);
            sec_dP=diff(sec_P);
            false_n=find(sec_dP>(1.5*dP_med));
        end
    end
    
    
    if f_flag==1
        if s==1
            P(1:secL+pad)=sec;
        elseif s==sections
            P((s-1)*secL+1-pad:end)=sec;
        else
            P((s-1)*secL+1-pad:s*secL+pad)=sec;
        end
    end
    
    %---------------------------------------------------------------
    %update wait bar
    percentdone=floor(s*100/sections);
    if floor(percentdone)>=barth
        if percentdone>=25 & Flag25==0
            fprintf('25%% ')
            Flag25=1;
        elseif percentdone>=50 & Flag50==0
            fprintf('50%% ')
            Flag50=1;
        elseif percentdone>=75 & Flag75==0
            fprintf('75%% ')
            Flag75=1;
        elseif percentdone==100
            fprintf('100%%\n')
        else
            fprintf('.')
        end

        while barth<=percentdone
            barth=barth+barth_step;
        end
        if barth>100
            barth=100;
        end
    end
end

barth=2;
Peaks=find(P==1);
peaksL=length(Peaks);
barth=2;
% Align Peaks (2)
%----------------
for p=2:peaksL-1
    
    %wait bar
    if p==2
        barth=5;
        barth_step=barth;
        Flag25=0;
        Flag50=0;
        Flag75=0;
         fprintf('\nStage 5 of 5: Aligning QRS Peaks  (2)\n');
    end
  
    %-----------------calc------------------------------------------
    
    ppn=1;
    for B=Peaks(p)-searchw:Peaks(p)+searchw
        C(ppn)=prcorr2(Ecg(B-hwinl:B+hwinl),m_qrs);
        ppn=ppn+1;
	end
	[CV,CP]=max(C);
	Beta=CP-(searchw+1);
	Peaks(p)=Peaks(p)+Beta;
    
    %---------------------------------------------------------------
    %update wait bar
    percentdone=floor((p-1)*100/(peaksL-2));
    if floor(percentdone)>=barth
        if percentdone>=25 & Flag25==0
            fprintf('25%% ')
            Flag25=1;
        elseif percentdone>=50 & Flag50==0
            fprintf('50%% ')
            Flag50=1;
        elseif percentdone>=75 & Flag75==0
            fprintf('75%% ')
            Flag75=1;
        elseif percentdone==100
            fprintf('100%%\n')
        else
            fprintf('.')
        end

        while barth<=percentdone
            barth=barth+barth_step;
        end
        if barth>100
            barth=100;
        end
    end    
end

% for p=2:peaksL-1
%     qrs2(:,p-1)=Ecg(Peaks(p)-hwinl:Peaks(p)+hwinl);
% end
% figure;plot(qrs2);
return;
