function [arg0,s0,ODI,x,EDF]=desatur(arg1,arg2,arg3)
% calculates desaturation of the Oxygen in the blood
% [V,s,ODI,x]=desatur(FN [,CH]) 
% [V,s,ODI,x]=desatur(s,Fs)
%
% FN    filename (of a file in EDF format) 
% CH    channel number, default 16.
% s 	recorded signal,
% Fs    sampling rate
%
% Automatic evaluation of oxygen saturation can deliver:
% V(1) - mean oxygen saturation
% V(2) - time spent with oxygen saturation below 90% in minutes
% V(3) - time spent with oxygen saturation below 80% in minutes
% V(4) - time spent with oxygen saturation below 70% in minutes
% V(5) - number of oxygen desaturations
% V(6) - mean value of oxygen desaturation in percent
% V(7) - mean duration of oxygen desaturation in seconds
%
% V(8) - median oxygen saturation
% V(9) - median value of oxygen desaturation in percent
% V(10)- median duration of oxygen desaturation in seconds
% 
% s   original Sa02 signal
% ODI Oxygen de-saturation index (number of detected de-saturations per hour)
% x   detector output 
%
%  see also: SDF-tools
%
%  Version 1.12
%  Copyright 2000 by A. Schloegl <a.schloegl@ieee.org>
%  Department for Medical Informatics, University of Technology Graz, AUSTRIA
%
% References:
%  [1] simple and robust method
%  [2] simple and robust method, detection of invalid samples
%  [3] personal correspondance with Georg Gruber from University hospital Vienna. 
%  [4] Taha BH, Dempsey JA, Weber SM, Badr MS, Skatrud JB, Young TB, Jacques AJ, Seow KC.
%      Automated detection and classification of sleep-disordered breathing from conventional polysomnography data.
%      Sleep. 1997 Nov;20(11):991-1001.% 
%  [5] J-C. Vazquez, W.H. Tsai, W.W. Flemons, A. Masuda, R. Brant, E.Hajduk, W.A. Whitelaw, J.E. Remmers, 
%      Automated analysis of digital oximetry in the diagnosis of obstractive sleep apnoea, 
%      Thorax 2000; 55:302-307.

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



%global k Fs x
        
if nargin>2
        Mode=arg3;
else
        Mode=1;
end;


%%%%% Load data 
if isstr(arg1),  % if string, assume arg1 is a EDF-Filename and SaO2 is at channel 16
        if nargin<2, CH=16; else CH=arg2; end;
        %EDF=sdfopen(arg1,'r',0);CH=find(kanaltyp(EDF)=='S');
        [pf,FN,ext,v]=fileparts(arg1);
        if exist(['/home/home8/siesta/SaO2/' FN '.int16'],'file')
                load(['/home/home8/siesta/t300/' FN 'res.mat'],'EDF');
                EDF=sedfchk(EDF);
                fid=fopen(['/home/home8/siesta/SaO2/' FN '.int16'],'r');
                [sa,c]=fread(fid,inf,'int16');
                fclose(fid);
                s=sa*EDF.Cal(16)+EDF.Off(16);
                EDF.PhysDim=EDF.PhysDim(CH,:);
        else
                EDF=sdfopen(arg1,'r',CH,'UCAL');
                [sa,EDF] = sdfread(EDF,inf);
                EDF=sdfclose(EDF);
                
                fid=fopen(['/home/home8/siesta/SaO2/' FN '.int16'],'w+');
                c=fwrite(fid,sa,'int16');
                fclose(fid);
                s=sa*EDF.Cal(16)+EDF.Off(16);
                %save(['/home/home8/siesta/t300/' FN 'res.mat'],'EDF','-append');
        end;
        Fs = EDF.SPR(CH)/EDF.Dur;   
        if isempty([findstr(EDF.PhysDim,'%'),findstr(lower(EDF.PhysDim),'percent')]),     %check Physical Dimension  
                if ~isfield(EDF.FILE,'stderr'), EDF.FILE.stderr=2; end; % define error output
                fprintf(EDF.FILE.stderr,'Warning DESATUR: Invalid Physical dimension\n');
        end;
        EDF.AS.SaO2.Fs=Fs;
        EDF.AS.SaO2.Chan=CH;
        
else              % assume arg1 is the SaO2 signal with given sampling rate
        s = arg1;   
        if isstruct(arg2); % if is struct, assume EDF struct
                Fs=EDF.AS.SaO2.Fs;
        else               % arg2 is the sampling rate
                Fs =arg2;
        end;
end;

s0=s; % output argument returns original data

if size(s,1)==1; s=s'); end; % column vector

% subsampling to 4Hz
if Fs==100,
        t=25;        
elseif Fs==256,
        t=64;
else
        t=1;
end;
s=rs(s,ones(t,1)/t);
Fs=Fs/t;

N1 = inf; %30*Fs; % window length of baseline filter, inf means baseline is the median of the all-nigth recording
N2 = 10*Fs; % window length of thresholding
TH = 4;     % Threshold in %
LowerBound = 51; % defines lower bound of valid range
UpperBound = 100; % defines upper bound of valid range

if Mode>1, % correct invalid values, corrently no correction is performed. 
        invalid=((s<LowerBound) | (s>UpperBound))'; % identify indeces with invalids;
        tmp=find(~invalid); 
        if ~isempty(tmp),
                tmp=tmp(1); % find first valid value
                for k=find(invalid(tmp+1:length(invalid)))+tmp;  s(k)=s(k-1); end; 
                for k=tmp-find(invalid(tmp-1:-1:1)); s(k)=s(k+1); end;
        end;
end;

m0=median(s);   
if ~isinf(N1),
	m = medfilt1(s,N1); %median filter and
	m = filter([zeros(1,N1/2),1],1,m); % delay
else
        m=m0;
end;

switch Mode
case {1,2}
        x = filter(ones(1,N2),1,(m-TH)>s)>=N2; %thresholding, smoothing, thresholding 
        x = [x(N2/2+1:length(x)); zeros(N2/2,1)]; % correction of the delay
case 3 
        x = filter(ones(1,N2),1,(m-medfilt1(s,N2))>TH)>=N2; %thresholding, smoothing, thresholding 
        x = [x(N2/2+1:length(x)); zeros(N2/2,1)]; % correction of the delay
case 4  % Taha (1997), 
	if (Fs>2) & (Fs==fix(Fs/2)*2), % downsampling to 2 Hz, if possible 
	    s=rs(s,ones(2,1)/2);
	    Fs=Fs/2;
	else
	    fprintf(1,'DESATUR: sampling rate is not 2Hz but %i\n',Fs);    
	end;
	
	x=repmat(0,size(s));
	s(length(s)+1)=-inf;
	ak0=nan; k=1;
	while k<=length(s)-Fs,
		while isnan(ak0)	% while a-point not found
			%(s(k)-s(k+Fs)),
			if (s(k)-s(k+Fs)) > 0.1
				a=(s(k)+s(k+Fs))/2; 	 % keep a-value in s(k0)	
				ak0=k+Fs/2;
			else k=k+1;
			end 
                end;
                bk0=ak0;
	        while (k<length(s)-1) & (s(k) > s(k+[0 2])),
			k=k+1;
			bk0=k;
                        %tmp=s(k+[-1 1]);
		end;
                if s(bk0) > a-2;   % not lower than 2%
		        bk0=nan; 		% not a valid b-point
			ak0=nan;
		end;  
		Q=1;   
		while ~isnan(ak0) & Q      
			k = k+1;
			Q = s(k) > max(a-1,s(bk0)+3);
		end;
		ck0=k;
		if abs((ck0-ak0-1)/Fs-35)<=25, % 10<=ck0-ak0<=60
			x(ak0:ck0-1)=1;		
		end;
		k=k+1;
	end;	    	
case 5
        if (Fs>1) % downsampling to 1 Hz, 
	    s=rs(s,ones(Fs,1)/Fs);
	    Fs=1;
	else
	    fprintf(1,'DESATUR: sampling rate is not 1Hz but %i\n',Fs);    
	end;
        
	x=repmat(0,size(s));
        x2=x;
        
        x1=filter([3 -1 -1 -1],1,sign(diff(s)))>5;
        
        if 1 
                x1(1:50)=0; %[]; x1(1:10),
                for k=find(diff(x1)>0)',
                        x(k)= any(s((-3:-1)+k) < mean(s(k+(-49:0))>prctile(s(k+(-49:0)),95))-4);
                end;
        else
                x=x1;
        end;
end;        
        
ODI = mean(diff(x)>0)/Fs*3600; % OD's per hour [1/h]

if ~any(x)
        arg0=[mean(s), sum(s<90)/Fs/60, sum(s<80)/Fs/60, sum(s<70)/Fs/60, sum(diff(x)>0), NaN, 0, m0, NaN, 0];
else
        arg0=[mean(s), sum(s<90)/Fs/60, sum(s<80)/Fs/60, sum(s<70)/Fs/60, sum(diff(x)>0), mean(s(find(x>0))), mean(find(diff(x)<0)-find(diff(x)>0))/Fs, m0, median(s(find(x>0))), median(find(diff(x)<0)-find(diff(x)>0))/Fs];
end;



function [y1]=rs(y1,T,f2)
% [y2] = rs(y1, T) resamples y1 to the target sampling rate y2 using T

        [yr,yc]=size(y1);
        LEN=yr/f1;
        for k=0:LEN-1
                y1(k*f2+(1:f2),:)=T'*y1(k*f1+(1:f1),:);
        end;
        y1=y1(1:LEN*f2,:);
end;