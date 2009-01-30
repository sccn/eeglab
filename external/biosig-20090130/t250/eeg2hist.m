function [HDR]=eeg2hist(FILENAME,CHAN);
% EEG2HIST histogram analysis based on [1]
%  It displays the histograms of the recorded data, 
%  and allows editing of the thresholds for 
%  the automated overflow detection. 
% 
% [HDR]=EEG2HIST(FILENAME); 
%
% input: FILENAME   EEG-File
%        CHAN       Channel select
% output: 
%	 HDR	   header information
%        HDR.HIS   histograms for each channel
%	 HDR.RES   summary statistics based on the histogram analysis [1]
% 	 HDR.THRESHOLD  Threshold values for overflow detection 
%
% 
% see also: SOPEN, SLOAD, HIST2RES
%
% REFERENCES:
% [1] A. Schlögl, B. Kemp, T. Penzel, D. Kunz, S.-L. Himanen,A. Värri, G. Dorffner, G. Pfurtscheller.
%   Quality Control of polysomnographic Sleep Data by Histogram and EntropyAnalysis. 
%   Clin. Neurophysiol. 1999, Dec; 110(12): 2165 - 2170.
%   http://dx.doi.org/10.1016/S1388-2457(99)00172-8
%
% [2] A. Schlögl, G. Klösch, W. Koele, J. Zeitlhofer, P.Rappelsberger, G. Pfurtscheller
% Qualitätskontrolle von Biosignalen,
% Jahrestagung der Österreichischen Gesellschaft für Klinische Neurophysiologie, 27. Nov. 1999, Vienna.
%
% [3] http://www.dpmi.tu-graz.ac.at/~schloegl/lectures/Q/index.htm
%
% [4] A. Schlögl, Time Series Analysis toolbox for Matlab. 1996-2003
% http://www.dpmi.tu-graz.ac.at/~schloegl/matlab/tsa/

% 	$Id: eeg2hist.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (C) 2002,2003,2006,2007 by Alois Schloegl <a.schloegl@ieee.org>		
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
%    BioSig is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BioSig is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BioSig.  If not, see <http://www.gnu.org/licenses/>.


MODE=0; 
if nargin<2, CHAN=0; end; 
if ischar(CHAN), 
	MODE = 1; 
	CHAN = 0; 
end; 
 

[s,HDR]=sload(FILENAME,0,'OVERFLOWDETECTION','OFF','UCAL','ON');
%H = histo2(s);

if 1, 
HDR = sopen(FILENAME,'r',CHAN,'UCAL');	% open EEG file in uncalibrated mode (no scaling of the data)
if HDR.FILE.FID<0,
        fprintf(2,'EEG2HIST: couldnot open file %s.\n',FILENAME); 
        return;
end;

HDR.FLAG.UCAL = 1; % do not scale 
HDR.FLAG.OVERFLOWDETECTION = 0; % OFF
if CHAN<1, CHAN=1:HDR.NS; end;

H.datatype='HISTOGRAM';
if all(HDR.GDFTYP==3) % strcmp(HDR.TYPE,'BKR') | strcmp(HDR.TYPE,'CNT') | strcmp(HDR.TYPE,'EEG'),
        [s,HDR]=sread(HDR);
        
        H.H = zeros(2^16,HDR.NS);
        for l = 1:HDR.NS,
		if exist('OCTAVE_VERSION') > 2,
                	for k = double(s(:,l)')+2^15+1, H.H(k,l) = H.H(k,l)+1;  end;
		else
                	H.H(:,l)=sparse(double(s(:,l)')+2^15+1,1,1,2^16,1);
		end;
        end;
        H.X = [-2^15:2^15-1]'; 	%int16
        %tmp = find(any(H.H,2));
        %H.X = tmp-2^15-1; 	%int16
        %H.H = H.H(tmp,:);

elseif all(HDR.GDFTYP==4) 
        [s,HDR]=sread(HDR);
        
        H.H = zeros(2^16,HDR.NS);
        for l = 1:HDR.NS,
		if exist('OCTAVE_VERSION') > 2,
                	for k = s(:,l)'+1, H.H(k,l) = H.H(k,l)+1;  end;
		else
                	H.H(:,l)=sparse(s(:,l)'+1,1,1,2^16,1);
		end;
        end;
        H.X = [0:2^16-1]';
        %tmp = find(any(H.H,2));
        %H.X = tmp-1; 	%uint16
        %H.H = H.H(tmp,:);

        
elseif all(HDR.GDFTYP==(255+24)) %strcmp(HDR.TYPE,'BDF'),
        [s,HDR] = sread(HDR);
        
        H.H = sparse(2^24,HDR.NS);
        for l = 1:HDR.NS,
                H.H(:,l)=sparse(s(:,l)'+2^23+1,1,1,2^24,1);
                %for k = s(:,l)'+2^23+1, H.H(k,l) = H.H(k,l)+1; end;
        end;
        tmp = find(any(H.H,2));
        H.X = tmp-2^23-1; 
        H.H = H.H(tmp,:);
        
        
elseif (strcmp(HDR.TYPE,'EDF') | strcmp(HDR.TYPE,'GDF') |  strcmp(HDR.TYPE,'ACQ')) & all(HDR.GDFTYP(1)==HDR.GDFTYP) & (HDR.GDFTYP(1)<=3)
        NoBlks=ceil(60*HDR.SampleRate/HDR.SPR);

	if isfield(HDR.AS,'SPR')
	        bi=[0;cumsum(HDR.AS.SPR)];     
	else
		bi=0:HDR.NS;
		HDR.AS.spb = HDR.SPR;
        end; 
        CHAN = HDR.InChanSelect;
        ns=length(CHAN);
        
        H.H = zeros(2^16,ns);
        
        k=0;
        while (k<HDR.NRec) & ~feof(HDR.FILE.FID)
                % READ BLOCKS of DATA
                [S, count] = fread(HDR.FILE.FID,[HDR.AS.spb,NoBlks],gdfdatatype(HDR.GDFTYP(1)));
                if 0, count < HDR.AS.spb*NoBlks
                        fprintf(2,'   Warning EEG2HIST: read error, only %i of %i read\n', count, HDR.AS.spb*NoBlks);
                end;

                %%%%% HISTOGRAM
                for l=1:ns,
	                h = zeros(2^16,1);
			if exist('OCTAVE_VERSION','builtin'),
	                        for k=reshape(S(bi(CHAN(l))+1:bi(CHAN(l)+1),:),1,HDR.SPR(l)*NoBlks)+2^15+1, h(k,l) = h(k,l)+1; end;     
                        else
				h = sparse(S(bi(CHAN(l))+1:bi(CHAN(l)+1),:)+2^15+1,1,1,2^16,1);	                        
                        end;
			H.H(:,l) = H.H(:,l) + h;
                end;
                
                k=k+NoBlks; 
        end; % WHILE     
        tmp = find(any(H.H,2));
        H.X = [-2^15:2^15-1]'; 
        % (tmp-2^15-1); 	%int16
        %H.H = H.H(tmp,:);
        
elseif ~strcmp(HDR.TYPE,'unknown')
	[s,HDR]=sread(HDR); 
	H = histo2(s); 
	if any(HDR.GDFTYP>6)
	        fprintf(2,'WARNING EEG2HIST: data is not integer.\n');
	end; 	
		
else
        fprintf(2,'EEG2HIST: format %s not implemented yet.\n',HDR.TYPE);
end;
HDR = sclose(HDR);
end; 

HDR.HIS = H; 

%%% complete histogram and display it 
H.N = sumskipnan(H.H);
%H.X = [ones(size(H.X,1),1),repmat(H.X,1,length(CHAN))]*HDR.Calib; 	%int16

if strcmp(HDR.TYPE,'GDF')
	HDR.THRESHOLD = [HDR.DigMin(:),HDR.DigMax(:)];
end; 

if ~isfield(HDR,'THRESHOLD')
	HDR.THRESHOLD = 1e5*ones(HDR.NS,1)./HDR.Cal(:) *[-1,1];
end; 

        N=ceil(sqrt(size(H.H,2)));
        for K = 1:size(H.H,2);
                t = H.X(:,min(K,size(H.X,2)))*HDR.Calib(K+1,K)+HDR.Calib(1,K);
                h = H.H(:,K);

                mu = sumskipnan(t.*h)/H.N(K);
                x  = t-mu; 
                sd2= sumskipnan((x.^2).*h)./H.N(K);

                if 0, 
                elseif isfield(HDR,'THRESHOLD'),
                        MaxMin=HDR.THRESHOLD(K,[2,1])*HDR.Calib(K+1,K)+HDR.Calib(1,K);
		else                
			%[max(t) min(t)],
                        MaxMin=[max(t) min(t)];
                end;
                xrange = [min(t(h>0)),max(t(h>0))]; 
                xrange = xrange + [-1,1]*(diff(xrange)/2+eps);
 
                a(K)= subplot(ceil(size(H.H,2)/N),N,K);
                tmp = diff(t);
                dQ  = min(tmp(abs(tmp)>eps));
                tmp = sqrt(sum(h(h>0))/sqrt(2*pi*sd2)*dQ);

		h2 = h; 
		h2((t>min(MaxMin)) & (t<max(MaxMin)))=NaN; 

                semilogy(t,[h+.01,exp(-((t-mu).^2)/(sd2*2))/sqrt(2*pi*sd2)*sum(h(h>0))*dQ],'-',t,h2+.01,'r',mu+sqrt(sd2)*[-5 -3 -1 0 1 3 5]',tmp*ones(7,1),'+-',MaxMin,tmp,'rx');
                %v=axis; v=[v(1:2) 1 max(h)]; axis(v);
                v=axis; v=[xrange 1-eps max(h)]; axis(v);
                title(HDR.Label{K});
	end; 

if MODE,	
	figure(gcf); %set(gcf,'CurrentFigure');
	drawnow; 
	fprintf(1,'\nEEG2HIST:>Now You can modify the thresholds with mouse clicks.'); 
	fprintf(1,'\nEEG2HIST:>When you are finished, PRESS ANY KEY on the keyboard.'); 
	fprintf(1,'\nEEG2HIST:> (make sure the focus is on the figure window).'); 
	
	% modify Threshold 
	if ~isfield(HDR,'THRESHOLD'),
		HDR.THRESHOLD = repmat(NaN,HDR.NS,2); 
	end;
	b  = 0; 
	K0 = 0; 
	while (b<2),
		[x,y,b]=ginput(1);
		if isempty(b) break;end; 

		v=axis; 
		K = find(a==gca); 
                %min(K,size(H.X,2))
                t = H.X(:,min(K,size(H.X,2)))*HDR.Calib(K+1,K)+HDR.Calib(1,K);
                %HISTO=hist2pdf(HISTO);
                h = H.H(:,K);

		if K~=K0, MaxMin = [NaN,NaN]; end; 
		MaxMin = [MaxMin,x];
		tmp = h; 
		tmp((t>min(MaxMin)) & (t<max(MaxMin)))=NaN; 
                semilogy(t,[h+.01],'b-',t,tmp+.01,'-r');
                ix = sort(MaxMin);
		HDR.THRESHOLD(K,1:2) = (ix(1:2)-HDR.Calib(1,K))/HDR.Calib(K+1,K);
		K0 = K; 
                v=[v(1:2) 1-eps max(h)]; axis(v);
                title(HDR.Label{K});
	end; 		
end; 
fprintf(1,' <FINISHED>\n'); 


%%% complete return argument
HDR.HIS = H; 
%HDR.RES = hist2res(H); 
HDR.datatype = 'qc:histo'; 

