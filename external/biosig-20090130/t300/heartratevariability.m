function [X] = heartratevariability(RRI,arg2)
% HeartRateVariability Analysis according to [1]
%
% X = heartratevariability(RRI [,units])
% X = heartratevariability(qrsindex [,Fs])
% X = heartratevariability(HDR)
% X = heartratevariability(HDR.EVENT)
% X = heartratevariability(filename)
% 
% 
% INPUT
%   RRI 	R-R-intervales [in seconds]
%   HDR		as defined in the header structure of BioSig 
%		and returned by QRS-detection
%   filename 	with event information
%   Fs 		sampling rate - used for conversion into time axis
%   units	time units e.g. 'ms' (default: 's')  
%
% OUTPUT
%   X  		struct containing the results as defined by [1]
%   X.meanRR      	meanRR = meanNN
%   X.SDRR		standard deviaation of RR intervales
%   X.RMSSD       	rmsSD = SDSD
%     NN50count1
%     NN50count2
%     NN50count
%	pNN50
%	SD1		width of Poincaré plot; equivalent to sqrt(2)*RMSSD [2]
%	SD2		length of Poincaré plot; i.e. 2SDRR²+SDSD²/2 [2]
%	r_RR 		correlation coefficient [2]
%   X.VLF               power of very low frequency band (< 0.04 Hz)
%   X.LF                power of low frequency band (0.04-0.15 Hz)
%   X.HF                power of high low frequency band (0.15-0.4 Hz)
%   X.TotalPower        power of high low frequency band (0.15-0.4 Hz)
%   X.LFHFratio         LF/HF-ratio
%   X.LFnu              normalized units of LF power (0.04-0.15 Hz)
%   X.HFnu              normalized units of HF power  (0.15-0.4 Hz)
%
%  semilogy(X.f,X.ASpectrum) shows the spectral density function  
%
% The spectral estimates are based on an autoregressive spectrum estimator 
% of the data which is oversampled by a factor of 4 using the Berger method.  
% The default model order is 15. In order to change these default settings, 
% change in the source code line 211 (oversampling factor) and/or 
% line 230 (order of the autoregressive model); 
%
% see also: QRSDETECT, BERGER, EVENTCODES.TXT
%
% Reference(s):
% [1] Heart Rate Variability
%       Standards of Measurement, physilogcial interpretation and clinical use.  
%       Taskforce of the European Society for Cardiology and the North Americal Society of Pacing and Electrophysiology.         
%       European Heart Journal (1996) 17, 354-381. 
% [2] M. Brennan, M.Palaniswami, P. Kamen
%	Do Existing Measures of Poincaré Plot Geometriy Reflect Nonlinear Features of Heart Rate Variablilty?
%	IEEE Trans Biomedical Eng. 48(11),2001, 
% [3] U. Rajendra Acharya, K. Paul Joseph,N. Kannathal, Choo Min Lim, Jasjit S. Suri. 
%	Heart rate variability: a review.
%	Med Bio Eng Comput (2006) 44:1031–1051

%	$Id: heartratevariability.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (C) 2005,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

 
% TODO: 
% TINN triangular interpolation of NN intervals 
% CD Correlation dimension
% DFA detrended fluctuation analysis 
% LLA Largest Lyaponow exponent
% ApEn Aproximate Entropy
% FD	Fractal Dimension 
% H	Hurst Exponent
% SampEn Sample Entropy 


%%%%%%%% check and convert the input arguments  %%%%%%%%%%
if ischar(RRI)
if exist(RRI,'file'); 
	RRI = sload(RRI); 
end;
end; 

Fs = 1; 
X.PhysDim = 's'; 
if nargin>1,
	if ischar(arg2)
		X.PhysDim = arg2;
	elseif isnumeric(arg2)
		Fs = arg2;	
	end; 
end;

	
if isstruct(RRI)
        Fs0 = NaN; 
        if isfield(RRI,'EVENT')
        	HDR = RRI; 
                EVENT = HDR.EVENT; 
		if isfield(HDR,'T0')
			X.T0 = HDR.T0; 
		end;	
		if isfield(HDR,'Patient')
			X.Patient = HDR.Patient;
		end;	
        else
                EVENT = RRI; 
        end;

        X.PhysDim = 's'; 
        if isfield(EVENT,'SampleRate')
                Fs0 = EVENT.SampleRate; 
        elseif isfield(RRI.SampleRate)
                Fs0 = RRI.EVENT.SampleRate; 
        else    
                warning('Invalid input argument causes unknown source sampleing rate.');
                Fs0 = 1; 
                X.PhysDim = '?'; 
        end;
        if isfield(EVENT,'POS') & isfield(EVENT,'TYP') & isfield(EVENT,'CHN') & isfield(EVENT,'DUR');
                ix = find(EVENT.TYP==hex2dec('0501'));
                if all(EVENT.CHN(ix(1)) == EVENT.CHN(ix));
                        on = EVENT.POS(EVENT.TYP==hex2dec('0501'))/Fs0;
                end;
        elseif isfield(EVENT,'POS') & isfield(EVENT,'TYP');
                on = EVENT.POS(EVENT.TYP==hex2dec('0501'))/Fs0;
        end;
        NN = diff(on);         
        
elseif ~any(diff(RRI)<0),	
        on = RRI(:)/Fs; 
        NN = diff(on);
%	X.PhysDim = 's'; 
else	
        NN = RRI(:); 
        on = [0;cumsum(NN)];
end;

if nargout<1,
	% scatterplot for testing of the distribution 
	plot(NN(1:end-1),NN(2:end),'x');drawnow;
	return;
end;

%%%%%%%%  convert from any time unit into ms %%%%%%%%%%%%%%%%
[PhysDimCode,scale1] = physicalunits(X.PhysDim); 
X.PhysDim = 'ms';
[X.PhysDimCode,scale2] = physicalunits(X.PhysDim); 
t_scale = scale2./scale1;
NN = NN/t_scale; 
on = on/t_scale;


X.datatype = 'HeartRateVariability'; 

%%%%%%%%% time domain parameters %%%%%%%%%%%%%%%%%%%%
X.meanNN   = mean(NN); 
X.SDNN     = std(NN); 
X.RMSSD    = rms(diff(NN));  % SDSD
X.NN50count1 = sum(-diff(NN)>0.050/t_scale);
X.NN50count2 = sum( diff(NN)>0.050/t_scale);
X.NN50count  = sum(abs(diff(NN))>0.050/t_scale);
X.pNN50    = X.NN50count/sum(~isnan(NN)); 


g = acovf(center(NN(:)'),2);
X.SD1 	= sqrt(g(1)-g(2));
X.SD2 	= sqrt(g(1)+g(2));
X.SD1 	= sqrt(g(1)-g(2));
X.r_RR	= g(2)/g(1);

t   = round(NN * 128)/128; 
HIS = histo3(t(:)); 
X.HRVindex128 = HIS.N/max(HIS.H); 

%still missing 
%X.SDANN = 0;
%X.SDNNindex = 0; 

%SDNNindex
%SDSD  = RMSSD


%%%%%%% AR-based spectrum  analysis %%%%%%%%%%%%%
OS = 1; 
if 0, 
	y = log(NN); 
	[y,m]=center(y); 
	f0= exp(-m)/t_scale
elseif 0, 
	y = NN; 
	[y,m] = center(y); 
	f0= 1/(m*t_scale);
else
	%% factor 1000 because data is converted to [ms], berger expects [s]
	OS = 4; 
	f0  = OS*1000/(X.meanNN);%% four-times oversampling
        [hrv,y] = berger(on/1000,f0); % resampleing 
	[y,m] = center(y*1000); 
end;

pmax = 100;
% choose levinson-durbin or Burg algorithm 
[mx,pe]=durlev(acovf(y(:)',pmax));
%[mx,pe]=lattice(y(:)',pmax); 

n = sum(~isnan(y));
[FPE,AIC,BIC,SBC,MDL,CAT,PHI,optFPE,optAIC,optBIC,optSBC,optMDL,optCAT,optPHI]=selmo(pe/pe(1),n);
X.mops = [optFPE,optAIC,optBIC,optSBC,optMDL,optCAT,optPHI]; 


% select model order - vary the model order in order to check how robust the results are with respect to the model order  
X.mop = optBIC;
%X.mop = optAIC; 
%X.mop = 15;
[a,r] = arcext(mx,X.mop);

[h,f] = freqz(sqrt(pe(X.mop+1)/f0),[1,-a],[0:512]/512*f0/OS,f0);
X.ASpectrum = abs(h);
X.f = f; 
ix = f<0.04;
try
find(ix),f0,OS,
X.VLF = trapz(f(ix),abs(h(ix)).^2);
ix = (f>0.04) & (f<0.15);
X.LF = trapz(f(ix),abs(h(ix)).^2);
ix = (f>0.15) & (f<0.40);
X.HF = trapz(f(ix),abs(h(ix)).^2);
ix = (f>0.15) & (f<0.40);
X.TotalPower = trapz(f,abs(h).^2);
X.LFHFratio = X.LF./X.HF; 
X.LFnu = X.LF./(X.TotalPower-X.VLF);
X.HFnu = X.HF./(X.TotalPower-X.VLF);
end;

%%%%%%% FFT-based spectrum  analysis %%%%%%%%%%%%%
% currently not implemented 


%%%%%%% slope %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if diff(on([1,end])*t_scale)>5*3600, 	%calculate only if more than 5 hours 
	t = (on(1:end-1) + on(2:end))/2;
	t = t * t_scale;
	f1 = 0.0001; f2 = 0.04; 
	[x1,n1] = sumskipnan(NN.*exp(-2*pi*j*f1*t));
	[x2,n2] = sumskipnan(NN.*exp(-2*pi*j*f2*t));
	X.alpha = 2*log(abs((x1*n2)/(x2*n1)))./log(f1/f2);
end; 
