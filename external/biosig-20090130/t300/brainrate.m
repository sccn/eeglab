function [br,sef90,sef95,br2] = brainrate(s,Fs,UC,A)
% BRAINRATE estimates the weighted mean frequency according to [1] 
%  Other (similar) parameters are the spectral edge frequency or
%  the Hjorth's Mobility parameter or Barlow's center frequency.   
%           
%  [BRAINRATE, SEF90, SEF95] = brainrate(...)
%  
%  [...] = brainrate(S,Fs)
%  [...] = brainrate(S,Fs,0)
%       calculate stationary brainrate parameter 
%  [...] = brainrate(S,Fs,UC) with 0<UC<1,
%       calculates time-varying brainrate parameter using 
%       exponential window 
% 
% Input: 
%	S	data (each channel is a column)
%	UC	update coefficient (0<UC<1) 
%
% Output: 
%	BRAINRATE	weighted mean frequency [1] 
% 	SEF90	spectral edge frequency (power below SEF90 is 90% of total power)        
% 	SEF95	spectral edge frequency (power below SEF95 is 95% of total power)        
%
% see also: HJORTH, BARLOW 
%
% REFERENCE(S):
% [1] Nada Pop-Jordanova and Jordan Pop-Jordanov
%     Spectrum-weighted EEG frequency ("Brainrate") as a quantitative 
%     indicator of arousal 
%     Contributions, Sec. Biol. Med. Sci., MASA, XXVI, 2, p. 35–42 (2005)
%     ISSN 0351–3254, UDK: 616.831-073.97

%	$Id: brainrate.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (C) 2006,2007 by Alois Schloegl <a.schloegl@ieee.org>
%	This is part of the BIOSIG-toolbox http://biosig.sf.net/

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


[N,K] = size(s); 	% number of electrodes K, number of samples N
MOP = 15; 		% order of autoregressive model 

if nargin<3, 
	UC = 0; 
	B = []; 
elseif nargin==3
	if UC==0,
		br = size(1,K); 
		sef90 = size(1,K); 
		sef95 = size(1,K); 
	elseif (UC>0) & (UC<1),
		B = UC; 
		A = [1, UC-1];
		br = size(N,K); 
		sef90 = size(N,K); 
		sef95 = size(N,K);
	else
		fprintf(2,'ERROR BrainRate: update coefficient UC(%f) is outside of interval ]0,1[\n',UC); 
		br =[]; sef90=[]; sef95=[];
		return; 	 
	end;
end; 

s = center(s,1);
for k2 = 1:K,
	if ~UC,
		[a,r,pe] = lattice(s(:,k2)',MOP); 	
	else
		X = tvaar(s(:,k2),MOP,UC); 
		X = tvaar(s(:,k2),X);
		a = X.AAR;  
		pe= X.PE; 
	end;
	
	fi = [2:4:18];	% f_i list of center frequencies
	for k1 = 1:size(a,1),
		%h = freqz(sqrt(pe(k1)/Fs),[1,-a(k1,:)],f,Fs); 
		% scaling not important because relative spectral distribution is used
		h = freqz(1,[1,-a(k1,:)],fi,Fs); 
		h2= abs(h); 
		br(k1,k2) = sum(fi.*h2)/sum(h2); 
	end; 

	if nargout>1,
		for k1 = 1:size(a,1),
			f = (0:256)/512*Fs; 
			h = freqz(1,[1,-a(k1,:)],f,Fs);
			h2= cumsum(real(h.*conj(h)),2);
			sef90(k1,k2) = f(min(find(h2>=.90*h2(end)))); 
			sef95(k1,k2) = f(min(find(h2>=.95*h2(end)))); 
		end; 
	end; 
end;
