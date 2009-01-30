function [Y,RRI] = berger(RRI,Fs)
% Resampling with the Berger algorithm
% [HRV,RRI] = berger(RRI, Fs)
% [HRV,RRI] = berger(ONSET, Fs)
% [HRV,RRI] = berger(HDR, Fs)
% 
% RRI 	R-to-R interval 
% ONSET onset time QRS-complex
% Fs	target sampling rate
% HDR	header struct as defined by SOPEN, SLOAD. 
% 	HDR.EVENT must contain the QRS events
% HRV 	heart rate variability sampled with Fs
% RRI	R-to-R interval sampled with Fs
%
% Reference(s):
% [1] Berger RD, Akselrod S, Gordon D, Cohen RJ. 
%    An efficient algorithm for spectral analysis of heart rate variability.
%    IEEE Trans Biomed Eng. 1986 Sep;33(9):900-4.

%       $Id: berger.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (c) 1997-2005, 2007 by Alois Schloegl <a.schloegl@ieee.org>	
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

if isstruct(RRI)
        Fs0 = NaN; 
        if isfield(RRI,'EVENT')
                EVENT = RRI.EVENT; 
        else
                EVENT = RRI; 
        end;
        if ~isfield(EVENT,'SampleRate')
                warning('Invalid input argument causes unknown source sampleing rate.');
        end;
        if isfield(EVENT,'POS') && isfield(EVENT,'TYP') && isfield(EVENT,'CHN') && isfield(EVENT,'DUR');
                ix = find(EVENT.TYP==hex2dec('0501'));
                if all(EVENT.CHN(ix(1)) == EVENT.CHN(ix));
                        on = EVENT.POS(EVENT.TYP==hex2dec('0501'))/EVENT.SampleRate;
                end;
        elseif isfield(EVENT,'POS') && isfield(EVENT,'TYP');
                on = EVENT.POS(EVENT.TYP==hex2dec('0501'))/EVENT.SampleRate;
        end;
        
elseif ~any(diff(RRI)<0),	
        on = RRI(:);%-min(RRI);
else		% calculate onset of bursts 
        on = cumsum([0;RRI(:)]);
end;

if ~exist('on','var')
        error('Invalid input argument'); 
end;        
%T=1/Fs;

N  = ceil(max(on)*Fs);
Y  = zeros(N,1);

on = [on; NaN] * Fs;

IX = 1;
for n = 1:N,
	while (on(IX+1)<(n-1)), IX=IX+1; end;
        ix = IX+1;
        CONT = 1;
        while CONT,
               
		% 1st case: RR interval is totally encompassed by sample window
	        if 0, 
		elseif n<=on(1)
                        Y(n) = nan;		
		elseif (on(ix-1) >= (n-1)) & (on(ix) <= (n+1)),
                        Y(n) = Y(n) + 1;
		% 2nd case: RR interval partially overlaps beginning of sample window
	        elseif (on(ix-1) <  (n-1)) & (on(ix) <  (n+1)),
                        Y(n) = Y(n) + (on(ix)-(n-1))/(on(ix)-on(ix-1));
		% 3rd case: RR interval partially overlaps the end of sample window
	        elseif (on(ix-1) >  (n-1)) & (on(ix) >  (n+1)),
                        Y(n) = Y(n) + ((n+1) - on(ix-1))/(on(ix)-on(ix-1));
		% 4th case: RR interval totally overlaps the sample window
	        elseif (on(ix-1) <= (n-1)) & (on(ix) >= (n+1)),
                        Y(n) = Y(n) + 2/(on(ix)-on(ix-1));
		elseif isnan(on(ix)),
                        Y(n) = nan;
                else
                        disp([ix,on(ix-1),on(ix), n, (n+1),Y(n)])
                        fprintf(2,'Warning BERGER: Invalid program state');
                        Y(n) = nan;
                end;
                if (on(ix)<(n+1)),%(on(ix-1)<(n+1)*T), 
                        ix = ix+1;
                else
                        CONT=0;
                end;
        end;
end;        

RRI = (2*max(on)/N/Fs)./Y;
%Y=Fs/2*N/max(on(1:length(on)-1))*Y;
Y = 60./RRI;


