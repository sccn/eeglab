function [Y,RRI]=berger(RRI,Fs)
% Resampling with the Berger algorithm
% [HRV,RRI] = berger(RRI, Fs)
% [HRV,RRI] = berger(ONSET, Fs)
% 
% RRI 	R-to-R interval 
% ONSET onset time QRS-complex
% Fs	target sampling rate
% HRV 	heart rate variability sampled with Fs
% RRI	R-to-R interval sampled with Fs
%

%	Version 0.25
%	12 Feb 2001
%	Copyright (c) 1997-2001 by  Alois Schloegl
%	a.schloegl@ieee.org	

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


if ~any(diff(RRI)<0),	
        on = RRI(:)-min(RRI);
else		% calculate onset of bursts 
        on = cumsum([0;RRI(:)]);
end;

%T=1/Fs;

N=ceil(max(on)*Fs);
Y=zeros(N,1);

on = [on; nan] * Fs;

IX = 1;
%on,break,
for n=1:N,
	while (on(IX+1)<(n-1)), IX=IX+1; end;
        ix=IX+1;
        CONT=1;
        while CONT,
                
		% 1st case: RR interval is totally encompassed by sample window
	        if     ((on(ix-1) >= (n-1)) & (on(ix) <= (n+1))),
                        Y(n) = Y(n) + 1;
		% 2nd case: RR interval partially overlaps beginning of sample window
	        elseif ((on(ix-1) <  (n-1)) & (on(ix) < (n+1))),
                        Y(n) = Y(n) + (on(ix)-(n-1))/(on(ix)-on(ix-1));
		% 3rd case: RR interval partially overlaps the end of sample window
	        elseif ((on(ix-1) > (n-1)) & (on(ix) >  (n+1))),
                        Y(n) = Y(n) + ((n+1) - on(ix-1))/(on(ix)-on(ix-1));
		% 4th case: RR interval totally overlaps the sample window
	        elseif ((on(ix-1) <=  (n-1)) & (on(ix) >=  (n+1))),
                        Y(n) = Y(n) + 2/(on(ix)-on(ix-1));
                else
                        disp([ix,on(ix-1)/Fs,on(ix)/Fs, n, (n+1),Y(n)])                       
                        fprintf(2,'Warning %s: Invalid program state',upper(mfilename));
                        Y(n) = nan;
                end;
                if (on(ix)<(n+1)),%(on(ix-1)<(n+1)*T), 
                        ix = ix+1;
                else
                        CONT=0;
                end;
        end;
end;        

RRI=(2*max(on)/N/Fs)./Y;
%Y=Fs/2*N/max(on(1:length(on)-1))*Y;
Y=60./RRI;


