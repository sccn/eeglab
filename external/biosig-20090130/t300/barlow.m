function [A,F,S]=barlow(S,UC,A)
%  BARLOW calculates ACTIVITY, MOBILITY, COMPLEXITY
%           
%  [Amplitude, Frequency, SPI] = barlow(...)
%  
%  [...] = barlow(S,0)
%       calculates stationary Hjorth parameter 
%  [...] = barlow(S,UC) with 0<UC<1,
%       calculates time-varying Hjorth parameter using 
%       exponential window 
%  [...] = barlow(S,N) with N>1,
%       calculates time-varying Hjorth parameter using 
%       rectangulare window of length N
%  [...] = barlow(S,B,A) with B>=1 oder length(B)>1,
%       calulates time-varying Hjorth parameters using 
%       transfer function B(z)/A(z) for windowing 
%
%       S       data (each channel is a column)
%       UC      update coefficient 
%       B,A     filter coefficients (window function) 
%              
% see also: TDP, HJORTH, WACKERMANN
%
% REFERENCE(S):
% [1] Goncharova II, Barlow JS.
%   Changes in EEG mean frequency and spectral purity during spontaneous alpha blocking.
%   Electroencephalogr Clin Neurophysiol. 1990 Sep;76(3):197-204. 

%	$Id: barlow.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (C) 2005 by Alois Schloegl <a.schloegl@ieee.org>
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

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

[N,K] = size(S); 	% number of electrodes K, number of samples N

d0 = S;
d1 = diff([zeros(1,K);S ],[],1);
d2 = diff([zeros(1,K);d1],[],1);

FLAG_ReplaceNaN = 0;

if nargin<2, 
        UC = 0; 
end;
if nargin<3;
        if UC==0,
                                
        elseif UC>=1,
                B = ones(1,UC);
                A = UC;
        elseif UC<1,
                FLAG_ReplaceNaN = 1;
                B = UC; 
                A = [1, UC-1];
        end;
else
        B = UC;    
end;


if ~UC,
        m0 = mean(abs(d0));
        m1 = mean(abs(d1));
        m2 = mean(abs(d2));
else
        if FLAG_ReplaceNaN;
                d0(isnan(d0)) = 0;
                d1(isnan(d1)) = 0;
                d2(isnan(d2)) = 0;
        end;
        m0 = filter(B,A,abs(d0))./filter(B,A,double(~isnan(d0)));
        m1 = filter(B,A,abs(d1))./filter(B,A,double(~isnan(d1)));
        m2 = filter(B,A,abs(d2))./filter(B,A,double(~isnan(d2)));
end;

A = m0; 
F = m1./m0; 
S = (F.*F)./(m2.*A);

