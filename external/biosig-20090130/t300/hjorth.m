function [ACTIVITY, MOBILITY, COMPLEXITY,m0,m1,m2] = hjorth(S,UC,A)
%  HJORTH calculates ACTIVITY, MOBILITY, COMPLEXITY
%           
%  [ACTIVITY, MOBILITY, COMPLEXITY] = hjorth(...)
%  
%  [...] = hjorth(S,0)
%       calculates stationary Hjorth parameter 
%  [...] = hjorth(S,UC) with 0<UC<1,
%       calculates time-varying Hjorth parameter using 
%       exponential window 
%  [...] = hjorth(S,N) with N>1,
%       calculates time-varying Hjorth parameter using 
%       rectangulare window of length N
%  [...] = hjorth(S,B,A) with B>=1 oder length(B)>1,
%       calulates time-varying Hjorth parameters using 
%       transfer function B(z)/A(z) for windowing 
%
%       S       data (each channel is a column)
%       UC      update coefficient 
%       B,A     filter coefficients (window function) 
%              
% see also: TDP, BARLOW, WACKERMANN
%
% REFERENCE(S):
% [1] B. Hjorth, 
%   EEG analysis based on time domain properties
%   Electroencephalography and Clinical Neurophysiology, vol. 29, no. 3, pp. 306–310, September 1970.
% [2] B. Hjorth, 
%   Time Domain Descriptors and their Relation to particulare Model for Generation of EEG activity. 
%   in G. Dolce, H. Kunkel: CEAN Computerized EEG Analysis, Gustav Fischer 1975, S.3-8. 


%	$Id: hjorth.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (C) 2004,2008 by Alois Schloegl <a.schloegl@ieee.org>
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

%m0 = mean(sumsq(S,2));
d0 = S;
%m1 = mean(sumsq(diff(S,[],1),2));
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
        m0 = mean(d0.^2);
        m1 = mean(d1.^2);
        m2 = mean(d2.^2);
else
        if FLAG_ReplaceNaN;
                d0(isnan(d0)) = 0;
                d1(isnan(d1)) = 0;
                d2(isnan(d2)) = 0;
        end;
        m0 = filter(B,A,d0.^2)./filter(B,A,double(~isnan(d0)));
        m1 = filter(B,A,d1.^2)./filter(B,A,double(~isnan(d1)));
        m2 = filter(B,A,d2.^2)./filter(B,A,double(~isnan(d2)));
end;

ACTIVITY   = m0;
MOBILITY   = sqrt(m1./m0); 
COMPLEXITY = sqrt(m2./m1); 
