function [SIGMA,PHI,OMEGA,m0,m1,m2] = wackermann(S,UC,A)
% Wackermann calculates SIGMA, PHI and OMEGA. 
%   
%  [SIGMA,PHI,OMEGA] = wackermann(...)
%  
%  [...] = wackermann(S,0)
%       calculates stationary Wackermann parameter 
%  [...] = wackermann(S,UC) with 0<UC<1,
%       calculates time-varying Wackermann parameter using 
%       exponential window 
%  [...] = wackermann(S,N) with N>1,
%       calculates time-varying Wackermann parameter using 
%       rectangulare window of length N
%  [...] = wackermann(S,B,A) with B>=1 oder length(B)>1,
%       calulates time-varying Wackermann parameters using 
%       transfer function B(z)/A(z) for windowing 
%
%       S       data (each channel is a column)
%       UC      update coefficient 
%       B,A     filter coefficients (window function) 
%
% see also: TDP, BARLOW, HJORTH
%
% REFERENCE(S):
% [1] Jiri Wackermann, Towards a quantitative characterization of
% functional states of the brain: from the non-linear methodology to the
% global linear descriptor. International Journal of Psychophysiology, 34 (1999) 65-80.
%
%

%	$Id: wackermann.m,v 1.1 2009-01-30 06:04:44 arno Exp $
%	Copyright (C) 2004,2005 by Alois Schloegl <a.schloegl@ieee.org>
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

if nargin<2, 
        UC = 0; 
end;
FLAG_ReplaceNaN = 0;
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

d0 = sum(S.*S,2);
d1 = sum(diff([zeros(1,K);S],[],1).^2,2);

[k1,k2] = meshgrid(1:size(S,2));
k1 = k1(:); k2 = k2(:); 
d2 = S(:,k1).*S(:,k2);  

if UC==0,
        m0 = mean(d0);
        m1 = mean(d1);
        m2 = mean(d2);
elseif 0, UC>=1
        m0 = rs(d0,UC,1);
        m1 = rs(d1,UC,1);
        m2 = rs(d2,UC,1);
else
        if FLAG_ReplaceNaN;
                d0(isnan(d0)) = 0;
                d1(isnan(d1)) = 0;
                d2(isnan(d2)) = 0;
        end;
        %m0 = mean(sumsq(S,2));
        m0 = filter(B,A,d0)./filter(B,A,real(~isnan(d0)));
        %m1 = mean(sumsq(diff(S,[],1),2));
        m1 = filter(B,A,d1)./filter(B,A,real(~isnan(d1)));
        m2 = filter(B,A,d2)./filter(B,A,real(~isnan(d2)));
end;

SIGMA = sqrt(m0/K);

PHI   = sqrt(m1./m0)/(2*pi);

OMEGA = repmat(NaN,size(m0));
for k = 1:size(m2,1),
        if all(finite(m2(k,:))), 
                L = eig(reshape(m2(k,:), [K,K]));
                L = L./sum(L);
                if all(L>0)
                        OMEGA(k) = -sum(L.*log(L));
                end;
        end;
end;        

