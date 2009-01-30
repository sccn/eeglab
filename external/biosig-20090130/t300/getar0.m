function [a0,A0] = getar0(S,P,NTR,NN);
% GETAR0 calculates average AR parameters for initialization of 
% AAR estimation
%
% C0 = getar0(S,P,NTR,NN);
%    S  signal 
%    P  list of model orders
%    NTR  number of realizations
%    NN   length each segment
%    C0 extended covarance matrix (contains a0 and A0)
%
% [a0,A0] = getar0(...); % in future this will become obsolete 
%    a0 average AR-parameters	
%    A0 covariance 

%	$Revision: 1.1 $
%	$Id: getar0.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (c) 1996-2003 by Alois Schloegl
%	e-mail: a.schloegl@ieee.org	


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


S  = S(:);
tmp= floor(rand(1,NTR)*(length(S)-NN));
z0 = zeros(NTR,NN);
for k = 1:length(tmp),
        z0(k,:) = S(tmp(k) + (1:NN))';
end;
[MX,pe] = durlev(acovf(z0,max(P)));

for k=1:length(P),
        a = MX(:,P(k)*(P(k)-1)/2+(1:P(k)));
	if nargout<2,
	        a0{k} = covm(a,'E');
        else	% OBSOLETE 
		a0{k} = mean(a);
    		A0{k} = covm(a,'D');
	end;
end;	

