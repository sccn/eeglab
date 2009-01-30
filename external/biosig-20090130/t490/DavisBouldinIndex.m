function [DBI]=DavisBouldinIndex(d,c,kk);
% Davis-Bouldin-Index is a Cluster separation index (CSI)

% KAPPA.M estimates Cohen's kappa coefficient 
%
% [kap,sd,H,z,OA,SA,MI] = kappa(d1,d2);
% [kap,sd,H,z,OA,SA,MI] = kappa(H);
%
% d1    data of scorer 1 
% d2    data of scorer 2 
%
% kap	Cohen's kappa coefficient point
% se	standard error of the kappa estimate
% H	data scheme (Concordance matrix or confusion matrix)
% z	z-score
% OA	overall agreement 
% SA	specific agreement 
% MI 	Mutual information or transfer information (in [bits])
%
% Reference(s):
%   [1] Bezdek, J.C.; Pal, N.R.; Some new indexes of cluster validity
% Systems, Man and Cybernetics, Part B, IEEE Transactions on Volume 28, Issue 3, June 1998 Page(s):301 - 315 
%
%  

% http://www.ucl.ac.uk/oncology/MicroCore/HTML_resource/distances_popup.htm

%	$Id: DavisBouldinIndex.m,v 1.1 2009-01-30 06:04:50 arno Exp $
%	Copyright (c) 2006 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% version 2 of the License, or (at your option) any later version.
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
%


CL = unique(c); 
M = length(CL);

t = 2; q = 2; 
D = zeros(M,M); 

for k=1:M,
        %v(k,:)=mean(dk,1);    % center of each cluster
        [dk, v(k,:)] = center(d(c==CL(k),:),1); 
        S(k) = sum(sqrt(sum(dk.^2,2)).^q).^(1/q);
end;         

for k1=1:M,
for k2=1:M,
        d(k1,k2)= sum(abs(v(k1,:)-v(k2,:)).^t).^(1/t);  % Minkowski Distance of order t   
end;
end; 

[x,y] = meshgrid(S);
DBI = mean(max((x+y)./D + diag(repmat(NaN,1,M))));


