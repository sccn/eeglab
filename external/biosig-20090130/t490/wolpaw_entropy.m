function B=wolpaw_entropy(P,N)
%  WOLPAWENTROPY calculates the Mutual Information according to [1]
%   
%   B = wolpaw_entropy(acc,M)
%
%   B = log2(M) + acc*log2(acc) + (1-acc)*log2((1-acc)/(M-1))
%
%   acc Accuracy [0..1]  1 = 100%; 
%   N   number of classes
%   B	mutual information 	 	 
% 	
%
% Reference(s): 
% [1] J.R. Wolpaw, N. Birbaumer, W.J. Heetderks, D.J. McFarland, P. Huntereckham, 
% 	G. Schalk, E. Donchin, L.A. Quatrano, C-J. Robinson and T.M. Vaughan.
% 	Brain-Computer Interface Technology: A Review of the First Inernational Meeting. 
% 	IEEE Transactions on Rehabilitation Engineering 8(2) (Jun 2000) 164-173.


%    $Id: wolpaw_entropy.m,v 1.1 2009-01-30 06:04:50 arno Exp $
%    Copyright (C) 2003-2005 by Alois Schloegl <a.schloegl@ieee.org>	
%    This is part of the BIOSIG-toolbox http://biosig.sf.net/


%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


B=log2(N)+P.*log2(P)+(1-P).*log2((1-P)./(N-1));

ix=find(P==0);
B(ix)=log2(N)+(1-P(ix)).*log2((1-P(ix))./(N-1));

ix=find(P==1);
B(ix)=log2(N)+P(ix).*log2(P(ix));

