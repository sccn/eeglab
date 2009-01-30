function [perm]=perm(x1,x2)
% PERM gives a vector containing all possible sums of two vectors
% [perm]=perm(x1,x2)
%
% It is useful for identifying segments of trigger data. 


%	$Revision: 1.1 $
%	$Id: perm.m,v 1.1 2009-01-30 06:04:49 arno Exp $
%	Copyright (c) 1998-2003 by Alois Schloegl
%	a.schloegl@ieee.org	


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


[r1,c1]=size(x1);
[r2,c2]=size(x2);

perm=ones(r2*c2,1)*reshape(x1,1,r1*c1)+reshape(x2,r2*c2,1)*ones(1,r1*c1);