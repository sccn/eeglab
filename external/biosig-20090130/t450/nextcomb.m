function [a,y]=nextcomb(y,N,k);
%
%The function[a,y] = nextcomb(y,N,k) is a subfunction in the main-function fdp.m.
%For help use the main-function.
%
%
%Copyright (C) 2006 by Claudia Hemmelmann <claudia.hemmelmann@mti.uni-jena.de>
%Institute of Medical Statistics, Computer Sciences and Documantation
%University of Jena
%This work was supported by DFG Project VO 683/2-1
%This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
%
%%***
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
%
%--------------------------------------------------------------------------


[u,j]=max(([y(2:k),N+1]-y(1:k))>ones(1,k));
%j=suche(1);

%suche=y+ones(1,k);
%for j=1:k
%i=find(y==suche(j));
%if isempty(i),break;end;
%end;
%j
a=0;
if u==0
a=1;return;
end;

y(j)=y(j)+1;
y(1:j-1)=[1:j-1];

