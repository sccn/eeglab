function t=ttest3(Input,n1)
%
%The function t=ttest3(Input,n1) is a subfunction in the main-function fdp.m.
%For help use the main-function.
%
%Copyright (C) 2006 by Claudia Hemmelmann <claudia.hemmelmann@mti.uni-jena.de>
%Institute of Medical Statistics, Computer Sciences and Documantation
%University of Jena
%This work was supported by DFG Project VO 683/2-1
%This is part of the BIOSIG-toolbox http://biosig.sf.net/
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




[n,k]=size(Input);
n2=n-n1;
X1=Input(1:n1,:);
X2=Input(n1+1:n,:);
z1=sum(X1,1)/n1;
z2=sum(X2,1)/n2;
%eigener Test
s1=X1-ones(n1,1)*sum(X1,1)/n1;s1=sum(s1.*s1,1)/(n1-1);
s2=X2-ones(n2,1)*sum(X2,1)/n2;s2=sum(s2.*s2,1)/(n2-1);
s=sqrt((n1-1)*s1+(n2-1)*s2);
t=(z1-z2).*(ones(1,k)*sqrt(n1*n2*(n-2)/n)./s);



