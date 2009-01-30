function [AaH,pwert]=homhof(p,k,alpha,u);
%
%The function[AaH,pwert] = homhof(p,k,alpha,u)) is a subfunction in the main-function fdp.m.
%Use the main-function for help.
%
%Erweiterung von Holm für u=1,2,3... nach Hommel und Hoffmann (1987)
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



 AaH=0;
 pwert=zeros(k,1);
 
 % t-Test für Einstichprobenfall
 
 c=[ones(1,u+1)*(u+1)*alpha/k, (u+1)*alpha./[k-1:-1:u+1]];
 
 for j=1:k
    if p(j)>=c(j),
       AaH = j-1;
       return;
    end;
    if j==k,
       AaH = k;
       return;
    end;
 end; 
