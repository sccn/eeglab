function [AaH,pwert]=lehrom(p,k,alpha,gamma,vari);
%
%The function [AaH,pwert] = lehrom(p,k,alpha,gamma,vari)
%is a procedure used in the main-function fdp.m.
%There are two variants of this procedure: LR1 and LR2.
%Use the main-function for help.
%
%
%
%Copyright (C) 2006 by Claudia Hemmelmann <claudia.hemmelmann@mti.uni-jena.de>
%Institute of Medical Statistics, Computer Sciences and Documantation
%University of Jena
%This work was supported by DFG Project VO 683/2-1
%This is part of the BIOSIG-toolbox http://biosig.sf.net/
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


% Erweiterung von Holm für gamma in (0, 1) nach Lehmann und Romano (2004)


 AaH=0;
 pwert=zeros(k,1);
 
 % t-Test für Einstichprobenfall
 
 if vari==1,
   c=(floor(gamma.*[1:k])+1)*alpha./(k+floor(gamma.*[1:k])+1-[1:k]);
 else
   sum_bound = floor(k*gamma)+1;
   summe = sum(1./[1:1:sum_bound]);
   c=((floor(gamma.*[1:k])+1)*alpha./(k+floor(gamma.*[1:k])+1-[1:k]))./summe;
 end;
 
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
