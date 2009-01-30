function [AaH,pup,adpval,q]=vereinM_A(Input,n1,samp,u,B,test,tail,alpha,k);
%
%The function [AaH,pup,adpval,q] = vereinM_A(Input,n1,samp,u,B,test,tail,alpha,k)
%is a procedure used in the main-function fdp.m.
%The function vereinM_A is called 'Av' for the procedure A (conservative) of Korn et al. (2004) in the main-function fdp. 
%Use the main-function for help.
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




format short;
AaH=0;
pup=0;



M=B;

%adjustierter p-Wer==-1 ==> Automatische Ablehnung
adpval=-ones(1,k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[P,M,q,Richtung]=pwerte(Input,n1,samp,M,test,tail,alpha);

schranke=floor(alpha*M)+1;

%schranke,M

d=[1:k]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Q,d]=sort(P);	% ordne die Original p-Werte der Größe nach
d=d(:,1);%d ist Indexvektor der original p-Werte
p=Q(:,1);

if u<=k

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=5*M*k/4;
    	       
mini=Q(u+1,:);
sortiert=sort(mini);   % sortiert diese Minima
y=sortiert(schranke);%[y,u,p(r)]

pup(1:u)=d(1:u);
AaH=u;	    
	    
    for r=u+1:k,	    
	    %Adj-p-Wert
            adpval(r)= sum(mini(1:M)<=p(r)*ones(1,M))/M;
            adpval(r)= max(adpval(1:r));
            
        if p(r)<=y,
       
          %disp(sprintf('Lehnen H(%g) ab (Parameter %g, p= %g)',r,d(r),q(d(r))));
          AaH=AaH+1;
          pup(AaH)=d(r);
          %pup(2,AaH)=q(d(r));
        else %disp(sprintf('Können H(%g) nicht ablehnen (Parameter %g, p= %g)',r,d(r),q(d(r))));
        return;
        end;%if	 
    end;%for

else
    
AaH=k;pup=d;
end;%if u<k

                      

