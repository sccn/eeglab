function [AaH,pup,adpval,q]=vereinM_B(Input,n1,samp,gamma,B,test,tail,alpha,k);
%
%The function [AaH,pup,adpval,q] = vereinM_B(Input,n1,samp,gamma,B,test,tail,alpha,k)
%is a procedure used in the main-function fdp.m.
%The function vereinM_B is called 'Bv' for the procedure B (conservative) of Korn et al. (2004) in the main-function fdp. 
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
%disp(sprintf('\nProzedur B - vereinfachte Methode'));
AaH=0;
pup=0;

M=B;

%adjustierter p-Wer==-1 ==> Automatische Ablehnung
adpval=-ones(1,k);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
[P,M,q,Richtung]=pwerte(Input,n1,samp,M,test,tail,alpha);

%size(P),P
schranke=floor(alpha*M)+1;

%schranke,M

d=[1:k]; 


[Q,d]=sort(P);	% ordne die Original p-Werte der Größe nach
d=d(:,1);%d ist Indexvektor der original p-Werte


p=Q(:,1);
mini=Q(1,:);			% Zeilenvektor der spaltenweisen Minima der p-Werte
sortiert=sort(mini);		% sortiert diese Minima
y=sortiert(schranke);

if p(1)<y,
   %disp(sprintf('Lehnen H(1) ab (Parameter %g, p= %g)',d(1),q(d(1))));
   AaH=AaH+1;
   pup(AaH)=d(1);
   %pup(2,AaH)=q(d(1));
   
   %Adj-p-Wert
   adpval(1)=sum(mini(1:M)<=p(1)*ones(1,M))/M;
   
   %disp(sprintf('Komponente %g ist signifikant, p-Wert=%6.5f, adjustierter p-Wert=%6.5f',d(1),q(d(1)),adpval(1) ));
else %disp('Können keine Hypothese ablehnen');
     %[AaH1,pup1] = vereinMalt1(P,M,n,k,pWert,alpha,gamma);
return;
end;%if

v=0;
for r=2:k,
    u=floor(r*gamma);
    if u>floor((r-1)*gamma),
      %disp(sprintf('Lehnen H(%g) automatisch ab (Parameter %g, p= %g)',r,d(r),q(d(r))));
      %disp(sprintf('Komponente %g ist signifikant(Automatische Ablehnung), p-Wert=%6.5f',d(r),q(d(r)) ));
        AaH=AaH+1;
   	    pup(AaH)=d(r);
   	  %pup(2,AaH)=q(d(r));
    else   
   
N=5*M*k/4;
         
    if u~=v,
        
	       
	    mini=Q(u+1,:);
	    sortiert=sort(mini);   % sortiert diese Minima
	    y=sortiert(schranke);%[y,u,p(r)]
	    v=u;
	    
    end;%if u~=v
	    
	    %[y,u,p(r)]
	     
	    %Adj-p-Wert
            adpval(r)= sum(mini(1:M)<=p(r)*ones(1,M))/M;
            adpval(r)= max(adpval(1:r));
                     
    if p(r)<=y,
	    %disp(sprintf('Komponente %g ist signifikant, p-Wert=%6.5f, adjustierter p-Wert=%6.5f',d(r),q(d(r)),adpval(r) ));
        %disp(sprintf('Lehnen H(%g) ab (Parameter %g, p= %g)',r,d(r),q(d(r))));
          AaH=AaH+1;
          pup(AaH)=d(r);
          %pup(2,AaH)=q(d(r));
	  
    else %disp(sprintf('Können H(%g) nicht ablehnen (Parameter %g, p= %g)',r,d(r),q(d(r))))
	  %[AaH1,pup1] = vereinMalt1(P,M,n,k,pWert,alpha,gamma);
	  return;
    end;%if
    end;
end;%for

%[AaH1,pup1] = vereinMalt1(P,M,n,k,pWert,alpha,gamma);
			
                    
                      

