function [AaH,pup,adpval,q]=exakteM_B(Input,n1,samp,gamma,B,test,tail,alpha,k);
%
%The function [AaH,pup,adpval,q] = exakteM_B(Input,n1,samp,gamma,B,test,tail,alpha,k)
%is a procedure used in the main-function fdp.m.
%The function exakteM_B is called 'Be' for the procedure B (conservative) of Korn et al. (2004) in the main-function fdp. 
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



format short;
%disp(sprintf('\nProzedur B - exakte Methode'));
AaH=0;
pup=0;

adpval=-ones(1,k);


M=B;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[P,M,q,Richtung]=pwerte(Input,n1,samp,M,test,tail,alpha);

schranke=floor(alpha*M)+1;

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if k>100

[Q,d]=sort(P(:,1));

Q=P;
    for j=1:k
        P(j,:)=Q(d(j),:);        % ordne die Original p-Werte der Größe nach
    end;%for


			% d merkt sich die Stelle der Original p-Werte

p=P(:,1);
mini=min(P);			% Zeilenvektor der spaltenweisen Minima der p-Werte
sortiert=sort(mini);		% sortiert diese Minima
y=sortiert(schranke);

	


else

[Q,d]=sort(P(:,1));
Q=P;
    for j=1:k
        P(j,:)=Q(d(j),:);        % ordne die Original p-Werte der Größe nach
    end;%for

p=P(:,1);
mini=min(P);			% Zeilenvektor der spaltenweisen Minima der p-Werte
sortiert=sort(mini);		% sortiert diese Minima
y=sortiert(schranke);
end;%if


if p(1)<y,
   %disp(sprintf('Lehnen H(1) ab (Parameter %g, p= %g)',d(1),q(d(1))));
   AaH=AaH+1;
   pup(AaH)=d(1);
   
    %Adj-p-Wert
    adpval(1)= sum(mini<=p(1)*ones(1,M))/M;
    %disp(sprintf('Komponente %g ist signifikant, p-Wert=%6.5f, adjustierter p-Wert=%6.5f',d(1),q(d(1)),adpval(1)));
    %pup(2,AaH)=q(d(1));
else %disp('Können keine Hypothese ablehnen');
    %size(P1),[anzahl1,AaH1,pup1,q1,adpval1]=exakteM(P1,M,n,k,pWert,fall,alpha,gamma);
return;
end;%if

for r=2:k,
   u=floor(r*gamma);
   if u>floor((r-1)*gamma),
        %disp(sprintf('Lehnen H(%g) automatisch ab (Parameter %g, p= %g)',r,d(r),q(d(r))));
        %disp(sprintf('Komponente %g ist signifikant(Automatische Ablehnung), p-Wert=%6.5f',d(r),q(d(r)) ));
      AaH=AaH+1;
   	pup(AaH)=d(r);
   	%pup(2,AaH)=q(d(r));
   else 
   
        if u==0,
            R=P(r:k,:);
     		mini=min(R,[],1);
            sortiert=sort(mini);			% sortiert diese Minima
		
		    %Adj-p-Wert
		    adpval(r)= sum(mini<=p(r)*ones(1,M))/M;
		    adpval(r)= max(adpval(1:r));
		
		    y=sortiert(schranke);
		    if p(r)<y,
                %disp(sprintf('Lehnen H(%g) ab (Parameter %g, p= %g)',r,d(r),q(d(r)))); 
   %             disp(sprintf('Komponente %g ist signifikant, p-Wert=%6.5f, adjustierter p-Wert=%6.5f',d(r),q(d(r)),adpval(r) ));
                AaH=AaH+1;
   	            pup(AaH)=d(r);
   	             %pup(2,AaH)=q(d(r));
  	        else 
	            %disp(sprintf('Können H(%g) nicht ablehnen (Parameter %g, p= %g)',r,d(r),q(d(r))));
   	       
                %size(P1),[anzahl1,AaH1,pup1,q1,adpval1]=exakteM(P1,M,n,k,pWert,fall,alpha,gamma);
            return;
         	end;%if
        else %if u~=0
	 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 j=nchoosek(r-1,u);%r-1,u
	
N=j*u;
	
	 grenze=1000;
	 
	     % wenn Anzahl Kombinationen zu groß, dann iterativ 
	     % Kombinationen bestimmen
	     
	        if j<=grenze
		        I=[nchoosek([1:r-1],u);ones(1,u)];
		    else
		        clear I; 
	        end;%if
	        
		%erste Kombination ist immer [1:u]
		x=[1:u];
		
		
		
		pval=0;y=1;
		
		    
	    
	            for i=1:j
		 
		            Q=zeros(k-r+1+u,M); 
		                for l=1:u
		                    Q(l,:)=P(x(l),:);
		                end;%for l=1:u
                  
		        %Nächste Kombination wird eingelesen
		            if j<grenze
		                x=I(i+1,:); 
		            else
		                [a,x]=nextcomb(x,r-1,u);
		            end;%if
		 
		        Q(u+1:k-r+1+u,:)=P(r:k,:);		  		 
		        Q=sort(Q);
                mini=Q(u+1,:);
		        sortiert=sort(mini); 
		        pval=max(sum(mini<=p(r)*ones(1,M))/M,pval);
		        y=min(sortiert(schranke),y); 
		        end;%for
	      
	        adpval(r)=max([adpval(1:r-1),pval]);
	      
	        
        
		%[u,y,p(r)]
				
            if p(r)<y
	        %disp(sprintf('Komponente %g ist signifikant, p-Wert=%6.5f, adjustierter p-Wert=%6.5f',d(r),q(d(r)),adpval(r)));
               
           AaH=AaH+1;
   	       pup(AaH)=d(r);
   		    %pup(2,AaH)=q(d(r));
				
	       else
	       %size(P1),[anzahl1,AaH1,pup1,q1,adpval1]=exakteM(P1,M,n,k,pWert,fall,alpha,gamma);
	       return;
	       end;%if p(r)<y			
        end;%if
    end;%if
    
    


end;%for
%size(P1),[anzahl1,AaH1,pup1,q1,adpval1]=exakteM(P1,M,n,k,pWert,fall,alpha,gamma);
