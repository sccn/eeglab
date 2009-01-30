function [AaH,pup,adpval,q]=exakteM_A(Input,n1,samp,u,B,test,tail,alpha,k);
%
%The function [AaH,pup,adpval,q] = exakteM_A(Input,n1,samp,u,B,test,tail,alpha,k)
%is a procedure used in the main-function gFWE.m.
%The function exakteM_A is called 'Ae' for the procedure A-exact in the main-function gFWE. 
%Use the main-function for help.
% Korn et al., exakte Methode
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
%disp(sprintf('\nProzedur A - exakte Methode'));
AaH=0;
pup=0;

adpval=-ones(1,k);


M=B; %number of permutations


 [P,M,q,Richtung]=pwerte(Input,n1,samp,M,test,tail,alpha);
schranke=floor(alpha*M)+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P1=P;%size(P1)

[Q,d]=sort(P(:,1));
Q=P;
for j=1:k
P(j,:)=Q(d(j),:);        % ordne die Original p-Werte der Größe nach
end;%for

p=P(:,1);

if u<=k
    if u==0
        for r=1:k
        Q=P(r:k,:);
        mini=min(Q,[],1);%size(mini),M
        sortiert=sort(mini); 
        adpval(r)=sum(mini<=p(r)*ones(1,M))/M;
        y=sortiert(schranke); 
            if r>1
            adpval(r)=max(adpval(r),adpval(r-1));
            end;%if
            if p(r)<y
            %disp(sprintf('Komponente %g ist signifikant, p-Wert=%6.5f, adjustierter p-Wert=%6.5f',d(r),q(d(r)),adpval(r)));
            AaH=AaH+1;
            pup(AaH)=d(r);
            else
            return;
            end;%if
        end;%for
    else
        %for i=1:u
        %disp(sprintf('Komponente %g ist signifikant(Automatische Ablehnung), p-Wert=%6.5f',d(i),q(d(i)) ));
        %end;%for
    pup(1:u)=d(1:u);
    AaH=u;
        for r=u+1:k,
	 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        j=nchoosek(r-1,u);%r-1,u
        N=j*u;
        groesser=(N*M>=10000);
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
		
            if groesser
                for i=1:j
		 
                Q=zeros(k-r+1+u,M); 
                    for l=1:u
                    Q(l,:)=P(x(l),:);
                    end;%for
		
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
            else %if groesser
	      
                for i=1:j
		 
                Q=zeros(k-r+1+u,M); 
                    for l=1:u
                    Q(l,:)=P(x(l),:);
                    end;%for
		 
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
	      
            end;%if groesser
            
            if p(r)<y
            %disp(sprintf('Komponente %g ist signifikant, p-Wert=%6.5f, adjustierter p-Wert=%6.5f',d(r),q(d(r)),adpval(r))); 
            AaH=AaH+1;
            pup(AaH)=d(r);	
            else
            %size(P1),[anzahl1,AaH1,pup1,q1,adpval1]=exakteM(P1,M,n,k,pWert,fall,alpha,gamma);
	        return;
            end;%if			     
        end;%for r
    end;%if u==0    
else   

    %for i=1:k
    %disp(sprintf('Komponente %g ist signifikant(Automatische Ablehnung), p-Wert=%6.5f',d(i),q(d(i)) ));
    %end;%for
AaH=k;pup=d;
end;%if u<k





%size(P1),[anzahl1,AaH1,pup1,q1,adpval1]=exakteM(P1,M,n,k,pWert,fall,alpha,gamma);
