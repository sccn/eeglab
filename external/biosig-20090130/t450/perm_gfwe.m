function [p]=perm_gfwe(Input,n1,samp,B,tstat,tail)
%
%The function[P,M,q,Richtung] = pwerte(Input,n1,samp,M,test,tail,alpha) 
% is a subfunction in the main-function gFWE.m.
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         BEGINN Einstichprobenproblem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=B;
[n,k]=size(Input);
switch samp
    case 'single'

        X = Input;
        Y = zeros(n,k);
        D = X - Y;
        t=ttestC(D); 
	
	switch tstat
	   case 'tsum'
	       t0sum=t*ones(k,1);
	       switch tail
            	case '~='
            	 t0sum = abs(t0sum);	% t-Werte zweiseitig
            	case '>'
            	 t0sum = t0sum;		% t-Werte einseitig
            	case '<'
            	 t0sum = -t0sum;		% t-Werte einseitig
               end;%switch
               t0=t0sum;
	   case 'tsumabs'
		t0sumbetrag=abs(t)*ones(k,1);
		t0=t0sumbetrag;
	   case 'tmax'
	       switch tail
                case '~='
            	 t = abs(t);	% t-Werte zweiseitig
                case '>'
            	 t = t;		% t-Werte einseitig
                case '<'
            	 t = -t;		% t-Werte einseitig
               end;%switch
	       [y,lidx]=max(t);	% Maximum der Beträge der t-Werte
   	       t0max=t(lidx);	% max. t-Wert
	       t0=t0max;
    	end%switch tstat
    	M=min(M,2^n1);

    	%Auswahl der Zahlen 1:M zufaellig
    	g=zahlen(M-1,2^n1-1);

        for j=1:M-1,						% erzeugt M Permutationsmatrizen P
            h=dec2bin(g(j),n1);			% Umwandlung der Zahl g(k) in eine binaere
            Z=D;
                for i=1:n1,
                    if h(i)=='1',
                        Z(i,:)=-D(i,:);
                    end;%if
                end;%for
            t=ttestC(Z); 
	    switch tstat
	     case 'tsum'
	       t0sum=t*ones(k,1);
	       switch tail
            	case '~='
            	 t0sum = abs(t0sum);	% t-Werte zweiseitig
            	case '>'
            	 t0sum = t0sum;		% t-Werte einseitig
            	case '<'
            	 t0sum = -t0sum;		% t-Werte einseitig
               end;%switch
	       tp(j)=t0sum;
	     case 'tsumabs'
		t0sumbetrag=abs(t)*ones(k,1);
		tp(j)=t0sumbetrag;
	     case 'tmax'
	       switch tail
                case '~='
            	 t = abs(t);	% t-Werte zweiseitig
                case '>'
            	 t = t;		% t-Werte einseitig
                case '<'
            	 t = -t;		% t-Werte einseitig
               end;%switch
	       [y,lidx]=max(t);	% Maximum der Beträge der t-Werte
   	       t0max=t(lidx);	% max. t-Wert
	       tp(j)=t0max;
    	    end%switch tstat                  
        end%for j=1:M-1  
      p = (sum(tp>=t0)+1)/M;   
    
    case 'paired'
        
       X = Input(1:n1,1:k);
       Y = Input(n1+1:n,1:k);
       D = X - Y;
       t=ttestC(D); 
	
	switch tstat
	   case 'tsum'
	       t0sum=t*ones(k,1);
	       switch tail
            	case '~='
            	 t0sum = abs(t0sum);	% t-Werte zweiseitig
            	case '>'
            	 t0sum = t0sum;		% t-Werte einseitig
            	case '<'
            	 t0sum = -t0sum;		% t-Werte einseitig
               end;%switch
               t0=t0sum;
	   case 'tsumabs'
		t0sumbetrag=abs(t)*ones(k,1);
		t0=t0sumbetrag;
	   case 'tmax'
	       switch tail
                case '~='
            	 t = abs(t);	% t-Werte zweiseitig
                case '>'
            	 t = t;		% t-Werte einseitig
                case '<'
            	 t = -t;		% t-Werte einseitig
               end;%switch
	       %[tmax1,lidx]=max(t);	% Maximum der Beträge der t-Werte
   	       %t0max=tmax1(1);	% max. t-Wert
	       t0max=max(t);
	       t0=t0max;
    	end%switch tstat
    	M=min(M,2^n1);

    	%Auswahl der Zahlen 1:M zufaellig
    	g=zahlen(M-1,2^n1-1);

        for j=1:M-1,						% erzeugt M Permutationsmatrizen P
            h=dec2bin(g(j),n1);			% Umwandlung der Zahl g(k) in eine binaere
            Z=D;
                for i=1:n1,
                    if h(i)=='1',
                        Z(i,:)=-D(i,:);
                    end;%if
                end;%for
            t=ttestC(Z); 
	    switch tstat
	     case 'tsum'
	       t0sum=t*ones(k,1);
	       switch tail
            	case '~='
	       	 t0sum = abs(t0sum);	% t-Werte zweiseitig
            	case '>'
            	 t0sum = t0sum;		% t-Werte einseitig
            	case '<'
            	 t0sum = -t0sum;		% t-Werte einseitig
               end;%switch
	       tp(j)=t0sum;
	     case 'tsumabs'
		t0sumbetrag=abs(t)*ones(k,1);
		tp(j)=t0sumbetrag;
	     case 'tmax'
	       switch tail
                case '~='
            	 t = abs(t);	% t-Werte zweiseitig
                case '>'
            	 t = t;		% t-Werte einseitig
                case '<'
            	 t = -t;		% t-Werte einseitig
               end;%switch
	      % [tmax1,lidx]=max(t);	% Maximum der Beträge der t-Werte
   	      % t0max=tmax1(1);
	       t0max=max(t);
	       tp(j)=t0max;
    	    end%switch tstat                  
        end%for j=1:M-1 
      p = (sum(tp>=t0)+1)/M;    
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   ENDE Eistichprobenproblem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% beginn Zweistichprobenproblem indept
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'indept'
    
      [n,k]=size(Input);
       nk=nchoosek(n,n1);
       M=min(M,nk);
       X = Input(1:n1,1:k);
       Y = Input(n1+1:n,1:k);
       if  strcmp(tstat,'tsum') | strcmp(tstat,'tsumabs') | strcmp(tstat,'tmax')
          t=ttest3(Input,n1);
       end;
       switch tstat
	   case 'tsum'
	       t0sum=t*ones(k,1);
	       switch tail
            	case '~='
            	 t0sum = abs(t0sum);	% t-Werte zweiseitig
            	case '>'
            	 t0sum = t0sum;		% t-Werte einseitig
            	case '<'
            	 t0sum = -t0sum;		% t-Werte einseitig
               end;%switch
               t0=t0sum;
	   case 'tsumabs'
		t0sumbetrag=abs(t)*ones(k,1);
		t0=t0sumbetrag;
	   case 'tmax'
	       switch tail
                case '~='
            	 t = abs(t);	% t-Werte zweiseitig
                case '>'
            	 t = t;		% t-Werte einseitig
                case '<'
            	 t = -t;		% t-Werte einseitig
               end;%switch
	       [y,lidx]=max(t);	% Maximum der Beträge der t-Werte
   	       t0max=t(lidx);	% max. t-Wert
	       t0=t0max;
	   case 'ta'
	       n2=n-n1;
	       s1=X-ones(n1,1)*sum(X,1)/n1;s1=sum(s1.*s1,1)/(n1-1);
	       s2=Y-ones(n2,1)*sum(Y,1)/n2;s2=sum(s2.*s2,1)/(n2-1);
               %s=ones(n,1)*sqrt(s1.*s1/n1+s2.*s2/n2);

               sx=ones(n1,1)*sqrt(s1+s2/n2);
               sy=ones(n2,1)*sqrt(s1/n1+s2);
               mx=sum(X,1)/n1;my=sum(Y,1)/n2;
           
               %Normierung
               v=(Y-ones(n2,1)*mx)./sy;
               w=(X-ones(n1,1)*mx)./sx;

               z_v=sum(v,1);
               z_w=sum(w,1);

               %t_a
               t_a1=v-sum(z_v,2)*ones(n2,k)/(n2*k);
               t_a2=w-sum(z_w,2)*ones(n1,k)/(n1*k);
               ta=[t_a1/n2;t_a2/n1];
               ta=sum(ta.*ta,1);
	       t0=max(ta,[],2);
    	end%switch tstat         
 
    
%Nun werden M Permutationen erzeugt und die zugehï¿½rigen
%p-Werte in Matrix P gespeichert
    
        I=[1:n1];
        Liste=[0];
        if  M>(nk/2) & M<nk
            Z=1;
            g=setdiff([1:nk-1],zahlen(M-1,nk-1));
         end;%if
    
        for j=1:M-1
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%BEGINN     %Erzeugung von disjunkten Kombination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        
            if M<=(nk/2)
         
    %würfelt zufällig eine Kombination, testet ob Wiederholung 
         
                a=1;
                while a
                    [zahl,U]=umord(Input,n1);
          
                        if isempty(find(Liste==zahl))
                            Liste=[Liste;zahl];a=0;
                        end;%if
                end;%while
           
    %Wenn M==nk, dann werden alle Permutationen betrachtet
            elseif  M==nk
                [b,I]=nextcomb(I,n,n1);
                J=[I,setdiff([1:n],I)];
            %Umordnen
                for r=1:n
                    U(r,:)=Input(J(r),:); 
                end;%for
        
            else  %M<=(nk/2)
        
    %Ermittel der Kombinationen indem alle iterativ erzeugt werden
    %Dort wo E<=sprung, diese Kombinationen werden zufï¿½llig eleminiert
         
   
                a=1;
                while a
                    [b,I]=nextcomb(I,n,n1);
     
                        if isempty(find(g==Z))
                            a=0;
                        end;%if 
                    Z=Z+1;
                end;%while
   
                J=[I,setdiff([1:n],I)];
   
            %Umordnen
                for r=1:n
                    U(r,:)=Input(J(r),:); 
                end;%for
   
   
            end;%if M<=(nk/2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%BEGINN     %Erzeugung von disjunkten Kombination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

   %Berechnung der Teststatistik und p-Werte
   	      if  strcmp(tstat,'tsum') | strcmp(tstat,'tsumabs') | strcmp(tstat,'tmax')
                 t=ttest3(U,n1);
              end;
             switch tstat
	      case 'tsum'
	       t0sum=t*ones(k,1);
	       switch tail
            	case '~='
            	 t0sum = abs(t0sum);	% t-Werte zweiseitig
            	case '>'
            	 t0sum = t0sum;		% t-Werte einseitig
            	case '<'
            	 t0sum = -t0sum;		% t-Werte einseitig
               end;%switch
	       tp(j)=t0sum;
	     case 'tsumabs'
		t0sumbetrag=abs(t)*ones(k,1);
		tp(j)=t0sumbetrag;
	     case 'tmax'
	       switch tail
                case '~='
            	 t = abs(t);	% t-Werte zweiseitig
                case '>'
            	 t = t;		% t-Werte einseitig
                case '<'
            	 t = -t;		% t-Werte einseitig
               end;%switch
	       [y,lidx]=max(t);	% Maximum der Beträge der t-Werte
   	       t0max=t(lidx);	% max. t-Wert
	       tp(j)=t0max;
	     case 'ta'
	       X=U(1:n1,:);Y=U(n1+1:end,:); 
	       n2=n-n1;
	       s1=X-ones(n1,1)*sum(X,1)/n1;s1=sum(s1.*s1,1)/(n1-1);
	       s2=Y-ones(n2,1)*sum(Y,1)/n2;s2=sum(s2.*s2,1)/(n2-1);
               %s=ones(n,1)*sqrt(s1.*s1/n1+s2.*s2/n2);

               sx=ones(n1,1)*sqrt(s1+s2/n2);
               sy=ones(n2,1)*sqrt(s1/n1+s2);
               mx=sum(X,1)/n1;my=sum(Y,1)/n2;
           
               %Normierung
               v=(Y-ones(n2,1)*mx)./sy;
               w=(X-ones(n1,1)*mx)./sx;

               z_v=sum(v,1);
               z_w=sum(w,1);

               %t_a
               t_a1=v-sum(z_v,2)*ones(n2,k)/(n2*k);
               t_a2=w-sum(z_w,2)*ones(n1,k)/(n1*k);
               ta=[t_a1/n2;t_a2/n1];
              
               ta=sum(ta.*ta,1);
               
               
	       tp(j)=max(ta,[],2);
	       
    	    end%switch tstat   
        end;%for j =1:M-1
	p = (sum(tp>=t0)+1)/M;    
end;%switch samp 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ende Zweistichprobenproblem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
