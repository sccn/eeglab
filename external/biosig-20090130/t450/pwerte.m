function [P,M,q,Richtung]=pwerte(Input,n1,samp,B,test,tail,alpha)
%
%The function[P,M,q,Richtung] = pwerte(Input,n1,samp,M,test,tail,alpha) is a subfunction in the main-function fdp.m.
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
%    BioSig is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BioSig is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BioSig.  If not, see <http://www.gnu.org/licenses/>.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         BEGINN Einstichprobenproblem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = B;
[n,k]=size(Input);
switch samp
    case 'single'

        X = Input;
        Y = zeros(n,k);
        D = X - Y;
        switch test           %chose test   
            case 'ttest'        %t-Test
                t=ttestC(D); 
                FG=n1-1;
                q=tcdf(t,n1-1);
            case 'wilcox'     %Wilcoxen signed rank test
                q=zeros(1,k);
                    for j=1:k,
                        [q(j),h]= wilcoxon_test(X(:,j),Y(:,j),alpha,tail); % Wilcoxentest fuer gepaarte Stichproben
                    end;%for
            case 'sign'       %sign test
                q=zeros(1,k);
                    for j=1:k,
                        [q(j),h]= signtest(X(:,j),Y(:,j),alpha,tail); %Vorzeichentest
                    end;%for
        end %switch test sing

    Richtung=[q<=0.5];
    
    if strcmp(test,'ttest')
        switch tail
            case '~='
            q=2*min(q,1-q);	% p-Werte zweiseitig
            case '>'
            q=1-q;	% p-Werte einseitig
            case '<'
            q=q;	% p-Werte einseitig
        end;%switch
    end%if
    P(:,1)=q';
    M=min(M,2^n1);

    %Auswahl der Zahlen 1:M zufaellig
    g=zahlen(M-1,2^n1-1);


        for j=1:M-1,						% erzeugt M Permutationsmatrizen P
            h=dec2bin(g(j),n1);			% Umwandlung der Zahl g(k) in eine binaere
            Z=Input; 
            X1=X;
            Y1=Y;
                for i=1:n1,
                    if h(i)=='1',
                        Z(i,:)=-Input(i,:);
                        H(i,:)=X1(i,:);
                        X1(i,:)=Y1(i,:);
                        Y1(i,:)=H(i,:);    
                    end;%if
                end;%for
        
                switch test
                    case 'ttest'
                    t=ttestC(Z);
                    p=tcdf(t,n1-1);
      
                    case 'wilcox'     %Wilcoxen signed rank test
                    p=zeros(1,k);
                        for e=1:k,
                        [p(e),h]= wilcoxon_test(X1(:,e),Y1(:,e),alpha,tail); % Wilcoxentest fuer gepaarte Stichproben
                        end;%for
                    case 'sign'       %sign test
                    p=zeros(1,k);
                        for e=1:k,
                        [p(e),h]= signtest(X1(:,e),Y1(:,e),alpha,tail); %Vorzeichentest
                        end;%for
                end %switch test sing     
         
                if strcmp(test,'ttest')
                switch tail
                    case '~='
                    p=2*min(p,1-p);	% p-Werte zweiseitig
                    case '>'
                    p=1-p;	% p-Werte einseitig
                    case '<'
                    p=p;	% p-Werte einseitig
                end;%switch
            end%if
            
            P(:,j+1)=p';
   
        end%for j=1:M-1  
   
    
    case 'paired'
        
    X = Input(1:n1,1:k);
    Y = Input(n1+1:n,1:k);
    D = X - Y;
        
        switch test          %chose test  
            case 'ttest'    %t-Test
                t=ttestC(D); 
                FG=n1-1;
                %einseitiger p-Wert des t-Tests  
                q=tcdf(t,FG);
            case 'wilcox'   %Wilcoxen signed rank test
                q=zeros(1,k);
                    for j=1:k,
                        [q(j),h]= wilcoxon_test(X(:,j),Y(:,j),alpha,tail); % Wilcoxentest fuer gepaarte Stichproben
                    end;%for
            case 'sign'     %sign test
                q=zeros(1,k);
                    for j=1:k,
                        [q(j),h]= signtest(X(:,j),Y(:,j),alpha,tail); %Vorzeichentest
                    end;%for
        end %switch test pair
                  
    Richtung=[q<=0.5];
    
        if strcmp(test,'ttest')
        switch tail
            case '~='
            q=2*min(q,1-q);	% p-Werte zweiseitig
            case '>'
            q=1-q;	% p-Werte einseitig
            case '<'
            q=q;	% p-Werte einseitig
        end;%switch
    end%if
    
    P(:,1)=q';



    M=min(M,2^n1);

    %Auswahl der Zahlen 1:M zufaellig
    g=zahlen(M-1,2^n1-1);


        for j=1:M-1,						% erzeugt M Permutationsmatrizen P
            h=dec2bin(g(j),n1);			% Umwandlung der Zahl g(k) in eine binaere
            Z=D;
            X1=X;
            Y1=Y;
                for i=1:n1,
                    if h(i)=='1',
                        Z(i,:)=-D(i,:);
                        H(i,:)=X1(i,:);
                        X1(i,:)=Y1(i,:);
                        Y1(i,:)=H(i,:);                
                    end;%if
                end;%for
                
                switch test
                    case 'ttest'
                    t=ttestC(Z);
                    p=tcdf(t,n1-1);
                    case 'wilcox'   %Wilcoxen signed rank test
                    p=zeros(1,k);
                        for e=1:k,
                        [p(e),h]= wilcoxon_test(X1(:,e),Y1(:,e),alpha,tail); % Wilcoxentest fuer gepaarte Stichproben
                        end;%for
                    case 'sign'     %sign test
                    p=zeros(1,k);
                        for e=1:k,
                        [p(e),h]= signtest(X1(:,e),Y1(:,e),alpha,tail); %Vorzeichentest
                        end;%for
                end %switch test pair
                      
                if strcmp(test,'ttest')
                switch tail
                    case '~='
                    p=2*min(p,1-p);	% p-Werte zweiseitig
                    case '>'
                    p=1-p;	% p-Werte einseitig
                    case '<'
                    p=p;	% p-Werte einseitig
                end;%switch
            end%if
   
            P(:,j+1)=p';
   
        end%for j=1:M-1  
            
        
        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   ENDE Eistichprobenproblem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% beginn Zweistichprobenproblem indept
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'indept'
    
        
    
    [n,k]=size(Input);
    n2=n-n1;
    FG=n-2;
    nk=nchoosek(n,n1);
    M=min(M,nk);
    %n,n1
        X = Input(1:n1,1:k);
        Y = Input(n1+1:n,1:k);
        
        switch test %chose test
            case 'ttest'        %t-Test 
	    
                t=ttest3(Input,n1);
                FG=n-2;
                %einseitiger p-Wert des t-Tests 
                q=tcdf(t,FG);
            case 'wilcox'       %Wilcoxen rank sum test
                q=zeros(1,k);	
                    for j=1:k,
                        [q(j),h] = u_test(X(:,j),Y(:,j),tail);
                    end;
        end %switch test indept
    
    Richtung=[q<=0.5];
    
    if strcmp(test,'ttest')
        switch tail
            case '~='
            q=2*min(q,1-q);	% p-Werte zweiseitig
            case '>'
            q=1-q;	% p-Werte einseitig
            case '<'
            q=q;	% p-Werte einseitig
        end;%switch
    end%if
        
    P(:,1)=q';
    
%Nun werden M Permutationen erzeugt und die zugehoerigen
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
         
    %wuerfelt zufaellig eine Kombination, testet ob Wiederholung 
         
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
    %Dort wo E<=sprung, diese Kombinationen werden zufaellig eleminiert
         
   
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
   
            switch test
                case 'ttest'
   
                t=ttest3(U,n1)';
                p=tcdf(abs(t),FG);
        
                case 'wilcox'       %Wilcoxen rank sum test
                X = U(1:n1,:);
                Y = U(n1+1:end,:);
                    
                    p=zeros(1,k);	
                    for e=1:k,
                        [p(e),h] = u_test(X(:,e),Y(:,e),tail);
                    end;
            end %switch test indept
  
        if strcmp(test,'ttest')
            switch tail
                case '~='
                p=2*min(p,1-p);	% p-Werte zweiseitig
                case '>'
                p=1-p;	% p-Werte einseitig
                case '<'
                p=p;	% p-Werte einseitig
            end;%switch
        end%if
   %%%%%%%%
  
        P(:,j+1)=p';
      
        end;%for j =1:M-1

end;%switch samp 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ende Zweistichprobenproblem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
