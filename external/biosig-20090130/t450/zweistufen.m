function  [AaH,pup,adpval]= zweistufen(psu,indexsu,k,alpha)
%
% The function [AaH,pup,adpval] = zweistufen(psd,indexsd,k,alpha) is a 
% subfunction in fdr.m.
% In fdr.m is zweistufen.m used as an procedure named 'BKY' for 
% Benjamini. Krieger and Yekutieli (2001).
% For help use the main-function fdr.m.
%
% REFERENCES:
%
% [1] Hemmelmann C, Horn M, Suesse T, Vollandt R, Weiss S.
%	New concepts of multiple tests and their use for evaluating 
%	high-dimensional EEG data.
%	J Neurosci Methods. 2005 Mar 30;142(2):209-17.
%
% [2]
%
% Copyright (C) 2006 Claudia Hemmelmann <claudia.hemmelmann@mti.uni-jena.de>
% Adapted by A Schloegl <a.schloegl@ieee.org> Dec 2006 
%
%***
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
pup=0;
adpval=0;
q=alpha/(alpha+1);
for j=1:k,
    c1=q*(k-j+1)/k;			% krit. Wert im 1.Schritt
    if psu(j) <= c1,
    
       if j == 1,
          AaH=k;
	  for i=1:k,
             pup(i)=indexsu(k-i+1);
	  end;
	  adpval=((k+1)*ones(1,k)-[1:k])/k;
          adpval=fliplr(psu./adpval);
          for i=k-1:-1:1
            adpval(i)=min(adpval(i),adpval(i+1));
	  end;
	  
	  return;
       end;
       
       %Berechnung der adjustierten p-Werte der bereits 
       %abgelehnten Hypothesen
       
       adpval=((k+1)*ones(1,k)-[1:k])/k;
       adpval=fliplr(psu./adpval);
       for i=k-j:-1:1
         adpval(i)=min(adpval(i),adpval(i+1));
       end;
       adpval=adpval(1:k-j+1);
       
       
       
       %wenn psu(j)<=c1 und j ~=1
       for i=1:k,
	   c2=q*(k-i+1)/(j-1);		% krit. Wert im 2.Schritt
	
	   if psu(i) <= c2 & psu(i) <= alpha,
	      AaH=k-i+1;
             for o=1:k-i+1,
                 pup(o)=indexsu(k-o+1);
	     end;
	      
	     if i<j
	      
	      pval=[k-i:-1:k-j+1]/(j-1);%j,i
	      pval=psu(j-1:-1:i)./pval;
	     
              adpval=[adpval,ones(1,j-i)*alpha];
	      
	     end;
	     return;
	   	  
	   
	   end;
       
       end;
    end;
end;
