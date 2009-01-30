function [AaH,pup,adpval] = bh95(psu,indexsu,k,alpha);
%
% The function [AaH,pup,adpval] = bh95(psd,indexsd,k,alpha) is a subfunction 
% in fdr.m.
% In fdr.m is bh95.m used as an procedure named 'B-H' for Benjamini and Hochberg.
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
%psu

for j=1:k,
    c=alpha*(k-j+1)/k;
    if psu(j)<=c,
       AaH=k-j+1;
       for i=1:k-j+1,
          pup(i)=indexsu(k-i+1);
       end;
       adpval=((k+1)*ones(1,k)-[1:k])/k;
       adpval=fliplr(psu./adpval);
       for i=k-j:-1:1
          adpval(i)=min(adpval(i),adpval(i+1));
       end;
       adpval=adpval(1:k-j+1);
       return;
    end;
end;
