function [AaH,pup,adpval] = bl01(psd,indexsd,k,alpha);
%
% The function [AaH,pup,adpval] = bl01(psd,indexsd,k,alpha) is a 
% subfunction in fdr.m.
% In fdr.m is bl01.m used as an procedure named 'B-L' for Benjamini and Liu
% (2001).
% For help use the main-function fdr.m.
%
% REFERENCES:
%
% [1] Hemmelmann C, Horn M, Suesse T, Vollandt R, Weiss S.
%	New concepts of multiple tests and their use for evaluating 
%	high-dimensional EEG data.
%	J Neurosci Methods. 2005 Mar 30;142(2):209-17.
%
%
% Copyright (C) 2006,2007 Claudia Hemmelmann <claudia.hemmelmann@mti.uni-jena.de>
% Adapted by A Schloegl <a.schloegl@ieee.org> 2006,2007 
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
%psd

for j=1:k
    c=min(1,alpha*k/(k-j+1)^2);	% krit. Wert für bl01
    if psd(j)<=c & psd(j) <= alpha,
       AaH=AaH+1;
       pup(AaH)=indexsd(j);
       
    else 
       return;
    end;
    adpval(j)=psd(j)/c*alpha;if j>1,adpval(j)=max(adpval(j),adpval(j-1));end;
end;

