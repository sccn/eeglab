function [I,ERR,SNR]=mutinfo(d,c);
% MUTINFO mutual information 
% [I,ERR,SNR]=mutinfo(d,c);
%   d   output DATA (each column is one trial) 
%   c   CLASS, vector with 0 and 1 
% 
% function [I,ERR,SNR]=mutinfo(d1,d0);
%   d1  DATA of class 1 
%   d2  DATA of class 0
% 
% OUTPUT:
%   I 	mutual information
%   ERR	Error rate;
%   SNR	signal to noise ratio  

% References:
% [1] Schloegl A., Neuper C. Pfurtscheller G.
% Estimating the mutual information of an EEG-based Brain-Computer-Interface
% Biomedizinische Technik 47(1-2): 3-8, 2002

%	$Id: mutinfo.m,v 1.1 2009-01-30 06:04:50 arno Exp $
%	Copyright (c) 1997-2003,2007 by  Alois Schloegl
%	a.schloegl@ieee.org	

%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% version 2 of the License, or (at your option) any later version.
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


c=c(:).';
ntr=length(c);
[nr,nc]=size(d);

CL=unique(c);
if length(CL)~=2, 
        error('invalid class label');
end;

% time course of the SNR+1, Mutual Information, and the Error rate [1]
SNRp1 	= 2*var(d,[],2)./(var(d(:,c==CL(1)),[],2)+var(d(:,c==CL(2)),[],2));
I 	= log2(SNRp1)/2;
ERR 	= 1/2 - mean(sign(d).*sign(c(ones(nr,1),:)-mean(CL)),2)/2;

SNR     = SNRp1-1;
