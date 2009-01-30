function [B,A]=lumped(K,Fs)
% Transfer functions of the lumped alpha model
% [B,A]=lumped(K [,Fs])
%
%  [B,A] are the denominator and nominator, resp. , for a 
%       the coupling parameter K and the sampling frequency Fs [default 128Hz].
%
% Reference)(s)
% [1]	Lopes da Silva FH, Hoeks A, Smits H, Zetterberg LH.
%	Model of brain rhythmic activity. The alpha-rhythm of the thalamus.
%       Kybernetik. 1974 May 31;15(1):27-37.
% [2]   P. Suffcynski, Thesis, 1999.
% [3]   Alois Schlögl (2000)
%       The electroencephalogram and the adaptive autoregressive model: theory and applications
%       Shaker Verlag, Aachen, Germany,(ISBN3-8265-7640-3). 


%       $Revision: 1.1 $
%	$Id: lumped.m,v 1.1 2009-01-30 06:04:44 arno Exp $
%	Copyright (C) 1999-2004 by Alois Schloegl <a.schloegl@ieee.org>
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

c1=6;
c2=10;
c3=15;
c4=10;

A=1.6;   % mV
B=3.2;   % mV

a1=55;   % 1/s
a2=605;  % 1/s
b1=27.5; % 1/s
b2=55;   % 1/s

lambda_g0=25;   % 1/s

q =1.5;  % 1/mV
Vd =7;   % mV

P1t_off=312; % pps
P2t_off=312; % pps

P1t_var=169; % pps²
P2t_var=169; % pps²

Mt = 8; 


if nargin<2
        Fs=128; % Hz
end;

if nargin<1
        K=c1*c2*q*q*(a2-a1)*(b2-b1)*A*B;
end;



x1=conv([1 b1],[1 b2]);
x2=conv([1 a1],[1 a2]);

Z = A*(a2-a1)*x1;
P = conv(x2,x1);
P(5)=P(5)+K;

[B,A] = bilinear(Z,P,Fs) ;


