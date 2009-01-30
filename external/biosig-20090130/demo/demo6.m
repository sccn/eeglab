%Demo6 - demonstrates the transfer function of the 
%        lumped circuit model for various feedback gains
%
%       see also LUMPED
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
%	$Id: demo6.m,v 1.1 2009-01-30 06:04:39 arno Exp $
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

LEG = [];
k = 0;
for K = 10.^[0:.25:9],
        k = k +1;
        LEG{k}=sprintf('%5f',K);
        [B,A]=lumped(K,128);
        [H,F]=freqz(B,A,0:.1:64,128); 
        semilogy(F,abs(H)); 
        hold on; 
end;
axis([0,64,1e-4,1])
%legend(LEG)
