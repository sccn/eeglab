function HDR = qc_histo(fn,arg2)
% QC_HISTO performs quality control using histogram and entropy analysis
%
%  R = qc_histo(filename [,CHAN]); 
%
%  plota(R); displays the result as in [1]
%
% References: 
% [1] A. Schlögl, B. Kemp, T. Penzel, D. Kunz, S.-L. Himanen, A. Värri, G. Dorffner, G. Pfurtscheller.
%       Quality Control of polysomnographic Sleep Data by Histogram and Entropy Analysis.
%       Clin. Neurophysiol. 1999, Dec; 110(12): 2165 - 2170.
%       http://dx.doi.org/10.1016/S1388-2457(99)00172-8


%  $Id: qc_histo.m,v 1.1 2009-01-30 06:04:43 arno Exp $ 
%  Copyright (C) 2006 by Alois Schloegl <a.schloegl@ieee.org>
%  This is part of the BIOSIG-toolbox http://biosig.sf.net/
%

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


MODE=[];
if nargin<2,
        CHAN=0; 
elseif ~isnumeric(arg2)
        CHAN=0;
        MODE=arg2;
else 
        CHAN = arg2; 
end;

HDR = sopen(fn,'r',CHAN,MODE); 
[s,HDR] = sread(HDR); 
sclose(HDR); 

HDR.HIS = histo3(s);
HDR.RES = hist2res(HDR.HIS); 
%HDR.RES = y2res(s); 
%HDR.Threshold = repmat([-1,1]/1e2,HDR.NS,1);
HDR.datatype = 'qc:histo'; 

[HDR.AR,HDR.RC,HDR.PE]  = lattice(center(s)',50);


