function Q=arspectrum(s,arg2,arg3,impedance)
% Spectral analysis using autoregressive method. 
%
%  X = arspectrum(s,Fs,PhysDim)
%       s  	signal data
%       Fs      Sampling Rate
%       PhysDim   physical units of s
%
%  X = arspectrum(filename)
%  X = arspectrum(filename,CHAN)
%       filename must be a signal format known by BIOSIG
%       
%  Result can be displayed with PLOTA(Q)          
%
%  see also: SLOAD, PLOTA, TSA/LATTICE, TSA/DURLEV

%	$Id: arspectrum.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (C) 2005 by Alois Schloegl <a.schloegl@ieee.org>
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

                
Mode='log';

if ischar(s),
        if nargin==2,
                CHAN = arg2; 
                %[s,H]=sload(s,CHAN,'OVERFLOWDETECTION:OFF');
                [s,H]=sload(s,CHAN);
                Q = H; 
                Q.Label = H.Label(CHAN,:);
                Q.PhysDim = H.PhysDim(CHAN,:);
                Q.PhysDimCode = H.PhysDimCode(CHAN);
        else
                %[s,H]=sload(s,0,'OVERFLOWDETECTION:OFF');
                [s,H]=sload(s,0);
                Q = H; 
        end;
        Fs = H.SampleRate;
elseif isnumeric(s) 
        H.NS = size(s,2); 
        Q.SampleRate = arg2; 
        Q.PhysDim = arg3; 
        Q.PhysDimCode = physicalunits(Q.PhysDim); 
        Q.Label = [repmat('#',H.NS,1),int2str([1:H.NS]')];
end;

if isnumeric(s)
        h = histo3(s);
        %Q.HISTO = h;
        %Q.stats = hist2res(h);
        %Q = y2res(s);
        q = hist2res(h);
        Q.MEAN = q.MEAN;
        Q.RMS = q.RMS;
        Q.STD = q.STD;
        if isempty(q.QUANT)
                Q.QUANT = full(diag(H.Calib(H.InChanSelect+1,:)));
        else
                Q.QUANT = q.QUANT;
        end;
        Q.ENTROPY = q.ENTROPY;
	[Q.AR,Q.RC,Q.PE]  = lattice(center(s)',50);
	%[Q.AR,Q.RC,Q.PE] = lattice(center(s)',500);
	%[Q.AR,Q.RC,Q.PE] = durlev(acovf(center(s)',50));
        Q.datatype   = 'spectrum';
	%[Q.VAR,Q.VRC,Q.VPE] = mvar(center(s),50);
end;

