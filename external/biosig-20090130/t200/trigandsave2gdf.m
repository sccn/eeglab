function [signal,H] = trigandsave2gdf(FILENAME,TI1,CHAN,FN2)
% TrigAndSave2GDF loads and triggers signal data and saves it to a GDF-file   
%
% [signal,HDR] = trig&save2gdf(SourceFilename, TI, [CHAN,], TargetFilename)
%
% S = reshape(signal,HDR.size) returns the corresponding 3-dim Matrix 
%
%
% see also: SLOAD, TLOAD, SVIEW, SOPEN, 


%	$Id: trigandsave2gdf.m,v 1.1 2009-01-30 06:04:42 arno Exp $
%	Copyright (C) 2004-2005 by Alois Schloegl <a.schloegl@ieee.org>
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

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


if nargin<3; CHAN=0; end;
if isempty(CHAN), CHAN = 0; end;
if (CHAN<1) | ~isfinite(CHAN),
        CHAN=0;
end;

HDR = sopen(FILENAME,'r',CHAN);
if HDR.FILE.FID<0, 
        fprintf(2,'Error TRIG2GDF: file %s not found\n',FILENAME);
        return; 
end; 
if ~isfield(HDR,'TRIG') 
        fprintf(2,'Error TRIG2GDF: No Trigger information found in file %s\n',FILENAME);
        HDR = sclose(HDR);
        return; 
end; 
HDR.FLAG.UCAL = 1; 
HDR.FLAG.OVERFLOWDETECTION = 0; 
[s,HDR] = sread(HDR,inf);
HDR = sclose(HDR);


TI1 = TI1*HDR.SampleRate;
if length(TI1)<3, TI1(3)=0; end; 
if HDR.FLAG.TRIGGERED & (any(TI1<1) | any(TI1>HDR.SPR))
	fprintf(2,'Warning TLOAD: data is already triggered - invalid trigger interval\n');
	[signal,sz] = trigg(s,HDR.TRIG,1,HDR.SPR,TI1(2)-TI1(1)-HDR.SPR+TI1(3));
	signal = [signal(:,1+end+TI1(1):end),signal(:,1:end+TI1(1))];
else
        [signal,sz] = trigg(s,HDR.TRIG,TI1(1)+1,TI1(2),TI1(3));
end;		
signal = signal'; 

H = HDR; 
H.size = sz([2,3,1]);
H.FLAG.TRIGGERED = 1; 
H.NS   = sz(1);
H.SPR  = sz(2);
H.Dur  = sz(2)/H.SampleRate;
H.AS = []; %ones(1,HDR.NS)*sz(2); 
H.NRec = sz(3); 
H.TRIG = (0:H.NRec-1)'*H.SPR; 

H.EVENT = [];
H.EVENT.POS = H.TRIG;
H.EVENT.TYP = repmat(hex2dec('0300'),length(H.TRIG),1); 
if isfield(HDR,'Classlabel');
        cl = HDR.Classlabel; 
        cl(isnan(cl)) = 15;
        H.EVENT.POS = [H.EVENT.POS;H.EVENT.POS]; 
        H.EVENT.TYP = [H.EVENT.TYP;hex2dec('0300')+HDR.Classlabel(:)]; 
end;

H.FileName = FN2;
H.TYPE = 'GDF'; 

%save matlab H HDR signal
save2gdf(H,signal);
