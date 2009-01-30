function [H2,HDR,s] = qrsdetect(fn,arg2,arg3)
% QRSDETECT - detection of QRS-complexes
%
%   HDR = qrsdetect(fn,chan,Mode)
%   HDR = qrsdetect(s,Fs,Mode)
%
% INPUT
%   fn	        filename
%   chan        channel number of ecg data
%   s           ecg signal data 
%   Fs          sample rate 
%   Mode        optional - default is 1
%               1: method [1] is used
%		2: method [2] is used	
% OUTPUT
%   HDR.EVENT  fiducial points of qrs complexes	
%
%
% see also: PROCESSING, EVENTCODES.TXT, SLOAD 
%
% Reference(s):
% [1] M.-E. Nygards, L. Sörnmo, Delineation of the QRS complex using the envelope of the e.c.g
%       Med. & Biol. Eng. & Comput., 1983, 21, 538-547.
% [2] V. Afonso, W. Tompkins, T. Nguyen, and S. Luo, "ECG beat detection using filter banks,"
% 	IEEE Trans. Biomed. Eng. 46(2):192-202, Feb. 1999.
%


%	$Id: qrsdetect.m,v 1.1 2009-01-30 06:04:44 arno Exp $
%	Copyright (C) 2000-2003,2006 by Alois Schloegl <a.schloegl@ieee.org>	
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


chan = 0; 
MODE = 1;       % default: 
if (nargin==2) 
	if isnumeric(arg2)
		chan = arg2;
	else
		MODE = arg3; 
	end;
elseif (nargin==3) 
	chan = arg2; 
	MODE = arg3;
end;

if isnumeric(fn),
	s = fn; 
	HDR.SampleRate = arg2;
	chan = 1:size(s,2); 
else
	HDR = sopen(fn,'r',chan);
	if chan==0;
		chan = strmatch('ecg',HDR.Label);
		if length(chan)==1,
			HDR = sclose(HDR);
			HDR = sopen(fn,'r',chan);
		else
			chan = 1:HDR.NS;	
		end;	
	end;	
	[s,HDR] = sread(HDR);
	HDR = sclose(HDR);
end;	

ET = [];
for k = 1:size(s,2),
	if 0,

        elseif MODE==2,     % QRS detection based on Afonso et al. 1999
                POS = nqrsdetect(s(:,k),HDR.SampleRate);

        elseif MODE==1,     % QRS detection based on Hilbert transformation. For details see [1]
                y   = processing({'ECG_envelope',HDR.SampleRate},s(:,k));
                TH  = quantile(y,.90);

                POS = gettrigger(y,TH);	% fiducial point

        else %if Mode== ???     % other QRS detection algorithms
		fprintf(2,'Error QRSDETECT: Mode=%i not supported',Mode);
		return; 
        end;

        % difference between R-peak and fiducial point
        [t,sz] = trigg(s(:,k),POS,floor(-HDR.SampleRate),ceil(HDR.SampleRate));
        [tmp,ix] = max(abs(mean(reshape(t,sz(2:3)),2)));
        delay = HDR.SampleRate + 1 - ix;

        ET = [ET; [POS-delay, repmat([hex2dec('0501'),chan(k),0], size(POS,1),1)]];
end;


[tmp,ix] = sort(ET(:,1));
if isfield(HDR,'T0')
	H2.T0 = HDR.T0; 
end;	
if isfield(HDR,'Patient')
	H2.Patient = HDR.Patient;
end;	
H2.EVENT.POS = ET(ix,1);
H2.EVENT.TYP = ET(ix,2);
H2.EVENT.CHN = ET(ix,3);
H2.EVENT.DUR = ET(ix,4);
H2.EVENT.SampleRate = HDR.SampleRate;
H2.TYPE = 'EVENT'; 
