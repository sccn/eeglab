function [HDR] = ssave(FILENAME,DATA,TYPE,Fs,gdftyp)
% SSAVE saves signal data in various data formats
% 
% Currently are the following data formats supported: 
%    EDF, BDF, GDF, BKR, SND/AU, (WAV, AIF)
%    and WSCORE event file
%
% HDR = ssave(HDR,data);
% HDR = ssave(FILENAME,data,TYPE,Fs);
%
% FILENAME      name of file
% data  signal data, each column is a channel
% TYPE 	determines dataformat
% Fs	sampling rate	
%
% see also: SSAVE, SOPEN, SWRITE, SCLOSE, doc/README
%

% $Id: ssave.m,v 1.1 2009-07-07 02:23:48 arno Exp $
% Copyright (C) 2003,2004,2007 by Alois Schloegl <a.schloegl@ieee.org>	
% This file is part of the biosig project http://biosig.sf.net/

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



if isstruct(FILENAME),
        HDR = FILENAME;
        if isfield(HDR,'FileName'),
                FILENAME = HDR.FileName;
        else
                fprintf(2,'Error SSAVE: missing FileName.\n');	
                return; 
        end;
else
        HDR.FileName = FILENAME;
        HDR.SampleRate = Fs; 
        gdftyp = 16;
        HDR.TYPE = 'native';
end;

if (nargin > 1),
	[HDR.SPR, HDR.NS] = size(DATA); HDR.NRec = 1; 
%	HDR.AS = rmfield(HDR.AS,'SPR'); 
	if (strcmp(HDR.TYPE,'BDF') | strcmp(HDR.TYPE,'EDF') | strcmp(HDR.TYPE,'GDF')) & (~isfield(HDR,'DigMax') | ~isfield(HDR,'DigMin') |~isfield(HDR,'PhysMax') | ~isfield(HDR,'PhysMin'))
		HDR.PhysMax = max(DATA,[],1);
		HDR.PhysMin = min(DATA,[],1);
		ix = find(HDR.PhysMax == HDR.PhysMin);
		HDR.PhysMin(ix) = HDR.PhysMin(ix) - 1;
		if strcmp(HDR.TYPE,'BDF')
		   	[datatyp,HDR.THRESHOLD,datatypes,HDR.bits,HDR.GDFTYP] = gdfdatatype(511+24*ones(HDR.NS,1));
		   	HDR.DigMax = HDR.THRESHOLD(:,2)';
		   	HDR.DigMin = HDR.THRESHOLD(:,1)';
		elseif strcmp(HDR.TYPE,'EDF')
		   	[datatyp,HDR.THRESHOLD,datatypes,HDR.bits,HDR.GDFTYP] = gdfdatatype(3*ones(HDR.NS,1));
		   	HDR.DigMax = HDR.THRESHOLD(:,2)';
		   	HDR.DigMin = HDR.THRESHOLD(:,1)';
		elseif strcmp(HDR.TYPE,'GDF')
		   	[datatyp,HDR.THRESHOLD,datatypes,HDR.bits,HDR.GDFTYP] = gdfdatatype(16*ones(HDR.NS,1));
			HDR.DigMax = HDR.PhysMax;
			HDR.DigMin = HDR.PhysMin;
		else 	
			HDR.DigMax = HDR.PhysMax;
			HDR.DigMin = HDR.PhysMin;
		end; 
	end;    	

	if (nargin > 2),
        	if strcmp(TYPE,'GDF2'),
        		HDR.TYPE = 'GDF';
	        	HDR.VERSION = 2;
        	elseif strncmp(TYPE,'GDF',3),
        		HDR.TYPE = 'GDF';
      	 	 	HDR.VERSION = 1.25;
        	else	
	        	HDR.TYPE = TYPE; 	% type of data format
		end;        
	end;
	HDR = sopen(HDR,'w');
	HDR = swrite(HDR,DATA);
	HDR = sclose(HDR);
end;

% Convert EVENT into WSCORE event format
if all([length(HDR.EVENT.POS), length(HDR.EVENT.TYP)]),
	p = which('sopen'); [p,H,e] = fileparts(p);
	H = sload(fullfile(p,'../doc/eventcodes.txt'));

	HDR.EVENT.CodeDesc  = H.CodeDesc;
	HDR.EVENT.CodeIndex = H.CodeIndex;
	if isfield(HDR.EVENT,'DUR')
	        HDR.EVENT.POS = [HDR.EVENT.POS; HDR.EVENT.POS + HDR.EVENT.DUR];
	        HDR.EVENT.TYP = [HDR.EVENT.TYP; HDR.EVENT.TYP + hex2dec('8000')];
	end;
	OnOff = {'On','Off'};
	
	[HDR.EVENT.POS, ix] = sort(HDR.EVENT.POS);
	HDR.EVENT.TYP       = HDR.EVENT.TYP(ix);
	[TYP, IX, IY]       = unique(HDR.EVENT.TYP);

	% write "free form" scoring file for WSCORE
	fid   = fopen(fullfile(HDR.FILE.Path,[HDR.FILE.Name,'.C07']),'w');
	for k = 1:length(TYP), 
    		fprintf(fid,'%2i %s (%s)\r\n', k, HDR.EVENT.CodeDesc(mod(TYP(k),2^15)==HDR.EVENT.CodeIndex), OnOff{(TYP(k)>=2^15)+1});
	end;
	fclose(fid);

	% write "free form" scoring file for WSCORE
	fid = fopen(fullfile(HDR.FILE.Path,[HDR.FILE.Name,'.007']),'w');
	fprintf(fid,'%8i %i\r\n', [round(HDR.EVENT.POS(:)),IY(:)]');
	fclose(fid);
end;
