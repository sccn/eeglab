function [argout,H1,h2] = hdr2ascii(source,dest)
% HDR2ASCII converts the header information into ASCII text. 
%
%   HDR2ASCII(HDR [, ...]);
%	converts file header HDR 
%   HDR2ASCII(file [, ...]);
%	converts header of file 
%   HDR2ASCII(arg,dest_file);
%	converts file header HDR and writes it into dest_file
%   HDR=HDR2ASCII(...);
%	returns header HDR 
%  
% see also: SLOAD, SOPEN

%	$Id: hdr2ascii.m,v 1.1 2009-07-07 02:23:46 arno Exp $
%	Copyright (C) 2007,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

if nargin<2,
	if isstruct(source) 
		HDR = source; 
	elseif ischar(source)
		HDR = sopen(source);
		HDR = sclose(HDR); 
	else
		'not implemented yet',
	end;	
	dest = [tempname,'.dlm']; 
elseif isstruct(source) && ischar(dest); 
	HDR = source; 
elseif ischar(source) && ischar(dest); 
	HDR = sopen(source); 
	HDR = sclose(HDR); 
else
	'not implemented yet',
end; 
if isnan(HDR.NS) && exist('mexSLOAD','file'),
	[s,HDR]=mexSLOAD(HDR.FileName); 
end

if nargin>1,
	fid = fopen(dest,'wt'); 
	if fid<0,
		fprintf(2,'ERROR HDR2ASCII: could not open file %s\n',dest);
		return; 
	end; 
else
	fid = 1; 
end; 


%%%%%%%%% FIXED HEADER %%%%%%%%%%%%%%		
fprintf(fid,'[BioSig Header]\n\n'); 
fprintf(fid,'Version=0.10\n'); 
fprintf(fid,'generated=%04i-%02i-%02i %02i:%02i:%04.1f\n',datevec(now)); 
if fid>2;
	fprintf(fid,'\n;This is a TAB-delimiter file. When you edit this file, make sure not to corrupt the TABs ASCII(9)!\n\n'); 
	fprintf(fid,'ThisFile=%s\n',dest); 
end; 


fprintf(fid,'\n[Fixed Header]\n'); 
fprintf(fid,'Filename\t= %s\n',HDR.FileName); 
fprintf(fid,'Format  \t= %s\n',HDR.TYPE); 
if isfield(HDR.FILE,'size'), fprintf(fid,'SizeOfFile\t= %i\n',HDR.FILE.size); end;
if ~isfield(HDR,'NS')
	HDR.NS = 0; 	
end;	 
fprintf(fid,'NumberOfChannels\t= %i\n',HDR.NS);
if isfield(HDR,'SampleRate')
	fprintf(fid,'SamplingRate    \t= %i\n',HDR.SampleRate); 
end;	 
if isfield(HDR,'NRec') && isfield(HDR,'SPR')
	fprintf(fid,'Number_of_Samples\t= %i\n',HDR.NRec*HDR.SPR); 
end;	 
T0 = zeros(1,6);
if isfield(HDR,'T0')
	switch length(HDR.T0)
	case 1, T0 = datevec(HDR.T0);
	case 6, T0 = HDR.T0;
	end; 
	fprintf(fid,'RecordingDateTime\t= %04i-%02i-%02i %02i:%02i:%06.3f\n',T0);
end; 	 
if isfield(HDR,'Patient')
	fprintf(fid,'Patient.\n'); 
	if isfield(HDR.Patient,'Name')
		fprintf(fid,'\tName      \t= %s\n',HDR.Patient.Name); 
	end;
	if isfield(HDR.Patient,'Id')
		fprintf(fid,'\tId\t\t= %s\n',HDR.Patient.Id); 
	end;
	if isfield(HDR.Patient,'Sex')
		if (HDR.Patient.Sex==1)
			fprintf(fid,'\tGender   \t= male\n'); 
		elseif (HDR.Patient.Sex==2)
			fprintf(fid,'\tGender   \t= female\n'); 
		else	
			fprintf(fid,'\tGender   \t= unknown\n');
		end;	
	end;
	T1 = zeros(1,6);
	if isfield(HDR.Patient,'Birthday')
		switch length(HDR.Patient.Birthday)
		case 1,    T1 = datevec(HDR.Patient.Birthday);
		case 6,    T1 = HDR.Patient.Birthday;
		end; 
		if ~any(isnan(T0))
			fprintf(fid,'\tAge\t\t= %4.1f years\n',(datenum(T0)-datenum(T1))/(365.25));
		end;	 
		fprintf(fid,'\tBirthday\t= %04i-%02i-%02i %02i:%02i:%06.3f\n',T1); 
	end;
end;

if isfield(HDR,'Manufacturer')
	fprintf(fid,'Manufacturer.\n');
	if isfield(HDR.Manufacturer,'Name')
		fprintf(fid,'\tName\t\t= %s\n',HDR.Manufacturer.Name); 
	end; 	
	if isfield(HDR.Manufacturer,'Model')
		fprintf(fid,'\tModel\t\t= %s\n',HDR.Manufacturer.Model); 
	end; 	
	if isfield(HDR.Manufacturer,'Version')
		fprintf(fid,'\tVersion \t= %s\n',HDR.Manufacturer.Version); 
	end; 	
	if isfield(HDR.Manufacturer,'SerialNumber')
		fprintf(fid,'\tSerialNumber \t= %s\n',HDR.Manufacturer.SerialNumber); 
	end; 	
end;


%%%%%%%% CHANNEL DATA %%%%%%%%%%%%%%%
if ~isfield(HDR,'AS') && isfield(HDR,'SampleRate')
	HDR.AS.SampleRate = repmat(HDR.SampleRate,HDR.NS,1); 
end;
if ~isfield(HDR.AS,'SPR'),
	HDR.AS.SPR = repmat(HDR.SPR,1,HDR.NS);
end;
if ~isfield(HDR.AS,'SampleRate'),
	HDR.AS.SampleRate = HDR.AS.SPR/HDR.SPR*HDR.SampleRate;  
end; 
if ~isfield(HDR,'THRESHOLD')
	HDR.THRESHOLD = repmat(NaN,HDR.NS,2); 
end;
if ~isfield(HDR,'PhysDimCode') 
	if isfield(HDR,'PhysDim')
		HDR.PhysDimCode = physicalunits(HDR.PhysDim);
	else
		HDR.PhysDimCode = zeros(1,HDR.NS);
	end	 
end;
if ~isfield(HDR,'LeadIdCode')
	HDR = leadidcodexyz(HDR); 
end;
if ~isfield(HDR,'REC')
	HDR.REC.Impedance = repmat(NaN,HDR.NS,1); 
end;
if ~isfield(HDR.REC,'Impedance')
	HDR.REC.Impedance = repmat(NaN,HDR.NS,1); 
end;
if ~isfield(HDR,'InChanSelect')
	InChanSelect = 1:HDR.NS;
else	
	InChanSelect = HDR.InChanSelect;
end
if ~isfield(HDR,'Off')
	HDR.Off = zeros(HDR.NS,1); 
	HDR.Cal(InChanSelect) = diag(HDR.Calib(2:end,:));
end;
if ~isfield(HDR,'Cal') && isfield(HDR,'Calib')
	HDR.Cal = ones(HDR.NS,1); 
	HDR.Cal(InChanSelect) = diag(HDR.Calib(2:end,:));
end;
if HDR.NS,
if length(HDR.Filter.HighPass)==1,
	HDR.Filter.HighPass = repmat(HDR.Filter.HighPass,HDR.NS,1); 
end;
if length(HDR.Cal)==1,
	HDR.Cal = repmat(HDR.Cal,HDR.NS,1); 
end;
if length(HDR.Filter.LowPass)==1,
	HDR.Filter.LowPass = repmat(HDR.Filter.LowPass,HDR.NS,1); 
end;
if length(HDR.Filter.Notch)==1,
	HDR.Filter.Notch = repmat(HDR.Filter.Notch,HDR.NS,1); 
end;
end; 


PhysDim = physicalunits(HDR.PhysDimCode); 
fprintf(fid,'\n[Channel Header]\n#No  LeadId  Label\tfs [Hz]\tGDFTYP\tTH-  TH+  Offset  Calib  PhysDim  HP[Hz]  LP[Hz]  Notch  R[kOhm]  x  y  z\n'); 
for k = 1:HDR.NS,
	Label = HDR.Label{k};
	Z = HDR.REC.Impedance(k)/1000; 
	gdftyp = HDR.GDFTYP(min(length(HDR.GDFTYP),k)); 
	Label(Label==9)=' '; % replace TAB's because TAB's are used as field delimiter
	fprintf(fid,'%3i  %i\t%-9s\t%6.1f %2i  %i\t%i\t%6e\t%6e %5s  %6.4f %5.1f  %i  %5.1f  %f %f %f\n',k,HDR.LeadIdCode(k),Label,HDR.AS.SampleRate(k),gdftyp,HDR.THRESHOLD(k,1:2),HDR.Off(k),HDR.Cal(k),PhysDim{k},HDR.Filter.HighPass(k),HDR.Filter.LowPass(k),HDR.Filter.Notch(k),Z,HDR.ELEC.XYZ(k,:)); 
end;

if ~isfield(HDR.EVENT,'SampleRate');
	HDR.EVENT.SampleRate = HDR.SampleRate;
end;	
%%%%%%%%%% EVENTTABLE %%%%%%%%%%%%%%%
fprintf(fid,'\n[Event Table]\n'); 
fprintf(fid,'NumberOfEvents=%i  SampleRate=%f\n   TYP\t   POS ',length(HDR.EVENT.POS),HDR.EVENT.SampleRate); 
if isfield(HDR.EVENT,'CHN')
	fprintf(fid,'\tCHN\tDUR/VAL'); 
end; 
fprintf(fid,'\tDescription\n'); 

% use global to improve speed
global BIOSIG_GLOBAL;
if ~isfield(BIOSIG_GLOBAL,'ISLOADED_EVENTCODES')
	BIOSIG_GLOBAL.ISLOADED_EVENTCODES = 0;
end; 
if ~BIOSIG_GLOBAL.ISLOADED_EVENTCODES,
	H=sopen('eventcodes.txt'); sclose(H); 
end;

for k = 1:length(HDR.EVENT.POS);
	fprintf(fid,'0x%04x\t%7i',[HDR.EVENT.TYP(k),HDR.EVENT.POS(k)]'); 
	if isfield(HDR.EVENT,'CHN')
		if ~isempty(HDR.EVENT.CHN)
			fprintf(fid,'\t%i\t%i',HDR.EVENT.CHN(k),HDR.EVENT.DUR(k)); 
		end;
	else 
		fprintf(fid,'\t-\t-');	
	end; 
	if HDR.EVENT.TYP(k)==hex2dec('7fff'),
		ch = HDR.EVENT.CHN(k);
		fprintf(fid,'\t%f %s',[1,HDR.EVENT.DUR(k)]*HDR.Calib([1,ch+1],ch),HDR.PhysDim{ch}); 
	elseif HDR.EVENT.TYP(k)==0,
		;
	elseif (isfield(HDR.EVENT,'CodeDesc') && (HDR.EVENT.TYP(k) <= length(HDR.EVENT.CodeDesc)))
		fprintf(fid,'\t%s',HDR.EVENT.CodeDesc{HDR.EVENT.TYP(k)});
	else
		ix = find(HDR.EVENT.TYP(k)==BIOSIG_GLOBAL.EVENT.CodeIndex);
		if length(ix)==1,
			fprintf(fid,'\t%s',BIOSIG_GLOBAL.EVENT.CodeDesc{ix});
		end; 
	end; 
	fprintf(fid,'\n');
end; 

if fid>2,
	fclose(fid); 
end;
if nargout>0,
	argout=HDR;
end;	

