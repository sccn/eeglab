function [DAT,S]=readbdf(DAT,Records,Mode)
% [DAT,signal]=readedf(EDF_Struct,Records)
% Loads selected Records of an EDF File (European Data Format for Biosignals) into MATLAB
% <A HREF="http://www.medfac.leidenuniv.nl/neurology/knf/kemp/edf.htm">About EDF</A> 
%
% Records1	List of Records for Loading
% Mode		0 	default
%		1 	No AutoCalib
%		Mode+2	Concatanated (channels with lower sampling rate if more than 1 record is loaded)

%	Version 2.11
%	03.02.1998
%	Copyright (c) 1997,98 by Alois Schloegl
%	a.schloegl@ieee.org	
                      
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program has been modified from the original version for .EDF files
% The modifications are to the number of bytes read on line 53 (from 2 to
% 3) and to the type of data read - line 54 (from int16 to bit24). Finally the name
% was changed from readedf to readbdf
% T.S. Lorig Sept 6, 2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3 Mode=0; end;
 
EDF=DAT.Head; 
RecLen=max(EDF.SPR);

S=NaN*zeros(RecLen,EDF.NS);
DAT.Record=zeros(length(Records)*RecLen,EDF.NS);
DAT.Valid=uint8(zeros(1,length(Records)*RecLen));
DAT.Idx=Records(:)';
        
for nrec=1:length(Records),

	NREC=(DAT.Idx(nrec)-1);
	if NREC<0 fprintf(2,'Warning READEDF: invalid Record Number %i \n',NREC);end;

	fseek(EDF.FILE.FID,(EDF.HeadLen+NREC*EDF.AS.spb*3),'bof');
	[s, count]=fread(EDF.FILE.FID,EDF.AS.spb,'bit24');

	S(EDF.AS.IDX2)=s;

	%%%%% Test on  Over- (Under-) Flow
%	V=sum([(S'==EDF.DigMax(:,ones(RecLen,1))) + (S'==EDF.DigMin(:,ones(RecLen,1)))])==0;
	V=sum([(S(:,EDF.Chan_Select)'>=EDF.DigMax(EDF.Chan_Select,ones(RecLen,1))) + ...
	       (S(:,EDF.Chan_Select)'<=EDF.DigMin(EDF.Chan_Select,ones(RecLen,1)))])==0;
	EDF.ERROR.DigMinMax_Warning(find(sum([(S'>EDF.DigMax(:,ones(RecLen,1))) + (S'<EDF.DigMin(:,ones(RecLen,1)))]')>0))=1;
	%	invalid=[invalid; find(V==0)+l*k];
                             
	if floor(Mode/2)==1
		for k=1:EDF.NS,
			DAT.Record(nrec*EDF.SPR(k)+(1-EDF.SPR(k):0),k)=S(1:EDF.SPR(k),k);
		end;
	else
		DAT.Record(nrec*RecLen+(1-RecLen:0),:)=S;
	end;

	DAT.Valid(nrec*RecLen+(1-RecLen:0))=V;
end;
if rem(Mode,2)==0	% Autocalib
	DAT.Record=[ones(RecLen*length(Records),1) DAT.Record]*EDF.Calib;
end;                   

DAT.Record=DAT.Record';
return;         

