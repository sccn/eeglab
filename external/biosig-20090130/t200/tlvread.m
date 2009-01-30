function [TLV] = tlvread(fid)
% FEFOPEN opens and reads FEF file 
%
%       HDR=fefopen(HDR)
%
%
% see also: SOPEN 
%
% References: 
% [1] <A HREF="ftp://sigftp.cs.tut.fi/pub/eeg-data/standards/cenf060.zip ">About CEN/TC251</A> 
% [2] http://telecom.htwm.de/ASN1/ber.htm
% [3] http://asn1.elibel.tm.fr

%	$Id: tlvread.m,v 1.1 2009-01-30 06:04:42 arno Exp $
%	Copyright (C) 2004  Alois Schloegl <a.schloegl@ieee.org>	
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

if feof(fid), TLV =[]; return; end;

tag=fread(fid,1,'uint8');
if bitand(tag,31)==31,
        t  = fread(fid,1,'uint8');
	ac = bitand(t,127);
        while bitand(t,128),
                t  = fread(fid,1,'uint8');
		ac = [ac*128 + bitand(t,127)];
        end;
        TAG = [tag,ac];
else
        TAG = tag;
end;

LEN = fread(fid,1,'uint8');
if bitand(LEN,128),
        len = bitand(LEN,127);
        if len==0,
                % unbestimmte form
                LEN = NaN; 
        else
                % Langform 
                len = fread(fid,len,'uint8');
                LEN = sum(len.*256.^[length(len)-1:-1:0]');
        end
% else % Kurzform
end;

class = bitand(tag,192)/64;             % get class 
FLAG.Primitive = ~bitand(tag,32);  % check P/C bit 


if (tag == 0) & (LEN==0)
        TLV = []; 
        return;
end;

TLV.TAG = TAG; 
TLV.Class = class; 
TLV.PC = FLAG.Primitive;
TLV.LEN = LEN;

if FLAG.Primitive,
        if ~isnan(LEN),
                if class == 0,  % Universal 
			classtag = bitand(tag,31);
			if classtag==0, 
	    		        TLV.reservedBER = fread(fid,LEN,'uchar');
			elseif classtag==1, 
	    		        TLV.Boolean = fread(fid,LEN,'uchar');
			elseif classtag==2, 
	    		        TLV.Integer = fread(fid,LEN,'int8');
			elseif classtag==3, 
	    		        TLV.Bitstring = fread(fid,LEN,'uchar');
			elseif classtag==4, 
	    		        TLV.OctetString = dec2hex(fread(fid,LEN,'uchar'));
			elseif classtag==5, 
	    		        TLV.Null = fread(fid,LEN,'uchar');
			elseif classtag==6, 
	    		        TLV.OID = fread(fid,LEN,'uchar');
			elseif classtag==7, 
	    		        TLV.ObjectDesc = fread(fid,LEN,'uchar');
			elseif classtag==8, 
	    		        TLV.external = fread(fid,LEN,'uchar');
			elseif classtag==9, 
	    		        TLV.real = fread(fid,LEN,'uchar');
			elseif classtag==10, 
	    		        TLV.enum = fread(fid,LEN,'uchar');
			elseif classtag==11, 
	    		        TLV.UTF8String = fread(fid,LEN,'uchar');
			elseif classtag==12, 
	    		        TLV.relativeOID = fread(fid,LEN,'uchar');
			elseif classtag==13, 
	    		        TLV.OID = fread(fid,LEN,'uchar');
                                
			elseif classtag==14, 
	    		        TLV.reserved14 = fread(fid,LEN,'uchar');
			elseif classtag==15, 
	    		        TLV.reserved15 = fread(fid,LEN,'uchar');
			elseif classtag==16, 
	    		        TLV.sequenceof = fread(fid,LEN,'uchar');
			elseif classtag==17, 
	    		        TLV.setof = fread(fid,LEN,'uchar');
			elseif classtag==18, 
	    		        TLV.numericString = fread(fid,LEN,'uchar');
			elseif classtag==19, 
	    		        TLV.PrintableString = fread(fid,LEN,'uchar');
			elseif classtag==20, 
	    		        TLV.TeletexStringT61 = fread(fid,LEN,'uchar');
			elseif classtag==21, 
	    		        TLV.VALUE = fread(fid,LEN,'uchar');
			elseif classtag==22, 
	    		        TLV.IA5String = fread(fid,LEN,'uchar');
			elseif classtag==23, 
	    		        TLV.UTCtime = fread(fid,LEN,'uchar');
			elseif classtag==24, 
	    		        TLV.generalizedTime = fread(fid,LEN,'uchar');
			elseif classtag==25, 
	    		        TLV.GraphicString = fread(fid,LEN,'uchar');
			elseif classtag==26, 
	    		        TLV.VisibleStringISO646String = fread(fid,LEN,'uchar');
			elseif classtag==27, 
	    		        TLV.GeneralString = fread(fid,LEN,'uchar');
			elseif classtag==28, 
	    		        TLV.UniversalString = fread(fid,LEN,'uchar');
			elseif classtag==29, 
	    		        TLV.CharacterString = fread(fid,LEN,'uchar');
			elseif classtag==30, 
	    		        TLV.BMPString = fread(fid,LEN,'uchar');
			elseif classtag==31, 
	    		        TLV.VALUE = fread(fid,LEN,'uchar');
			else
	    		        TLV.VALUE = fread(fid,LEN,'uchar');
			end;	
                        
                elseif class == 1,  % Application
    		        %TLV.VALUE = fread(fid,LEN,'uchar');
	                [TLV.VALUE] = tlvread(fid);
                        
                elseif class == 2,  % Context-Specific
            		TLV.VALUE = fread(fid,LEN,'uchar');
                        
                elseif class == 3,  % Private
            		TLV.VALUE = fread(fid,LEN,'uchar');
                        
                end;
        else
                warning('unspecified length %i',ftell(fid));
                accu = 255;
                tmp = fread(fid,256,'uint8');
                TLV.VALUE = [];
                while ~any((tmp==0) & ([accu;tmp(1:end-1)]==0)) & ~feof(fid)
                        TLV.VALUE = [TLV.VALUE,tmp'];
                        accu = tmp(end);
                        tmp = fread(fid,256,'uint8');
                end;
                x = find((tmp==0) & ([accu;tmp(1:end-1)]==0)); 
                TLV.VALUE = [TLV.VALUE, tmp(end-min(x))'];
                status = fseek(fid,min(x)-length(tmp),0);
        end;
        
else
        K = 0; 
        VAL = tlvread(fid);
        while ~isempty(VAL), 
                K = K + 1; 
                TLV.VALUE{K} = VAL; 
                VAL = tlvread(fid);
        end;
end
