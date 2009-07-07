function [HDR]=swrite(HDR,data)
% SWRITE writes signal data. 
% HDR = swrite(HDR,data)
% Appends data to an Signal File 

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
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


%	$Id: swrite.m,v 1.1 2009-07-07 02:23:48 arno Exp $
%	Copyright (c) 1997-2005 by Alois Schloegl <a.schloegl@ieee.org>	
%       This file is part of the biosig project http://biosig.sf.net/


if HDR.FILE.OPEN < 2,
	fprintf(HDR.FILE.stderr,'Error SWRITE can not be applied, File %s is not opened in WRITE mode\n',HDR.FileName);
	return;
end;


if strcmp(HDR.TYPE,'EDF') | strcmp(HDR.TYPE,'GDF') | strcmp(HDR.TYPE,'BDF'),
        if ~all(HDR.GDFTYP==HDR.GDFTYP(1)) 
                fprintf(2,'Error SWRITE: different GDFTYPs not supported yet!\n');
                return;
        end;        
        if ~HDR.FLAG.UCAL,
	   	data = data - repmat(HDR.PhysMin(:)',size(data,1),1);
	   	data = data * sparse(1:HDR.NS,1:HDR.NS,(HDR.DigMax-HDR.DigMin)./(HDR.PhysMax-HDR.PhysMin));  % scale Phys->Dig
	   	data = data + repmat(HDR.DigMin(:)',size(data,1),1);
        end;

        if ~any(HDR.GDFTYP(1)==[16,17,18]),
                data(data< HDR.THRESHOLD(1,1)) = HDR.THRESHOLD(1,1); %underflow
                data(data> HDR.THRESHOLD(1,2)) = HDR.THRESHOLD(1,2); %overflow 
                data(isnan(data))  = HDR.THRESHOLD(1,3);        % missing value
        end
        
        if 0, HDR.SIE.RAW,
                if sum(HDR.AS.SPR)~=size(data,1)
                        fprintf(2,'Warning SWRITE: datasize must fit to the Headerinfo %i %i %i\n',HDR.AS.spb,size(data));
                        fprintf(2,'Define the Headerinformation correctly.\n',HDR.AS.spb,size(data));
                end;
                D = data; 
                
        elseif (HDR.SPR == 1), 
                D = data'; 
                
        else    
                % fill missing data with NaN
                tmp = rem(size(data,1),HDR.SPR);
		if tmp,
			fprintf(HDR.FILE.stderr,'Warning SWRITE: %i NaNs added to complete data block.\n',HDR.SPR-tmp);
	                data = [data;repmat(HDR.THRESHOLD(1,3),HDR.SPR-tmp,size(data,2))];
		end;	
                NRec = size(data,1)/HDR.SPR;
                D = repmat(NaN,sum(HDR.AS.SPR),NRec);
                for k = 1:HDR.NS;
                        if HDR.AS.SPR(k)>0,
                                D(HDR.AS.bi(k)+1:HDR.AS.bi(k+1),:) = rs(reshape(data(:,k),HDR.SPR,NRec),HDR.SPR/HDR.AS.SPR(k),1);
                        end;
                end;
        end;
	GDFTYP = HDR.GDFTYP(1);
	if ~exist('OCTAVE_VERSION','builtin')
	        count = fwrite(HDR.FILE.FID,D,gdfdatatype(GDFTYP));
	else
		if GDFTYP<256,
	    		count = fwrite(HDR.FILE.FID,D,gdfdatatype(GDFTYP));
		else
		    [datatyp,limits,datatypes] = gdfdatatype(GDFTYP);
		    [nr,nc] = size(D);
		    if 1, %GDFTYP>511, % unsigned bitN
			bits = GDFTYP - 511; 
			
			if 0, 
			elseif bits==4,
				X(2:2:2*nr,:) = floor(D/16);
				X(1:2:2*nr,:) = mod(D,16);
		    		count = fwrite(HDR.FILE.FID,X,'uint8');
			elseif bits==12,
				X(3:3:1.5*nr,:) = floor(D(2:2:nr,:)/16);
				X(1:3:1.5*nr,:) = mod(D([1:2:nr],:),256);
				X(2:3:1.5*nr,:) = mod(floor(D(1:2:nr,:)/256),16)+16*mod(D(2:2:nr,:),16);
		    		count = fwrite(HDR.FILE.FID,X,'uint8');
			elseif bits==24,
				X(3:3:3*nr,:) = mod(floor(D*2^-16),256);
				X(1:3:3*nr,:) = mod(D,256); 
				X(2:3:3*nr,:) = mod(floor(D/256),256);
		    		count = fwrite(HDR.FILE.FID,X,'uint8');
			else
				ix0 = ceil([1:nr*bits]'/8);
				ix1 =  mod([1:nr*bits]',8);
				ix2 = ceil([1:nr*bits]'/bits);
				ix3 =  mod([1:nr*bits]',bits);
				for k = 1:nr*bits,
				
				end;
			end;	
		    end;
		end;
	end;	
        HDR.FILE.POS  = HDR.FILE.POS  + count/HDR.AS.spb;
        
        
elseif strcmp(HDR.TYPE,'BKR'),
        count=0;
        if HDR.NS~=size(data,2) & HDR.NS==size(data,1),
                fprintf(2,'SWRITE: number of channels fits number of rows. Transposed data\n');
                data = data';
        end
        
        if any(HDR.NS==[size(data,2),0]),	% check if HDR.NS = 0 (unknown channel number) or number_of_columns 
                if HDR.NS==0, HDR.NS = size(data,2); end;   % if HDR.NS not set, set it. 
                if ~HDR.FLAG.UCAL,
                        if isnan(HDR.PhysMax)
                                HDR.PhysMax = max(abs(data(:)));
                        elseif HDR.PhysMax < max(abs(data(:))),
                                fprintf(2,'Warning SWRITE: Data Saturation. max(data)=%f is larger than HDR.PhysMax %f.\n',max(abs(data(:))),HDR.PhysMax);
                        end;
                        data = data*(HDR.DigMax/HDR.PhysMax);
                end;
                % Overflow detection
                data((data>2^15-1) | (data<-2^15)) = -2^15; 
                count = fwrite(HDR.FILE.FID,data','short');
        else
                fprintf(2,'Error SWRITE: number of columns (%i) does not fit Header information (number of channels HDR.NS %i)',size(data,2),HDR.NS);
                return;
        end;
        if HDR.SPR==0, 
                HDR.SPR=size(data,1);
        end;
        if HDR.FLAG.TRIGGERED > 1,
                HDR.FILE.POS = HDR.FILE.POS + size(data,1)/HDR.SPR;
	else
		HDR.FILE.POS = 1;		% untriggered data
	end;
        %HDR.AS.endpos = HDR.AS.endpos + size(data,1);
        
        
elseif strcmp(HDR.TYPE,'CFWB')
        count=0;
        if HDR.NS~=size(data,2) & HDR.NS==size(data,1),
                fprintf(2,'SWRITE: number of channels fits number of rows. Transposed data\n');
                data = data';
        end
        
	if HDR.GDFTYP==3,
		if ~HDR.FLAG.UCAL,
			data = data*diag(1./HDR.Cal)-HDR.Off(:,ones(1,size(data,1)))';
		end;
                % Overflow detection
                data(data>2^15-1)=  2^15-1;	
                data(data<-2^15) = -2^15;
                count = fwrite(HDR.FILE.FID,data','short');
        else
                count = fwrite(HDR.FILE.FID,data',gdfdatatype(HDR.GDFTYP));
        end;
        if HDR.NRec==0, 
                HDR.NRec=size(data,1);
        end;
        HDR.FILE.POS = HDR.FILE.POS + size(data,1);
        %HDR.AS.endpos = HDR.AS.endpos + size(data,1);
        
        
elseif strcmp(HDR.TYPE,'AIF') | strcmp(HDR.TYPE,'SND') | strcmp(HDR.TYPE,'WAV'),
        count = 0;
        if (HDR.NS ~= size(data,2)) & (HDR.NS==size(data,1)),
                fprintf(2,'Warning SWRITE: number of channels fits number of rows. Transposed data\n');
                data = data';
        end
        
        if strcmp(HDR.TYPE,'SND') 
                if (HDR.FILE.TYPE==1),
                        data = lin2mu(data);
                end;
        elseif strcmp(HDR.TYPE,'WAV'),
                if ~HDR.FLAG.UCAL,
                        data = round((data + HDR.Off) / HDR.Cal-.5);
                end;
	elseif strcmp(HDR.TYPE,'AIF'),
		if ~HDR.FLAG.UCAL,
                        data = data * 2^(HDR.Bits-1);
		end;
	end;

        count = fwrite(HDR.FILE.FID,data',gdfdatatype(HDR.GDFTYP));
	if HDR.NS==0, 
                HDR.NS=size(data,2);
        end;
        
        
elseif strcmp(HDR.TYPE,'MIT')
        count = 0;
        if HDR.NS~=size(data,2) & HDR.NS==size(data,1),
                fprintf(2,'SWRITE: number of channels do not fit number of columns but number of rows. Data transposed!?!\n');
                data = data';
        end
        
	if ~HDR.FLAG.UCAL,
		data = data - HDR.Off(ones(size(data,1),1),:);
		data = data * diag(1./HDR.Cal);
	end;
	if all(HDR.GDFTYP==3),
                % Overflow detection
                data(data>2^15-1)=  2^15-1;	
                data(data<-2^15) = -2^15;
                count = fwrite(HDR.FILE.FID,data','short');
        else
                fprintf(2,'ERROR SWRITE (MIT): datatype %i not supported\n',HDR.GDFTYP(1));
        end;
        if HDR.NRec==0, 
                HDR.NRec=size(data,1);
        end;
        HDR.FILE.POS = HDR.FILE.POS + size(data,1);
        
        
elseif strcmp(HDR.TYPE,'ET-MEG')
	count = fwrite(HDR.FILE.FID,data',gdfdatatype(HDR.GDFTYP));
        HDR.FILE.POS = HDR.FILE.POS + size(data,1);
        if ~(HDR.NS>0)
		HDR.NS = size(data,2);
		HDR.FILE.OPEN = 3;
	end;		
        if (HDR.FILE.POS~=HDR.NRec*HDR.SPR)
		HDR.SPR = size(data,1);
		HDR.FILE.OPEN = 3;
	end;		
        
	        

else
        fprintf(2,'Error SWRITE: file type %s not supported \n',HDR.TYPE);
        

end;                        
