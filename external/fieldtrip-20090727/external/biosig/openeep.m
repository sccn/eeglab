function HDR=openeep(HDR,arg2,arg3,arg4,arg5,arg6)
% OPENEEP opens EEProbe files (but does not read the data). 
% However, it is recommended to use SOPEN instead of OPENEEP.
%
% see also: SLOAD, SOPEN, SREAD, SCLOSE, SEOF, STELL, SSEEK.

%	$Id: openeep.m,v 1.1 2009-07-07 02:23:47 arno Exp $
%	Copyright (c) 2007 by Alois Schloegl <a.schloegl@ieee.org>	
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

if strcmp(HDR.TYPE,'EEProbe-CNT'),

       	HDR.FILE.FID = fopen(HDR.FileName,'rb');
	H = openiff(HDR.FILE.FID);		
        
	if isfield(H,'RIFF');
                HDR.RIFF      = H.RIFF;
                HDR.FILE.OPEN = 1; 
                HDR.Label     = {};
                HDR.PhysDim   = {};
                HDR.SPR       = inf; 
                HDR.NRec      = 1; 
                HDR.FILE.POS  = 0; 
                if ~isfield(HDR.RIFF,'CNT');
                	HDR.TYPE = 'unknown'; 
                elseif ~isfield(HDR.RIFF.CNT,'eeph') | ~isfield(HDR.RIFF.CNT,'LIST');
                	HDR.TYPE = 'unknown'; 
                else	
			s = char(HDR.RIFF.CNT.eeph);
			field = '';
			while ~isempty(s)
				[line,s] = strtok(s,[10,13]);
				if strncmp(line,'[Sampling Rate]',15); 
					field = 'SampleRate';
				elseif strncmp(line,'[Samples]',9); 
					field = 'SPR';
				elseif strncmp(line,'[Channels]',10); 
					field = 'NS';
				elseif strncmp(line,'[Basic Channel Data]',20); 
					k = 0; 
					while (k<HDR.NS),     
						[line,s] = strtok(s,[10,13]); 
						if ~strncmp(line,';',1);
							k = k+1;
							[num,status,sa]=str2double(line);
							HDR.Label{k} = sa{1};
							HDR.PhysDim{k,1} = sa{4};
							HDR.Cal(k) = num(2)*num(3);
						end;
					end;	
				elseif strncmp(line,'[History]',9); 
					[t,s] = strtok(s,[10,13]);
					while ~strncmp(t,'EOH',3)
						[t,s] = strtok(s,[10,13]);
					end;
				elseif strncmp(line,';',1);      
				elseif strncmp(line,'[',1);				
					field = '';      
				elseif ~isempty(field);
					[num,status,sa] = str2double(line);
					if ~status,
						HDR = setfield(HDR,field,num);     
						field = '';
					end;	
				end;     
			end;	
			
			% decode data block
			HDR.SPR = H.RIFF.CNT.LIST.raw3.ep(1);
			HDR.NRec = length(H.RIFF.CNT.LIST.raw3.ep)-2;
			
if 0, %try
	%%% ### FIXME ### %%%
	%%% decoding of data block
	%%% this is work in progress. 
			HDR.EEP.epoch_length = H.RIFF.CNT.LIST.raw3.ep(1);
			HDR.EEP.epoch_start  = [H.RIFF.CNT.LIST.raw3.ep(2:end),length(H.RIFF.CNT.LIST.raw3.data)];
			n = [HDR.EEP.epoch_length,length(H.RIFF.CNT.LIST.raw3.chan)];
%dec2bin(double(HDR.RIFF.CNT.LIST.raw3.data(HDR.EEP.epoch_start(1:end-1)+1)),8)
			%n=n([2,1]);
			accu  = zeros(1,HDR.NS); 
			for k = 1:HDR.NRec,
				ix1 = HDR.EEP.epoch_start(k)+1;
				bytes = double(H.RIFF.CNT.LIST.raw3.data(HDR.EEP.epoch_start(k)+1:HDR.EEP.epoch_start(k+1)));
				ix1 = 1;
			for ch= 1:HDR.NS,
				ix1 = ceil(ix1);	
				byte1 = H.RIFF.CNT.LIST.raw3.data(HDR.EEP.epoch_start(k)+ix1);
				meth  = double(bitshift(byte1,-4));

				if any(meth==[1:3])
					nbits = bitand(bytes(1),15);
					nexcbits = bitshift(bytes(2),4); 
					if ~nexcbits, nexcbits = 16; end; 
					y = repmat(NaN,HDR.SPR,1); 
					y(1) = bitshift(bitand(bytes(2),15)*256+bytes(3),4)+bitshift(bytes(4),-4);
					bitpos = 28;
					ix1 = ix1 + 3.5;
					for k1 = 2:HDR.SPR,
						%ix1 = bitpos/8;
						ix2 = ix1 + nbits/8;
						i1 = ceil(ix1);
						r1 = (ix1-i1+1)*8;	
						i2 = ceil(ix2);
						r2 = (ix2-i2+1)*8;	
%[r1,i1,r2,i2],
						a = bitand(bytes(i1),2^r1-1);
						a = bitshift(a,16) + bitshift(bytes(i1+1),8) + bytes(i1+2);
						a = bitshift(a,r2-16-r1);
						ix1 = ix2;
						if a>2^(nbits-1)
							a=a-2^nbits;
						elseif a==2^(nbits-1)
							ix2 = ix1 + nexcbits/8;
							i1 = ceil(ix1);
							r1 = (ix1-i1+1)*8;	
							i2 = ceil(ix2);
							r2 = (ix2-i2+1)*8;	
							a = bitand(bytes(i1),2^r1-1);
							a = bitshift(a,16) + bitshift(bytes(i1+1),8) + bytes(i1+2);
							a = bitshift(a,r2-16-r1);
							ix1 = ix2;
							if a>2^(nexcbits-1)
								a=a-2^nexcbits;
							end;	
						end;
						y(k1)=a;
					end; 
					%fprintf(HDR.FILE.stderr,'Warning SOPEN(EEProbe): decompression method %i not implented\n',meth);
				elseif any(meth==[9:11])
					nbits = bitshift(bitand(bytes(1),15),2) + bitshift(bytes(2),-6);
					nexcbits = bitand(bytes(2),63); 
					y = repmat(NaN,HDR.SPR,1); 
					y(1) = bytes(3:6)*(256.^[3,2,1,0]');
					bitpos = 48; 
					ix1 = ix1 + 6;
					for k1 = 2:HDR.SPR,
						%ix1 = bitpos/8;
						ix2 = ix1 + nbits/8;
						i1 = ceil(ix1);
						r1 = (ix1-i1+1)*8;	
						i2 = ceil(ix2);
						r2 = (ix2-i2+1)*8;	
%[r1,i1,r2,i2],
						a = bitand(bytes(i1),2^r1-1);
						a = bitshift(a,32) + bitshift(bytes(i1+1),24) + bitshift(bytes(i1+2),16) + bitshift(bytes(i1+3),8) + bytes(i1+4);
						a = bitshift(a,r2-32-r1);
						ix1 = ix2;
						if a>2^(nbits-1)
							a=a-2^nbits;
						elseif a==2^(nbits-1)
							ix2 = ix1 + nexcbits/8;
							i1 = ceil(ix1);
							r1 = (ix1-i1+1)*8;	
							i2 = ceil(ix2);
							r2 = (ix2-i2+1)*8;	
							a = bitand(bytes(i1),2^r1-1);
							a = bitshift(a,16) + bitshift(bytes(i1+1),8) + bytes(i1+2);
							a = bitshift(a,r2-16-r1);
							ix1 = ix2;
							if a>2^(nexcbits-1)
								a=a-2^nexcbits;
							end;	
						end;
						y(k1)=a;
					end; 
					fprintf(HDR.FILE.stderr,'Warning SOPEN(EEProbe): decompression method %i not implemented\n',meth);
					error('####');
				elseif meth==0,
					y = bytes(ix1+(1:2:HDR.SPR*2));
					y = y - 256*(y>127);
					y = 256*y + bytes(ix1+(2:2:HDR.SPR*2));
					ix1 = ix1 + 1 + 2 * HDR.SPR;
				elseif meth==8,
					y = bytes(ix1+(1:HDR.SPR*4));
					y(1:4:end) = y(1:4:end) - 256*(y(1:4:end)>127); 
					y = reshape(y,[4,length(y)/4])'*(256.^[3;2;1;0]);
					ix1 = ix1 + 1 + 4 * HDR.SPR;
				else 
					y = repmat(NaN,HDR.SPR,1);
					%fprintf(HDR.FILE.stderr,'ERROR SOPEN(EEProbe): decompression method %i not supported\n',double(meth));
					error(sprintf('EEProbe: decompression method %d not supported\n',double(meth))); 
				end;
				
				if any(meth==[1,9]);
					y = reshape(y,[HDR.SPR,1]);
					y(1,:) = accu(ch)+y(1,:); 
					y = cumsum(y,1); 
					accu(ch) = y(end,:); 
				elseif any(meth==[2,10]);
					y = reshape(y,[HDR.SPR,1]);
					y(2,:) = accu(ch)+y(1,:); 
					y(2:end,:) = cumsum(y(2:end,:),1); 
					accu(ch) = y(end,:); 
				elseif any(meth==[3,11]);
					y = cumsum(reshape(y,[HDR.SPR,1]),2);
				elseif any(meth==[0,8]);
					%y = reshape(y,n([2,1]))';
					y = reshape(y,[HDR.SPR,1]);
				end; 
				HDR.data(HDR.SPR*(k-1)+1:HDR.SPR*k,H.RIFF.CNT.LIST.raw3.chan(ch)+1) = y;
			end;
			end; 		
			HDR.Calib = sparse(2:HDR.NS+1,1:HDR.NS,HDR.Cal); 
			HDR.TYPE = 'native'; 
else %catch
			HDR.Calib = sparse(2:HDR.NS+1,1:HDR.NS,1); % because SREAD uses READ_EEP_CNT.MEX 
end; 			
                end
                HDR.SPR = HDR.SPR*HDR.NRec;
                HDR.NRec = 1; 
        end
end;

	% read event file, if applicable 
	fid = 0; 
	if strcmp(HDR.TYPE,'EEProbe-TRG'),
	        fid = fopen(HDR.FileName,'rt');
        elseif strcmp(HDR.TYPE,'EEProbe-CNT')
	        fid = fopen(fullfile(HDR.FILE.Path,[HDR.FILE.Name,'.trg']),'rt');
	end;
	if fid>0,
                tmp = str2double(fgetl(fid));
                if ~isnan(tmp(1))
                	HDR.EVENT.SampleRate = 1/tmp(1); 
	                N = 0; 
	                while ~feof(fid),
	                        tmp = fscanf(fid, '%f %d %s', 3);
	                        if ~isempty(tmp)
	                                N = N + 1; 
	                                HDR.EVENT.POS(N,1)  = round(tmp(1)*HDR.EVENT.SampleRate);
	                                HDR.EVENT.TYP(N,1)  = 0;
	                                %HDR.EVENT.DUR(N,1) = 0;
	                                %HDR.EVENT.CHN(N,1) = 0;
	                                
	                                HDR.EVENT.TeegType{N,1} = char(tmp(3:end));
	                                HDR.EVENT.TYP(N,1)  = str2double(HDR.EVENT.TeegType{N,1});		% numeric
	                        end
	                end;
	                HDR.EVENT.TYP(isnan(HDR.EVENT.TYP))=0;
	                HDR.TRIG = HDR.EVENT.POS(HDR.EVENT.TYP>0); 
	%                HDR.EVENT.POS = HDR.EVENT.POS(HDR.EVENT.TYP>0);
	%                HDR.EVENT.TYP = HDR.EVENT.TYP(HDR.EVENT.TYP>0);
	%                HDR.EVENT = rmfield(HDR.EVENT,'TeegType');
		end;
                fclose(fid);
        end;
                
	if strcmp(HDR.TYPE,'EEProbe-AVR'),
	        % it appears to be a EEProbe file with an averaged ERP
        	try
        	        tmp = read_eep_avr(HDR.FileName);
        	catch
        	        fprintf(HDR.FILE.stderr,'ERROR SOPEN (EEProbe): Cannot open EEProbe-file, because read_eep_avr.mex not installed. \n');
        	        fprintf(HDR.FILE.stderr,'ERROR SOPEN (EEProbe): see http://www.smi.auc.dk/~roberto/eeprobe/\n');
        	        return;
	        end

        	% convert the header information to BIOSIG standards
	        HDR.FILE.FID = 1;               % ?
        	HDR.FILE.POS = 0;
	        HDR.NS = tmp.nchan;             % number of channels
        	HDR.SampleRate = tmp.rate;      % sampling rate
	        HDR.NRec  = 1;                   % it is an averaged ERP, therefore one record
        	HDR.SPR   = tmp.npnt;             % total number of samples in the file
	        HDR.Dur   = tmp.npnt/tmp.rate;    % total duration in seconds
        	HDR.Calib = [zeros(1,HDR.NS) ; eye(HDR.NS, HDR.NS)];  % is this correct?
	        HDR.Label = char(tmp.label);
	        HDR.PhysDim   = 'uV';
	        HDR.FLAG.UCAL = 1;
	        HDR.FILE.POS  = 0; 
	        HDR.AS.endpos = HDR.SPR;
	        HDR.Label = tmp.label;
	        HDR.TriggerOffset = 0; 
        
	        HDR.EEP.data = tmp.data';
	end;        
