function [HDR] = save2gdf(arg1,arg2,arg3);
% SAVE2GDF loads EEG data and saves it in GDF format
%    It has been tested with data of the following formats:
%       Physiobank, BKR, CNT (Neurscan), EDF, 
%
%       HDR = save2gdf(sourcefile [, destfile [, option]]);  
%
%	HDR = save2gdf(HDR,data);
%
%
% see also: SLOAD, SOPEN, SREAD, SCLOSE, SWRITE

%
%   sourcefile	sourcefile wildcards are allowed
%   destfile	destination file in BKR format 
%	if destfile is empty or a directory, sourcefile but with extension .bkr is used.
%   options
%       gain            Gain factor for unscaled EEG data (e.g. old Matlab files) 
%       'removeDC'      removes mean
%       'autoscale k:l'	uses only channels from k to l for scaling
%       'detrend k:l'	channels from k to l are detrended with an FIR-highpass filter.
%       'PhysMax=XXX'	uses a fixed scaling factor; might be important for concanating BKR files 
%			+XXX and -XXX correspond to the maximum and minimum physical value, resp. 
% 		You can concanate several options by separating with space, komma or semicolon 
%
%   HDR		Header, HDR.FileName must contain target filename
%   data	data samples

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.


% 	$Id: save2gdf.m,v 1.1 2009-01-30 06:04:42 arno Exp $
%	Copyright (C) 2003-2005,2007,2008 by Alois Schloegl <a.schloegl@ieee.org>		
%       This file is part of the biosig project http://biosig.sf.net/


FLAG_REMOVE_DC = 0;
FLAG_AUTOSCALE = 0;
FLAG_DETREND = 0;
FLAG_PHYSMAX = 0;
FLAG_removeDrift = 0;

chansel = 0; 

if nargin<2, arg2=[]; end;

if isstr(arg1), 
        inpath = fileparts(arg1);
        infile = dir(arg1);	% input  file 
        if isempty(infile)
                fprintf(2,'ERROR SAVE2GDF: file %s not found.\n',arg1);
                return;
        end;
        outfile = arg2;
elseif isstruct(arg1) & isnumeric(arg2),
        HDR  = arg1;
        data = arg2;
else  %if isstruct(arg1) & isnumeric(arg2),
        fprintf(2,'Error SAVE2GDF: invalid input arguments\n');	        
        return;
end;		

if isstruct(arg1),
        %HDR.FileName 	= destfile;	% Assign Filename
        if isfield(HDR,'NS')
                if HDR.NS==size(data,2),
                        % It's ok. 
                elseif HDR.NS==size(data,1),
                        warning('data is transposed\n');
                        data = data';
                elseif HDR.NS==size(data,2)+1,
                	HDR.NS = size(data,2); 
                        %warning('data is transposed\n');
                        %data = data';
                else
                        fprintf(2,'HDR.NS=%i is not equal to number of data columns %i\n',HDR.NS,size(data,2));
                        return;
                end;				
        else
                HDR.NS = size(data,2);	% number of channels
        end;

        % THRESHOLD, GDFTYP -> Phys/Dig/Min/Max
	if isa(data,'single') & strncmp(version,'6',1)
		data = double(data);
	end;
        if HDR.NS,
	if isfield(HDR,'THRESHOLD') 
		HDR.DigMax  = HDR.THRESHOLD(1:HDR.NS,2)';
		HDR.DigMin  = HDR.THRESHOLD(1:HDR.NS,1)';

	else %if ~isfield(HDR,'THRESHOLD')
		fprintf(2,'Warning SAVE2GDF: no THRESHOLD value provided - automated overflow detection not supported\n');

	        HDR.DigMax = max(data,[],1);
		HDR.DigMin = min(data,[],1);
	end; 
	end; 

	HDR.TYPE = 'GDF';
	if ~isfield(HDR,'VERSION');
		HDR.VERSION = 2.0;	
	end;	
%	HDR.FLAG.UCAL = 0;              % data is de-calibrated, no rescaling within SWRITE 

	if strcmp(HDR.TYPE,'EVENT')
                HDR.SampleRate = HDR.EVENT.SampleRate;
        elseif isfield(HDR,'GDFTYP')
		if any((HDR.GDFTYP>=16) & (HDR.GDFTYP<=18)) 
			if (HDR.VERSION < 1.90)
				digmin = -(2^30);
				digmax = 2^30;
		                HDR.PhysMax = [1,HDR.DigMax]*HDR.Calib;
		        	HDR.PhysMin = [1,HDR.DigMin]*HDR.Calib;
				data = (data - repmat(HDR.DigMin(:)',size(data,1),1));
				data = data * diag((digmax-digmin)./HDR.Cal) + digmin;
			        HDR.DigMin(:) = digmin; 
			        HDR.DigMax(:) = digmax; 
			        HDR.Cal = (HDR.PhysMax-HDR.PhysMin)./(HDR.DigMax-HDR.DigMin);
	                	HDR.Off = HDR.PhysMin - HDR.Cal .* HDR.DigMin;
	                	HDR.Calib = [HDR.Off;diag(HDR.Cal)];
			else 
				% do nothing
			end
		else 
	                HDR.PhysMax = [1,HDR.DigMax]*HDR.Calib;
	        	HDR.PhysMin = [1,HDR.DigMin]*HDR.Calib;
    			%bits = ceil(log2(max(HDR.DigMax-HDR.DigMin+1))/8)*8;    % allows int8, int16, int24, int32, etc. 
    			bits1 = ceil(log2(HDR.DigMax-HDR.DigMin+1));
		        [datatyp,limits,datatypes] = gdfdatatype(HDR.GDFTYP);
			bits = log2(limits(:,2)-limits(:,1)+1);
	    		fprintf(1,'SAVE2GDF: %i bits needed, %i bits used for file %s\n',max(bits1),max(bits),HDR.FileName);
		end;
		
	        % re-scale data to account for the scaling factor in the header
	        %HIS = histo3(data); save HIS HIS
	else	
		tmp = sort(data,1);
		tmp = diff(tmp);
		tmp(tmp<8*eps) = NaN;
		dQ = min(tmp);
	
		digmax = HDR.DigMax; 
		digmin = HDR.DigMin; 

    		bits = ceil(log2(max(digmax-digmin+1)));        % allows any bit-depth
	        if (min(dQ)<1) & (HDR.VERSION>1.9),	GDFTYP = 16;  	% float32
	        elseif min(dQ)<1, GDFTYP = 5;  	% int32
		elseif bits==8,  GDFTYP = 1;	% int8
	        elseif bits==16, GDFTYP = 3;	% int16
	        elseif bits==32, GDFTYP = 5;	% int32
	        elseif bits==64, GDFTYP = 7;	% int64
	        elseif ~isempty(bits);	GDFTYP = 255+bits;	% intN
		else        	 GDFTYP = 3; 	% int8
	        end;
		HDR.GDFTYP = GDFTYP; 

	        if length(HDR.GDFTYP)==HDR.NS,
    	        elseif length(HDR.GDFTYP)==1,
    	                HDR.GDFTYP = HDR.GDFTYP*ones(1,HDR.NS);  % int16
    	        else
    	                %% PROBLEM 
    	        end

	        [datatyp,limits,datatypes] = gdfdatatype(HDR.GDFTYP);
		if 1, 
			% here, the data is forced to a different data type
			% this is useful if data is float
			% this can cause round-off errors    
		        HDR.DigMin = limits(:,1)'; 
		        HDR.DigMax = limits(:,2)'; 
		        HDR.FLAG.UCAL = 0;   % data is calibrated, rescaling within SWRITE 
	
		else
			% here, the data is of integer type
			% no round of errors occur. 
			
		        c0 = 0; 
		        while any(digmin'<limits(:,1)),
	        	        c = 2^ceil(log2(max(limits(:,1)-digmin')))
	                	digmin = digmin + c;
	        	        digmax = digmax + c;
		                c0 = c0 + c;
		        end;
	        	while any(digmax'>limits(:,2)),
		                c = 2^ceil(log2(max(digmax'-limits(:,2))))
		                digmin = digmin - c;
	                	digmax = digmax - c;
		                c0 = c0 - c;
	        	end;
		        while any(digmin'<limits(:,1)),
		                c = 2^ceil(log2(max(limits(:,1)-digmin')))
	                	digmin = digmin + c;
		                digmax = digmax + c;
	        	        c0 = c0 + c;
		        end;

                	data = data + c0;
		        
		        HDR.DigMax = digmax; %limits(:,2); %*ones(1,HDR.NS);
		        HDR.DigMin = digmin; %limits(:,1); %*ones(1,HDR.NS);
	        
		        %fprintf(1,'Warning SAVE2GDF: overflow detection not implemented, yet.\n');
	        end; 
                if isfield(HDR,'Calib') & ~isfield(HDR,'PhysMax');
               	        HDR.PhysMax = [1,HDR.DigMax]*HDR.Calib;
       	        	HDR.PhysMin = [1,HDR.DigMin]*HDR.Calib;
                end;
	end;

        if ~isfield(HDR,'Dur'); 
                HDR.Dur = 1/HDR.SampleRate;
                HDR.SPR = 1; 
        end;

%%	[HDR.PhysMax;HDR.PhysMin;HDR.DigMax;HDR.DigMin;max(data);min(data)],

        HDR = sopen(HDR,'w');
        if HDR.FILE.FID < 0,
                fprintf(1,'Error SAVE2GDF: couldnot open file %s.\n',HDR.FileName);
                return;
        end;

        if numel(data)>0,
	        HDR = swrite(HDR,data(:,1:HDR.NS));  	% WRITE GDF FILE
        end;
	HDR = sclose(HDR);

        % final test 
        try
                H2 = sopen(HDR.FileName,'r');
	        H2.FLAG.UCAL = 0; 
    		H2.FLAG.OVERFLOWDETECTION = 0; 
    		[y1,H2] = sread(H2,inf);
                H2 = sclose(H2);
		d2 = [ones(size(data,1),1),data]*H2.Calib;
                if all(all((d2==y1) | (isnan(d2) & isnan(y1)))),
                        fprintf(2,'SAVE2GDF: saving file %s OK.\n',HDR.FileName);
		else 
                        fprintf(2,'SAVE2GDF: file %s saved. Maximum relative roundoff error is %f.\n',HDR.FileName, max(max((d2-y1)./(abs(d2)+abs(y1)))) );
                end;
        catch
                fprintf(2,'Error SAVE2GDF: saving file %s failed\n',HDR.FileName);
        end;
        return;
end;

for k=1:length(infile);
        filename = fullfile(inpath,infile(k).name);
        [pf,fn,ext] = fileparts(filename);
        
        % load eeg data 
        %[data,HDR] = sload(filename);
        HDR = sopen(filename,'r',0);
        if HDR.FILE.FID<0, 
                fprintf(2,'Error SAVE2GDF: file %s not found\n',filename);
                return; 
        end; 
        HDR.FLAG.UCAL = 1; 
        HDR.FLAG.OVERFLOWDETECTION = 0; 
        [data,HDR] = sread(HDR,inf);

        HDR = sclose(HDR);
	if isfield(HDR,'EDF')
	if isfield(HDR.EDF,'Annotations')        
	if ~isempty(HDR.EDF.Annotations)
    		fprintf(2,'Warning SAVE2GDF: Annotations in EDF+ are not fully supported.\n'); 
	end;
	end;
	end;
        if ~isfield(HDR,'DigMax'),
	        HDR.DigMax = max(data,[],1);
	end;          
        if ~isfield(HDR,'DigMin'),
	        HDR.DigMin = min(data,[],1); 
	end;          
        if ~isfield(HDR,'DigMin'),
		HDR.PhysMax = [1,HDR.DigMax]*HDR.Calib;
	end;          
        if ~isfield(HDR,'DigMin'),
		HDR.PhysMin = [1,HDR.DigMin]*HDR.Calib;
	end; 	        	
        if ~isfield(HDR,'NS'),
                warning(['number of channels undefined in ',filename]);
                HDR.NS = size(data,2);
        end;
       
        if isempty(outfile), 	% default destination directory  
                ix = max(find(filename=='.'));
                HDR.FileName  = [HDR.FILE.Name,'.gdf'];     % destination directory is current working directory 
        elseif isdir(outfile),	% output file
                HDR.FILE.Path = outfile;            
                HDR.FileName  = fullfile(outfile,[HDR.FILE.Name,'.gdf']);
        else
                [HDR.FILE.Path,HDR.FILE.Name,Ext] = fileparts(outfile);
                HDR.FileName = fullfile(HDR.FILE.Path,[HDR.FILE.Name,Ext]);
        end;

        HDR=save2gdf(HDR,data);
end;
