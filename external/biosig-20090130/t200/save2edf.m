function [HDR] = save2edf(arg1,arg2,arg3);
% SAVE2GDF loads EEG data and saves it in GDF format
%    It has been tested with data of the following formats:
%       Physiobank, BKR, CNT (Neurscan), EDF, 
%
%       HDR = save2edf(sourcefile [, destfile [, option]]);  
%
%	HDR = save2edf(HDR,data);
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
%

% 	$Id: save2edf.m,v 1.1 2009-01-30 06:04:42 arno Exp $
%	Copyright (C) 2003-2005,2008 by Alois Schloegl <a.schloegl@ieee.org>		
%       This file is part of the biosig project http://biosig.sf.net/

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
	if isfield(HDR,'THRESHOLD') 
		HDR.DigMax  = HDR.THRESHOLD(1:HDR.NS,2)';
		HDR.DigMin  = HDR.THRESHOLD(1:HDR.NS,1)';

	else %if ~isfield(HDR,'THRESHOLD')
		fprintf(2,'Warning SAVE2GDF: no THRESHOLD value provided - automated overflow detection not supported\n');

	        HDR.DigMax = max(data,[],1);
		HDR.DigMin = min(data,[],1);
	end; 
	if strcmp(HDR.TYPE,'EVENT')
                HDR.SampleRate = HDR.EVENT.SampleRate;
        elseif isfield(HDR,'GDFTYP')
                HDR.PhysMax = [1,HDR.DigMax]*HDR.Calib;
        	HDR.PhysMin = [1,HDR.DigMin]*HDR.Calib;
    		%bits = ceil(log2(max(HDR.DigMax-HDR.DigMin+1))/8)*8;    % allows int8, int16, int24, int32, etc. 
    		bits1 = ceil(log2(HDR.DigMax-HDR.DigMin+1));
	        [datatyp,limits,datatypes] = gdfdatatype(HDR.GDFTYP);
		bits = log2(limits(:,2)-limits(:,1)+1);
    		fprintf(1,'SAVE2GDF: %i bits needed, %i bits used for file %s\n',max(bits1),max(bits),HDR.FileName);

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
	        if min(dQ)<1,    GDFTYP = 16;  	% float32
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
		% HDR.PhysMax = [1,digmax]*HDR.Calib; %max(data,[],1); %gives max of the whole matrix
	        % HDR.PhysMin = [1,digmin]*HDR.Calib; %min(data,[],1); %gives max of the whole matrix
	        [datatyp,limits,datatypes] = gdfdatatype(HDR.GDFTYP);
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
                if isfield(HDR,'Calib') & ~isfield(HDR,'PhysMax');
                        HDR.PhysMax = [1,HDR.DigMax]*HDR.Calib;
                	HDR.PhysMin = [1,HDR.DigMin]*HDR.Calib;
                end;
	end;

        HDR.FLAG.UCAL = 1;              % data is de-calibrated, no rescaling within SWRITE 
	HDR.TYPE = 'EDF';
        if ~isfield(HDR,'Dur'); 
                HDR.Dur = 1/HDR.SampleRate;
                HDR.SPR = 1; 
        end;

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
if 1,
                HDR = sopen(HDR.FileName,'r');
	        HDR.FLAG.UCAL = 1; 
    		HDR.FLAG.OVERFLOWDETECTION = 0; 
    		[y1,HDR] = sread(HDR,inf);
                HDR = sclose(HDR);
                if all(all((data==y1) | (isnan(data) & isnan(y1)))),
                        fprintf(2,'SAVE2EDF: saving file %s OK.\n',HDR.FileName);
                end;
else
                fprintf(2,'Error SAVE2EDF: saving file ### %s failed\n',HDR.FileName);
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
    		fprintf(2,'Warning SAVE2GDF: EDF+ not fully supported yet.\n'); 
	end;
	end;
	end;
        if ~isfield(HDR,'NS'),
                warning(['number of channels undefined in ',filename]);
                HDR.NS = size(data,2);
        end;
       
        if isempty(outfile), 	% default destination directory  
                ix = max(find(filename=='.'));
                HDR.FileName  = [HDR.FILE.Name,'.edf'];     % destination directory is current working directory 
        elseif isdir(outfile),	% output file
                HDR.FILE.Path = outfile;            
                HDR.FileName  = fullfile(outfile,[HDR.FILE.Name,'.edf']);
        else
                [HDR.FILE.Path,HDR.FILE.Name,Ext] = fileparts(outfile);
                HDR.FileName = fullfile(HDR.FILE.Path,[HDR.FILE.Name,Ext]);
        end;
        HDR=save2edf(HDR,data);
end;
