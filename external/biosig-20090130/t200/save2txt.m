function [HDR] = save2txt(arg1,arg2,arg3);
% SAVE2TXT loads biosignal data and saves it in ASCII-Text format
% The following data formats are supported:
%	SCP, CNT, EDF, BKR, MAT, etc. format
%
%       HDR = save2txt(sourcefile [, destfile [, option]]);  
%
%   sourcefile	sourcefile wildcards are allowed
%   destfile	destination file in TXT format 
%	if destfile is empty or a directory, sourcefile but with extension .txt is used.
%
% see also: SAVE2BKR
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

%	$Revision: 1.1 $
% 	$Id: save2txt.m,v 1.1 2009-01-30 06:04:42 arno Exp $
%	Copyright (C) 2004 by Alois Schloegl <a.schloegl@ieee.org>		
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


FLAG_REMOVE_DC = 0;
FLAG_AUTOSCALE = 0;
FLAG_DETREND = 0;
FLAG_PHYSMAX = 0;
FLAG_removeDrift = 0;
FORM = '%14g\t';

chansel = 0; 

if nargin<2, arg2=[]; end;

if nargin<3,
        cali = NaN;
elseif isnumeric(arg3)
        cali = arg3;        
else
        FLAG_REMOVE_DC = findstr(lower(arg3),'removedc');        
        tmp = findstr(arg3,'autoscale');
        if ~isempty(tmp);
                [chansel,tmp] = strtok(arg3(tmp+9:length(arg3)),' ;,+');
                tmp = str2num(chansel);
                if isempty(tmp),
                        fprintf(2,'invalid autoscale argument %s',chansel);
                        return;
                else
		        FLAG_AUTOSCALE = 1;
	                chansel = tmp;
                end;
        end;
        tmp = findstr(lower(arg3),'detrend');
        if ~isempty(tmp);
                [chansel_dt,tmp] = strtok(arg3(tmp+7:length(arg3)),' ;,+');
                tmp = str2num(chansel_dt);
                if isempty(tmp),
                        fprintf(2,'invalid detrend argument %s',chansel_dt);
                        return;
                else
                        FLAG_DETREND = 1;
	                chansel_dt = tmp;
                end;
        end;
        
        tmp = findstr(lower(arg3),'removedrift');
        if ~isempty(tmp);
                [chansel_dt2,tmp] = strtok(arg3(tmp+11:length(arg3)),' ;,+');
                tmp = str2num(chansel_dt2);
                if isempty(tmp),
                        fprintf(2,'invalid RemoveDrift argument %s',chansel_dt2);
                        return;
                else
                        FLAG_removeDrift = 1;
	                chansel_dt2 = tmp;
                end;
        end;
        tmp = findstr(lower(arg3),'physmax=');
        if ~isempty(tmp);
                [tmp,tmp1] = strtok(arg3(tmp+8:length(arg3)),' ;,');
                PHYSMAX = str2num(tmp);
                if isempty(PHYSMAX ),
                        fprintf(2,'invalid PhysMax argument %s',tmp);
                        return;
                else
                        FLAG_PHYSMAX = 1;
                end;
        end;
end;


if isstr(arg1), 
        inpath = fileparts(arg1);
        infile.name = arg1;%dir(arg1);	% input  file 
        if isempty(infile)
                fprintf(2,'ERROR SAVE2TXT: file %s not found.\n',arg1);
                return;
        end;
        outfile = arg2;
elseif isstruct(arg1) & isnumeric(arg2),
	HDR  = arg1;
	data = arg2;
else  %if isstruct(arg1) & isnumeric(arg2),
        fprintf(2,'Error SAVE2TXT: invalid input arguments\n');	        
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
		else
			fprintf(2,'HDR.NS=%i is not equal number of data columns %i\n',HDR.NS,size(data,2));
			return;
		end;				
	else
	        HDR.NS = size(data,2);	% number of channels
	end;
        
        if ~isfield(HDR,'NRec'),
                HDR.NRec = 1;		% number of trials (1 for continous data)
        end;	
        HDR.SPR 	= size(data,1)/HDR.NRec;	% number of samples per trial
        
        if FLAG_REMOVE_DC,
                %y = center(y,1);
                data = data - repmat(mean(data,1),size(data,1),1);
        end;
        if chansel == 0;
                chansel=1:HDR.NS;
        end;
        tmp = data(:,chansel);
        HDR.PhysMax = max(abs(tmp(:))); %gives max of the whole matrix
        for k = 1:HDR.NS,
                if any(k==chansel),
                        data(:,k) = data(:,k)*HDR.DigMax/HDR.PhysMax;
                else
                        mm = max(abs(data(:,k)));
                        data(:,k) = data(:,k)*HDR.DigMax/mm;
                end;
        end;
        
        format = '';
        for k1 = 1:HDR.NS,
                format = [format, FORM];
        end
        format = [format,'\n'];
        
        fid   = fopen(HDR.FileName,'w+');
        if fid < 0,
                fprintf('Error SAVE2TXT: couldnot open file %s.\n',HDR.FileName);
                return;
        end;
        count = fprintf(fid,format,y');
        fclose(fid);
        
        return;
end;

for k = 1:length(infile);
        filename = fullfile(inpath,infile(k).name);
        [pf,fn,ext] = fileparts(filename);
        
        % load eeg data 
        [y,HDR] = sload(filename);
        
        if isempty(y), 
                fprintf(2,'Error SAVE2TXT: file %s not found\n',filename);
                return; 
        end; 
        
        if ~isfield(HDR,'NS'),
                warning(['number of channels undefined in ',filename]);
                HDR.NS = size(y,2);
        end;
        if ~isfield(HDR,'NRec'),
                HDR.NRec = 1;
        end;
        if ~isfield(HDR,'SPR'),
                HDR.SPR = size(y,1)/HDR.NRec;
        end;
        
        if FLAG_REMOVE_DC,
                y = y - repmat(mean(y,1),size(y,1),1);
        end;
        if FLAG_DETREND,
                B = -ones(1,HDR.SampleRate)/HDR.SampleRate;
                B(HDR.SampleRate/2) = B(HDR.SampleRate/2)+1;
                HDR.Filter.B = B;
                HDR.Filter.A = 1;
                
                Delay = length(B)/2;        
                HDR.FLAG.FILT = 1;
                HDR.Filter.HighPass = .5;
                
                for k = chansel_dt,
                        tmp = filter(B,1,[y(:,k);zeros(length(B),1)]);
                        y(:,k) = tmp(Delay+1:size(y,1)+Delay);
	        end;                
        end;
        
        if FLAG_removeDrift,
                B = .5*(1 - cos(2*pi*(1:4*HDR.SampleRate+1)'/(4*HDR.SampleRate+2))); 
                B = -B/sum(B);
                B(2*HDR.SampleRate) = B(HDR.SampleRate)+1;
                
                B = -ones(1,HDR.SampleRate)/HDR.SampleRate;
                B(HDR.SampleRate/2) = B(HDR.SampleRate/2)+1;
                
                HDR.Filter.B = B;
                HDR.Filter.A = 1;
                
                Delay = (length(B)-1)/2;        
                HDR.FLAG.FILT = 1;
                HDR.Filter.HighPass = .5;
                
                for k = chansel_dt2,
                        y(:,k) = filtfilt(B,1,y(:,k));
	        end;                
        end;
        
        if chansel == 0;
                chansel=1:HDR.NS;
        end;

        % add event channel 
        if isfield(HDR,'EVENT')
                if HDR.EVENT.N > 0,
                        event = zeros(size(y,1),1);
                        event(HDR.EVENT.POS) = HDR.EVENT.TYP;        
                        HDR.NS = HDR.NS + 1;
                        y = [y, event];
                end;
        end;
        
        HDR.FILE.Ext = 'txt';
        
        if isempty(outfile), 	% default destination directory  
                ix = max(find(filename=='.'));
                %HDR.FileName = [filename(1:ix-1),'.bkr'];  % destination directory is same as source directory 
                HDR.FileName  = [HDR.FILE.Name,'.txt'];     % destination directory is current working directory 
        elseif 0,isdir(outfile),	% output file
                HDR.FILE.Path = outfile;            
	        HDR.FileName  = fullfile(outfile,[HDR.FILE.Name,'.txt']);
        else
                %[HDR.FILE.Path,HDR.FILE.Path,HDR.FILE.Path]=fileparts(outfile);
                HDR.FileName = outfile;
        end;
        
        format = '';
        for k1 = 1:HDR.NS,
                format = [format, FORM];
        end
        format = [format,'\n'];
        
        fid   = fopen(HDR.FileName,'w+');
        if fid < 0,
                fprintf('Error SAVE2TXT: couldnot open file %s.\n',HDR.FileName);
                return;
        end;
                
        count = fprintf(fid,format,y');
        fclose(fid);
        
end;
