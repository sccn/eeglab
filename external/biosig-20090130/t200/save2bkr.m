function [HDR] = save2bkr(arg1,arg2,arg3);
% SAVE2BKR loads EEG data and saves it in BKR format
% The following data formats are supported:
%	CNT, EDF, BKR, MAT, etc. format
%
%       HDR = save2bkr(sourcefile [, destfile [, option]]);  
%
%       HDR = eegchkhdr();
%	HDR = save2bkr(HDR,data);
%
%   sourcefile	sourcefile wildcards are allowed
%   destfile	destination file in BKR format 
%	if destfile is empty or a directory, sourcefile but with extension .bkr is used.
%   options
%       gain            Gain factor for unscaled EEG data (e.g. old Matlab files) 
%       'removeDC'      removes mean
%       'regressEOG k:l,m:n'     removes EOG (channels m:n) from EEG (channels k:l)  
%       'autoscale k:l'	uses only channels from k to l for scaling
%       'detrend k:l'	channels from k to l are detrended with an FIR-highpass filter.
%       'PhysMax=XXX'	uses a fixed scaling factor; might be important for concanating BKR files 
%			+XXX and -XXX correspond to the maximum and minimum physical value, resp. 
% 		You can concanate several options by separating with space, komma or semicolon 
%
%   HDR		Header, HDR.FileName must contain target filename
%   data	data samples
%
% Examples: 
%   save2bkr('/tmp/*.cnt',[],'autoscale 5:30');
%	converts all CNT-files from subdir /tmp/ into BKR files 
%       and saves them in the current directory 
%   save2bkr('/tmp/*.cnt','/tmp2/','autoscale 5:30, PhysMax=200');
%	converts all CNT-files from subdir /tmp/ into BKR files 
%       and saves them in the directory /tmp2/
%	
%
%
% see also: EEGCHKHDR, REGRESS_EOG, SLOAD

%	$Revision: 1.1 $
% 	$Id: save2bkr.m,v 1.1 2009-01-30 06:04:42 arno Exp $
%	Copyright (C) 2002-2003 by Alois Schloegl <a.schloegl@ieee.org>		
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


FLAG_REGRESS_EOG = 0;
FLAG_REMOVE_DC = 0;
FLAG_AUTOSCALE = 0;
FLAG_DETREND = 0;
FLAG_PHYSMAX = 0;
FLAG_removeDrift = 0;

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
                [chansel,tmp] = strtok(arg3(tmp+9:length(arg3)),' ;+');
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
                [chansel_dt2,tmp] = strtok(arg3(tmp+11:length(arg3)),' ;+');
                tmp = str2num(chansel_dt2);
                if isempty(tmp),
                        fprintf(2,'invalid RemoveDrift argument %s',chansel_dt2);
                        return;
                else
                        FLAG_removeDrift = 1;
                        chansel_dt2 = tmp;
                end;
        end;
        
        tmp = findstr(lower(arg3),'regresseog');
        if ~isempty(tmp);
                [chansel_dt3,tmp] = strtok(arg3(tmp+11:length(arg3)),' ;,+');
                [chansel_dt4,tmp] = strtok(tmp,' ;,+');
                tmp = str2num(chansel_dt3);
                FLAG_REGRESS_EOG = ~isempty(tmp);
                if isempty(tmp),
                        fprintf(2,'invalid REGRESSEOG argument %s',chansel_dt3);
                        return;
                else
                        FLAG_REGRESS_EOG = 1;
                        chansel_dt3 = tmp;
                end;
                tmp = str2num(chansel_dt4);
                FLAG_REGRESS_EOG = FLAG_REGRESS_EOG * ~isempty(tmp);
                if isempty(tmp),
                        fprintf(2,'invalid REGRESSEOG argument %s',chansel_dt4);
                        return;
                else
                        chansel_dt4 = tmp;
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

if 0,
elseif exist(arg1,'file')
	inpath='';
	infile.name=arg1;
        outfile = arg2;
elseif exist(arg1,'dir')	
        inpath = fileparts(arg1);
        infile = dir(arg1);	% input  file 
        if isempty(infile)
                fprintf(2,'ERROR SAVE2BKR: file %s not found.\n',arg1);
                return;
        end;
        outfile = arg2;
elseif isstruct(arg1) & isnumeric(arg2),
        HDR  = arg1;
        data = arg2;
else  %if isstruct(arg1) & isnumeric(arg2),
        fprintf(2,'Error SAVE2BKR: invalid input arguments\n');	        
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
        %HDR.SampleRate	= 100;		% Sampling rate
        %HDR.Label  	= hdr.Label; % Labels, 
        
        %HDR.PhysMax 	= max(abs(data(:)));	% Physical maximum 
        %HDR.DigMax 	= max(2^15-1);		% Digital maximum
        % --- !!! Previous scaling gave an error up to 6% and more !!! ---
        
        %HDR.Filter.LowPass  = 30;	% upper cutoff frequency
        %HDR.Filter.HighPass = .5;	% lower cutoff frequency
        HDR.FLAG.REFERENCE   = ' ';	% reference '', 'LOC', 'COM', 'LAP', 'WGT'
        %HDR.FLAG.REFERENCE  = HDR.Recording;
        HDR.FLAG.TRIGGERED   = HDR.NRec>1;	% Trigger Flag
        
        if FLAG_REMOVE_DC,
                %y = center(y,1);
                data = data - repmat(mean(data,1),size(data,1),1);
        end;
        if chansel == 0;
                chansel=1:HDR.NS;
        end;
        tmp = data(:,chansel);
        HDR.PhysMax = max(abs(tmp(:))); %gives max of the whole matrix
        HDR.DigMax = 2^15-1;            % maximum resulution
        for k = 1:HDR.NS,
                if any(k==chansel),
                        data(:,k) = data(:,k)*HDR.DigMax/HDR.PhysMax;
                else
                        mm = max(abs(data(:,k)));
                        data(:,k) = data(:,k)*HDR.DigMax/mm;
                end;
        end;
        HDR.FLAG.UCAL = 1;              % data is de-calibrated, no rescaling within SWRITE 
        %HDR = eegchkhdr(HDR);   
        HDR.TYPE = 'BKR';

        HDR = sopen (HDR,'w',0);     	% OPEN BKR FILE
        HDR = swrite(HDR,data);  	% WRITE BKR FILE
        %fwrite(HDR.FILE.FID,data','int16');  	% WRITE BKR FILE
        HDR = sclose(HDR);            % CLOSE BKR FILE
        
        % save Classlabels
        if isfield(HDR,'Classlabel'),
                fid = fopen([HDR.FileName(1:length(HDR.FileName)-4) '.par'],'wt');
                fprintf(fid, '%i\n', HDR.Classlabel);
                fclose(fid);
        end;
        if isfield(HDR,'ArtifactSelection'),
                fid = fopen([HDR.FileName(1:length(HDR.FileName)-4) '.sel'],'w');
                fprintf(fid, '%i\r\n', HDR.ArtifactSelection);
                fclose(fid);
        end;
        
        % final test 
        try
                HDR = sopen(HDR.FileName,'r');
                HDR = sclose(HDR);
        catch
                fprintf(2,'Error SAVE2BKR: saving file %s failed\n',HDR.FileName);
        end;
        return;
end;

for k=1:length(infile);
        filename = fullfile(inpath,infile(k).name);
        [pf,fn,ext] = fileparts(filename);
        
        % load eeg data 
        [y,HDR] = sload(filename);
        
        % load classlabels if the exist
        tmp = fullfile(HDR.FILE.Path,[HDR.FILE.Name,'.mat']);
        if exist(tmp),
                tmp=load(tmp);
                if isfield(tmp,'classlabel');
                        HDR.Classlabel = tmp.classlabel;
                end;
        end;
        
        if isempty(y), 
                fprintf(2,'Error SAVE2BKR: file %s not found\n',filename);
                return; 
        end; 
        
        if ~isfield(HDR,'NS'),
                warning(['number of channels undefined in ',filename]);
                HDR.NS = size(y,2);
        end;
        if ~HDR.FLAG.TRIGGERED,
                HDR.NRec = 1; 
                HDR.SPR = size(y,1);
        end;
        if ~isfield(HDR,'NRec'),
                HDR.NRec = 1;
        end;
        if ~isfield(HDR,'SPR'),
                HDR.SPR = size(y,1)/HDR.NRec;
        elseif length(HDR.SPR)>1,       % use just one sampling rate 
                HDR.SPR = HDR.AS.MAXSPR;
                HDR.SampleRate = HDR.AS.MAXSPR/HDR.Dur;
                FLAG_PHYSMAX = 1; 
                PHYSMAX = max(abs(y(:)));
                HDR.DigMax  = 2^15-1;
        end;

        if FLAG_REGRESS_EOG,
                fprintf(1,'\tREGRESS_EOG \n');
                [R,y] = regress_eog(y,chansel_dt3,chansel_dt4);
        end;

        if FLAG_REMOVE_DC,
                fprintf(1,'\tREMOVE_DC \n');
                y = y - repmat(mean(y,1),size(y,1),1);
        end;
        if FLAG_DETREND,
                B = -ones(1,HDR.SampleRate)/HDR.SampleRate;
                B(HDR.SampleRate/2) = B(HDR.SampleRate/2)+1;
                HDR.Filter.B = B;
                HDR.Filter.A = 1;
                %HDR.Filter.B=B;%conv(-B, HDR.Filter.B);
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
                %B(1) = B(1)+1;
                
                HDR.Filter.B = B;
                HDR.Filter.A = 1;
                %HDR.Filter.B=B;%conv(-B, HDR.Filter.B);
                Delay = (length(B)-1)/2;        
                HDR.FLAG.FILT = 1;
                HDR.Filter.HighPass = .5;
                
                for k = chansel_dt2,
                        y(:,k) = filtfilt(B,1,y(:,k));
                        %y(:,k) = tmp(Delay+1:size(y,1)+Delay);
                end;                
        end;
        
        if chansel == 0;
                chansel=1:HDR.NS;
        end;
        
        % add event channel 
        if isfield(HDR,'EVENT')
                if ~length(HDR.EVENT.TYP),
                elseif 0,
                        % TypeList = unique(HDR.EVENT.TYP); but ignores NaN's
                        [sY ,idx] = sort(HDR.EVENT.TYP(:));
                        TypeList  = sY([1;find(diff(sY,1)>0)+1]);
                        
                        event = zeros(size(y,1),length(TypeList));
                        for k2 = 1:length(TypeList),
                                tmp = (HDR.EVENT.TYP==TypeList(k2));
                                event(HDR.EVENT.POS(tmp),k2) = HDR.EVENT.TYP(tmp);        
                        end;
                        HDR.NS = HDR.NS + size(event,2);
                        y = [y, event];
                elseif all(HDR.EVENT.TYP < 256),  % only NeuroScan Events are converted into separate channels
                        K = 0; event = [];
                        for k1 = 0:7,
                                tmp = bitand(HDR.EVENT.TYP,2^k1);
                                if any(tmp),
                                        K = K+1;
                                        event(size(y,1),K) = 0;        
                                        event(HDR.EVENT.POS(tmp>0),K) = 1;        
                                end;				
                        end;    
                        if any(sum(event,2)>1),
                                fprintf(2,'Warning SAVE2BKR: simulateneous events occur. \n');
                        end;	
                        HDR.NS = HDR.NS + size(event,2);
                        y = [y, event];
                end;
        end;
        
        % re-scale data to account for the scaling factor in the header
        HDR.DigMax = 2^15-1;
        if FLAG_PHYSMAX,
                HDR.PhysMax = PHYSMAX;
        else
                tmp = y(:,chansel);
                HDR.PhysMax = max(abs(tmp(:))); %gives max of the whole matrix
        end;
        for k = 1:HDR.NS,
                if any(k==chansel),
                        y(:,k) = y(:,k)*HDR.DigMax/HDR.PhysMax; % keep correct scaling factor 
                else
                        mm = max(abs(y(:,k)));
                        y(:,k) = y(:,k)*HDR.DigMax/mm;          % scale to maximum resolution
                end;
        end;
        HDR.FLAG.UCAL = 1;              % data is de-calibrated, no rescaling within SWRITE 
        
        tmp = round(HDR.PhysMax);
        fprintf(1,'Rounding of PhysMax yields %f%% error.\n',abs((HDR.PhysMax-tmp)/tmp)*100);
        HDR.PhysMax = tmp;
        HDR.TYPE = 'BKR';
        HDR.FLAG.REFERENCE = ' ';
        HDR.FLAG.TRIGGERED = (HDR.NRec>1);
        
        if isempty(outfile), 	% default destination directory  
                ix = max(find(filename=='.'));
                %HDR.FileName = [filename(1:ix-1),'.bkr'];  % destination directory is same as source directory 
                HDR.FileName  = [HDR.FILE.Name,'.bkr'];     % destination directory is current working directory 
        elseif isdir(outfile),	% output file
                HDR.FILE.Path = outfile;            
                HDR.FileName  = fullfile(outfile,[HDR.FILE.Name,'.bkr']);
        else
                [HDR.FILE.Path,HDR.FILE.Name,Ext] = fileparts(outfile);
                HDR.FileName = fullfile(HDR.FILE.Path,[HDR.FILE.Name,Ext]);
        end;
        %HDR = eegchkhdr(HDR);
        
        HDR = sopen(HDR,'w');
        if HDR.FILE.FID < 0,
                fprintf(1,'Error SAVE2BKR: couldnot open file %s.\n',HDR.FileName);
                return;
        end;
        % writes data
        HDR = swrite(HDR,y(:,1:HDR.NS));  	% WRITE BKR FILE
        %count = fwrite(HDR.FILE.FID,y','short');
        HDR = sclose(HDR);
        
        % save classlabels
        if isfield(HDR,'Classlabel'),
                if ~isempty(HDR.Classlabel),
                        fid = fopen([HDR.FileName(1:length(HDR.FileName)-4) '.par'],'w');
                        fprintf(fid, '%i\r\n', HDR.Classlabel);
                        fclose(fid);
                end;
        end;
        if isfield(HDR,'ArtifactSelection'),
                fid = fopen([HDR.FileName(1:length(HDR.FileName)-4) '.sel'],'w');
                fprintf(fid, '%i\r\n', HDR.ArtifactSelection);
                fclose(fid);
        end;
        
        % final test 
        try
                HDR = sopen(HDR.FileName,'r');
                HDR = sclose(HDR);
        catch
                fprintf(2,'Error SAVE2BKR: saving file %s failed\n',HDR.FileName);
        end;
end;
