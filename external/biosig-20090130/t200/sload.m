function [signal,H] = sload(FILENAME,varargin)
% SLOAD loads signal data of various data formats
% [signal,header] = sload(FILENAME)
% [signal,header] = sload(FILENAME,CHAN)
%       read selected channels in list CHAN
%       CHAN=0 [default], reads all channels
%
% [signal,header] = sload(dir('f*.emg') ... )
% [signal,header] = sload('f*.emg' ...)
%  	loads all files 'f*.emg'
%
% [signal,header] = sload(..., PropertyName1,PropertyValue1,...)
%  	PropertyName(s)		PropertyValue
%	'UCAL'			-		data uncalibrated (not scaled)
%	'OVERFLOWDETECTION'	'On'		[default] 
%				'Off'		no overflow detection 
%	'OUTPUT'		'single'	single precision data [default: 'double'] 
%	'NUMBER_OF_NAN_IN_BREAK'   N		inserts N NaN's between two concatanated segments
%						default: N=100
%	'SampleRate'		Fs		target sampling rate (supports resampling)
%	'EOG_CORRECTION'	'On'		uses two-channel regression analysis - if possible 
%				'Off'		no correction of EOG artifacts [default]
%	
% The list of supported formats is available here: 
% http://hci.tugraz.at/~schloegl/biosig/TESTED
%
%    SLOAD loads all the data (of the selected channels)
%    at once. In case of large data files, This can be 
%    a problem. Instead, You can use 
%       HDR = sopen(FILENAME,'r');
%       [s,HDR]=sread(HDR,duration_segment1); 
%       [s,HDR]=sread(HDR,duration_segment2); 
%       ....
%       [s,HDR]=sread(HDR,duration_segmentM); 
%       HDR = sclose(HDR); 
%
% see also: SVIEW, SOPEN, SREAD, SCLOSE, SAVE2BKR, TLOAD
%
% In order to increase the speed, install mexSLOAD.mex from biosig4c++
%
% Reference(s):


%	$Id: sload.m,v 1.1 2009-01-30 06:04:42 arno Exp $
%	Copyright (C) 1997-2007,2008 by Alois Schloegl 
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

%
%    BioSig is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BioSig is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BioSig.  If not, see <http://www.gnu.org/licenses/>.

if length(varargin)<2; 
	MODE = ''; 
else	
	MODE = varargin{2};
end;

CHAN = 0; 
STATE.UCAL = 0; 
STATE.OVERFLOWDETECTION = 1; 
STATE.NUMBER_OF_NAN_IN_BREAK = 100; 
STATE.EOG_CORRECTION = 0 ; 
STATE.OUTPUT = 'double';

Fs = NaN; 
k = 1; 
while (k<=length(varargin))
	if isnumeric(varargin{k})
		if (k==1), CHAN=varargin{k}; end;
	elseif ~isempty(strfind(varargin{k},'UCAL'))
		MODE = varargin{k}; 
		if strcmpi(varargin{k+1},'on'); 
			STATE.UCAL = 1; 
			k=k+1;
		elseif strcmpi(varargin{k+1},'off'); 
			STATE.UCAL = 0; 
			k=k+1;
		else
			STATE.UCAL = 1; 
		end; 
	elseif ~isempty(strfind(varargin{k},'OVERFLOWDETECTION:OFF'))
		MODE = varargin{k}; 
		STATE.OVERFLOWDETECTION = 0;
	elseif strcmpi(varargin{k},'OVERFLOWDETECTION')
		if strcmpi(varargin{k+1},'on'); 
			STATE.OVERFLOWDETECTION = 1;
		elseif strcmpi(varargin{k+1},'off'); 
			STATE.OVERFLOWDETECTION = 0;
		elseif isnumeric(varargin{k+1})
			STATE.OVERFLOWDETECTION = varargin{k+1};
		end;			
		k=k+1;
	elseif ~isempty(strfind(varargin{k},'OUTPUT:SINGLE'))
		MODE = varargin{k}; 
		STATE.OUTPUT = 'single';
	elseif strcmpi(varargin{k},'OUTPUT')
		STATE.OUTPUT = varargin{k+1};
	elseif strcmpi(varargin{k},'EOG_CORRECTION')
		if strcmpi(varargin{k+1},'on'); 
			STATE.EOG_CORRECTION = 1;
		elseif strcmpi(varargin{k+1},'off'); 
			STATE.EOG_CORRECTION = 0;
		elseif isnumeric(varargin{k+1})
			STATE.EOG_CORRECTION = varargin{k+1};
		end;			
		k=k+1;
	elseif strcmpi(varargin{k},'Number_of_nan_in_break')
		STATE.NUMBER_OF_NAN_IN_BREAK = varargin{k+1};
		k=k+1;
	elseif strcmpi(varargin{k},'SampleRate')
		Fs = varargin{k+1};
		k=k+1;
	end; 	
	k=k+1;
end;


%%% resolve wildcards %%%
if ischar(FILENAME) 
if any(FILENAME=='*')
        p = fileparts(FILENAME);
        f = dir(FILENAME);
        EOGix = zeros(1,length(f));
        for k = 1:length(f);
                f(k).name = fullfile(p,f(k).name);
                [p,g,e]   = fileparts(f(k).name);
                lg = length(g);
                if (lg>2) & strcmp(upper(g(lg+(-2:0))),'EOG')
                        EOGix(k) = 1;
                end
        end;
        FILENAME={f([find(EOGix),find(~EOGix)]).name};
end;        
end;

if nargout>0,
	signal = [];
end;

if (iscell(FILENAME) && (length(FILENAME)==1)),
	FILENAME = FILENAME{k};
end; 	
	
if ((iscell(FILENAME) || isstruct(FILENAME))),
	signal = [];
	for k = 1:length(FILENAME),
		if iscell(FILENAME(k))
			f = FILENAME{k};
		else  
			f = FILENAME(k);
		end	

                [s,h] = sload(f,varargin{:});
		if k==1,
			H = h;
			signal = s;  
			H.SegLen = [0,size(s,1)];
                        %H.EVENT.POS = [H.EVENT.POS; 0];
                        %H.EVENT.TYP = [H.EVENT.TYP; hex2dec('7ffe')];
                        %if isfield(H.EVENT,'CHN');
                        %        H.EVENT.CHN = [H.EVENT.CHN; 0];
                        %end;
                        %if isfield(H.EVENT,'DUR');
                        %        H.EVENT.DUR = [H.EVENT.DUR; size(s,1)];
                        %end;
                        %if isfield(H.EVENT,'Desc');	% TFM-Excel-Beat-to-Beat
                        %        H.EVENT.Desc = [H.EVENT.Desc; {'New Segment'}; h.EVENT.Desc];
                        %end;
		else
			H.FILE(k) = h.FILE;
                        H.T0(k,1:6) = h.T0;
			if ~isnan(h.SampleRate) && (H.SampleRate ~= h.SampleRate),
				fprintf(2,'Warning SLOAD: sampling rates of multiple files differ %i!=%i.\n',H.SampleRate, h.SampleRate);
			end;
                        if ~isempty(h.EVENT.POS),
                                H.EVENT.POS = [H.EVENT.POS; size(signal,1); h.EVENT.POS+size(signal,1)];
                                H.EVENT.TYP = [H.EVENT.TYP; hex2dec('7ffe'); h.EVENT.TYP];
                                if isfield(H.EVENT,'CHN');
        	                	H.EVENT.CHN = [H.EVENT.CHN; 0; h.EVENT.CHN];
                                end;
                                if isfield(H.EVENT,'DUR');
      					H.EVENT.DUR = [H.EVENT.DUR; size(s,1); h.EVENT.DUR];
                                end;
                                if isfield(H.EVENT,'Desc');	% TFM-Excel-Beat-to-Beat
                                        H.EVENT.Desc = [H.EVENT.Desc; {'New Segment'}; h.EVENT.Desc];
                                end;
                        end;

                        if size(s,2)==size(signal,2), %(H.NS == h.NS) 
				signal = [signal; repmat(NaN,STATE.NUMBER_OF_NAN_IN_BREAK,size(s,2)); s];
				H.SegLen = [H.SegLen,size(signal,1)];
			else
				error('ERROR SLOAD: incompatible channel numbers %i!=%i of multiple files\n',H.NS,h.NS);
			end;
                        if ~isequal(H.Label,h.Label) || ~isequal(H.PhysDimCode,h.PhysDimCode)
				warning('Labels and PhysDim of multiple files differ!\n');
                                for k2 = 1:length(H.InChanSelect),
                                	k1 = H.InChanSelect(k2); 
                                        if ~strcmp(H.Label{k1},h.Label{k1}) || (H.PhysDimCode(k1)~=h.PhysDimCode(k1)),
                                                warning('#%02i:  %s | %s | (%i)%s | (%i)%s\n',k1, H.Label{k1},h.Label{k1}, H.PhysDimCode(k1),H.PhysDim{k1},h.PhysDimCode(k1),h.PhysDim{k1});
                                        end;
                                end;
                        end;
                        if ~isequal(H.Calib,h.Calib) && H.FLAG.UCAL,
				error('SLOAD: concatanating uncalibrated data with different scaling factors does not make sense!');
                        end;

                        if isfield(h,'TRIG'), 
                                if ~isfield(H,'TRIG'),
                                        H.TRIG = [];
                                end;
                                H.TRIG = [H.TRIG(:); h.TRIG(:)+size(signal,1)-size(s,1)];
                        end;
                        
                        if isfield(H,'TriggerOffset'),
                                if H.TriggerOffset ~= h.TriggerOffset,
                                        fprintf(2,'Warning SLOAD: Triggeroffset %f does not fit.%f \n',H.TriggerOffset,h.TriggerOffset);
                                end;
                        end;
                        if isfield(H,'Classlabel'), 
				if isfield(h,'ArtifactSelection'),
                                	if (any(h.ArtifactSelection>1) || (length(h.ArtifactSelection) < length(h.Classlabel)))
                                        	sel = zeros(size(h.Classlabel));
                                                sel(h.ArtifactSelection) = 1; 
                                        else
                                                sel = h.ArtifactSelection(:);
                                        end;
                                else 
                                	h.ArtifactSelection = repmat(logical(0),length(h.Classlabel),1);        
                                end;         

                                if isfield(H,'ArtifactSelection')
                                	H.ArtifactSelection = [H.ArtifactSelection; h.ArtifactSelection(:)];
                                else
                                	H.ArtifactSelection = [repmat(logical(0),length(H.Classlabel),1); h.ArtifactSelection(:)];
                                end;
                                H.Classlabel = [H.Classlabel(:); h.Classlabel(:)];
                        end;
			clear s;
                end;
	end;
        
	% fprintf(1,'SLOAD: data segments are concatenated with NaNs in between.\n');
	return;	
end;
%%% end of multi-file section 


%%%% start of single file section
%%%%%%%%%% --------- EOG CORRECTION -------------- %%%%%%
if STATE.EOG_CORRECTION,
try
        h = get_regress_eog(FILENAME,'REG fft16');
catch,
	fprintf(H.FILE.stderr,'Error: SLOAD (EOG_CORRECTION): %s\n',lasterr); 
	H.FLAG.EOG_CORRECTION = 0; 
end;
end; 


FlagLoaded = 0;
if exist('mexSLOAD','file')==3,
	try
		valid_rerefmx = 1;
		if all(size(CHAN)>1) | any(floor(CHAN)~=CHAN) | any(CHAN<0) | (any(CHAN==0) & (numel(CHAN)>1));
		        ReRefMx = CHAN; 
	        	CHAN = find(any(CHAN,2));
		elseif all(CHAN>0) & all(floor(CHAN)==CHAN), 
			[tmp,ix]= sort(CHAN);
		        ReRefMx = sparse(CHAN,1:length(CHAN),1);
		else    
		        ReRefMx = [];
			valid_rerefmx=0;
		end

		if STATE.EOG_CORRECTION,
			ReRefMx = h.r0*ReRefMx;
			valid_rerefmx=1;
		end
		if STATE.OVERFLOWDETECTION,
			arg1 = 'OVERFLOWDETECTION:ON';
		else
			arg1 = 'OVERFLOWDETECTION:OFF';
		end			
		if STATE.UCAL,
			arg2 = 'UCAL:ON';
		else
			arg2 = 'UCAL:OFF';
		end
		if ~valid_rerefmx,
			[signal,HDR] = mexSLOAD(FILENAME,0,arg1,arg2);
			FlagLoaded = isfield(HDR,'NS');
			HDR.InChanSelect = 1:HDR.NS;
		else
			InChanSelect = find(any(ReRefMx,2));
			[signal,HDR] = mexSLOAD(FILENAME,InChanSelect,arg1,arg2);
			FlagLoaded = isfield(HDR,'NS');
			HDR.InChanSelect = InChanSelect(InChanSelect <= HDR.NS);
			signal = signal*ReRefMx(HDR.InChanSelect,:); %% can be sparse if just a single channel is loaded
			signal = full(signal); 	%% make sure signal is not sparse 
		end; 
		
		HDR.T0 = datevec(HDR.T0);
		HDR.Patient.Birthday = datevec(HDR.Patient.Birthday);
		HDR.Calib = [HDR.Off(:)';diag(HDR.Cal)];
		[HDR.FILE.Path,HDR.FILE.Name,HDR.FILE.Ext] = fileparts(FILENAME);
		HDR.FileName = FILENAME;
		HDR = leadidcodexyz(HDR);
		HDR.EVENT.POS = HDR.EVENT.POS+1; % convert from 0-based to 1-based index
		
		H=HDR;
		H.FLAG.EOG_CORRECTION = STATE.EOG_CORRECTION; 

		if isfield(HDR,'Patient') && isfield(HDR.Patient,'Weight') && isfield(HDR.Patient,'Height')
			%% Body Mass Index 
       			HDR.Patient.BMI = HDR.Patient.Weight * HDR.Patient.Height^-2 * 1e4;

		       	%% Body Surface Area
			% DuBois D, DuBois EF. A formula to estimate the approximate surface area if height and weight be known. Arch Intern Medicine. 1916; 17:863-71.
			% Wang Y, Moss J, Thisted R. Predictors of body surface area. J Clin Anesth. 1992; 4(1):4-10.
		       	HDR.Patient.BSA = 0.007184 * HDR.Patient.Weight^0.425 * HDR.Patient.Height^0.725;
		end; 

		if strcmp(H.TYPE,'GDF');
                        % Classlabels according to 
                        % http://biosig.cvs.sourceforge.net/*checkout*/biosig/biosig/doc/eventcodes.txt
                        if (length(H.EVENT.TYP)>0)
                                ix = (H.EVENT.TYP>hex2dec('0300')) & (H.EVENT.TYP<hex2dec('030d'));
                                ix = ix | ((H.EVENT.TYP>=hex2dec('0320')) & (H.EVENT.TYP<=hex2dec('037f')));
                                ix = ix | (H.EVENT.TYP==hex2dec('030f')); % unknown/undefined cue
                                H.Classlabel = mod(H.EVENT.TYP(ix),256);
                                H.Classlabel(H.Classlabel==15) = NaN; % unknown/undefined cue
                        end;
                        
                        % Trigger information and Artifact Selection 
                        ix = find(H.EVENT.TYP==hex2dec('0300')); 
                        H.TRIG = H.EVENT.POS(ix);
                        ArtifactSelection = repmat(logical(0),length(ix),1);
                        for k = 1:length(ix),
                                ix2 = find(H.EVENT.POS(ix(k))==H.EVENT.POS);
                                if any(H.EVENT.TYP(ix2)==hex2dec('03ff'))
                                        ArtifactSelection(k) = logical(1);                
                                end;
                        end;
                        if any(ArtifactSelection), % define only if necessary
                                H.ArtifactSelection = ArtifactSelection; 
                        end;
                        
			% apply channel selections to EVENT table
			if valid_rerefmx && ~isempty(H.EVENT.POS) && (isfield(H.EVENT,'CHN')),	% only if channels are selected. 
				sel = (H.EVENT.CHN(:)==0);	% memory allocation, select all general events
				for k = find(~sel'),		% select channel specific elements
					sel(k) = any(H.EVENT.CHN(k)==InChanSelect);
				end;
				H.EVENT.POS = H.EVENT.POS(sel);
				H.EVENT.TYP = H.EVENT.TYP(sel);
				H.EVENT.DUR = H.EVENT.DUR(sel);	% if EVENT.CHN available, also EVENT.DUR is defined. 
				H.EVENT.CHN = H.EVENT.CHN(sel);
				% assigning new channel number 
				a = zeros(1,HDR.NS);
				for k = 1:length(InChanSelect),		% select channel specific elements
					a(InChanSelect(k)) = k;		% assigning to new channel number. 
				end;
				ix = H.EVENT.CHN>0;
				H.EVENT.CHN(ix) = a(H.EVENT.CHN(ix));	% assigning new channel number
			end;	
			
		elseif strcmp(H.TYPE,'BKR');
	                H.Classlabel = [];
	                H.TRIG = [];
		        tmp=fullfile(H.FILE.Path,[H.FILE.Name,'.mat']);
        		if ~exist(tmp,'file'),
        		        tmp=fullfile(H.FILE.Path,[H.FILE.Name,'.MAT']);
       			end
		        x = [];
        		if exist(tmp,'file'),
                		x = load('-mat',tmp);
        		end;

                        if isfield(x,'header'),
                                H.MAT  = x.header;
                                if isfield(x.header,'Setup'), 
                                        if isfield(x.header.Setup,'Bits'), 
                                                H.Bits = x.header.Setup.Bits;
                                                [datatyp, limits, datatypes] = gdfdatatype(H.Bits+255);
                                                % THRESHOLD for Overflow detection
                                                if ~isfield(H,'THRESHOLD')
                                                    H.THRESHOLD = repmat(limits, H.NS, 1);
                                                end;         
                                        end;
                                end;
                                if isfield(x.header,'Result') && isfield(x.header.Result,'Classlabel'),
                                        H.Classlabel = x.header.Result.Classlabel;
                                end;
                                if isfield(x.header,'Paradigm')
                                        if isempty(H.Classlabel) && isfield(x.header.Paradigm,'Classlabel')
                                                H.Classlabel = x.header.Paradigm.Classlabel;
                                        end;
                                        H.BCI.Paradigm = x.header.Paradigm;
                                        if isfield(H.BCI.Paradigm,'TriggerOnset');
                                                H.TriggerOffset = H.BCI.Paradigm.TriggerOnset;
                                        elseif isfield(H.BCI.Paradigm,'TriggerTiming');
                                            %    H.BCI.Paradigm.TriggerTiming,
                                                H.TriggerOffset = H.BCI.Paradigm.TriggerTiming;
                                                fprintf(2,'Warning BKROPEN: Paradigm.TriggerOnset is unknown. Paradigm.TriggerTiming= %f ms is used instead\n',H.TriggerOffset);
                                        end;
                                end;

                                if isfield(x.header,'PhysioRec'), % R. Leeb's data 
                                       H.Label = cellstr(x.header.PhysioRec);
                                end;
                                if isfield(x.header,'BKRHeader'), % R. Scherer Data 
                                        if isfield(x.header.BKRHeader,'TO'),
                                                H.T0 = x.header.BKRHeader.TO;
                                        end;
                                        if isfield(x.header.BKRHeader,'Label'),
                                                H.Label = cellstr(x.header.BKRHeader.Label);
                                                ns = H.NS-length(H.Label);
                                                if ns == 1;
                                                        H.Label = strvcat(H.Label,'TRIGGER');
                                                elseif ns > 1;
                                                        H.Label = strvcat(H.Label,char(repmat('n.a.',ns,1)));
                                                end;
                                        end;
                                end;
                                if isfield(x.header,'Model'), % More 
                                        if isfield(x.header.Model,'AnalogInput'), 
                                                for k = 1:length(x.header.Model.AnalogInput),
                                                        H.Filter.HighPass(k) = x.header.Model.AnalogInput{k}{5};
                                                        H.Filter.LowPass(k)  = x.header.Model.AnalogInput{k}{6};
                                                        H.Filter.Notch(k)    = strcmpi(x.header.Model.AnalogInput{k}{7},'on');

                                                        H.MAT.Cal(k) = x.header.Model.AnalogInput{k}{3};
                                                end
                                        end;
                                end;
                                if ~isempty(strmatch('TRIGGER',H.Label))
                                        H.AS.TRIGCHAN = H.NS; %strmatch('TRIGGER',H.Label); 
                                end;
                        end;
                        if 1; ~isfield(H,'Classlabel');
                                tmp=fullfile(H.FILE.Path,[H.FILE.Name,'.par']);
                                if ~exist(tmp,'file'),
                                    tmp=fullfile(H.FILE.Path,[H.FILE.Name,'.PAR']);
                                end
                                if exist(tmp,'file'),
                                        H.Classlabel = load(tmp);
                                end;
                        end;

                        %%% Artifact Selection files 
                        tmp1=fullfile(H.FILE.Path,[H.FILE.Name,'.sel']);
                        if ~exist(tmp1,'file'),
                                tmp1=fullfile(H.FILE.Path,[H.FILE.Name,'.SEL']);
                        end
                        tmp2 = fullfile(H.FILE.Path,[H.FILE.Name,'_artifact.mat']);
                        SW   = (exist(tmp1,'file')>0) + 2*(exist(tmp2,'file')>0);
                        if SW == 0, 
                        elseif SW == 1,
                                if exist('OCTAVE_VERSION','builtin')
                                        H.ArtifactSelection = load('-ascii',tmp1);
                                else
                                        H.ArtifactSelection = load(tmp1);
                                end;
                        elseif SW == 2,
                                if exist('OCTAVE_VERSION','builtin')
                                        tmp = load('-mat',tmp2);
                                else
                                        tmp = load(tmp2);
                                end;
                                H.ArtifactSelection = tmp.artifact(:);
                        elseif SW == 3,
                                fprintf(H.FILE.stderr,'Warning BKROPEN: more than one ArtifactSelection files. File %s is used.\n',tmp1);
                                if exist('OCTAVE_VERSION')>5
                                        H.ArtifactSelection = load('-ascii',tmp1);
                                else
                                        H.ArtifactSelection = load(tmp1);
                                end;
                        end;
                        if isfield(H,'ArtifactSelection'),
                                if any(H.ArtifactSelection>1) | (length(H.ArtifactSelection)<length(H.Classlabel))
                                        sel = zeros(size(H.Classlabel));
                                        sel(H.ArtifactSelection) = 1;
                                        H.ArtifactSelection = sel(:);
                                end;
                                H.ArtifactSelection = H.ArtifactSelection(:);
                        end;

                        if isfield(H.AS,'TRIGCHAN') % & isempty(H.EVENT.POS)
                                if H.AS.TRIGCHAN <= H.NS, %size(H.data,2),
                                        H.THRESHOLD(H.AS.TRIGCHAN,1:2) = [-1-2^15,2^15]; % do not apply overflow detection for Trigger channel 
                                        data = mexSLOAD(H.FileName,H.AS.TRIGCHAN,'UCAL:ON','OVERFLOWDETECTION:OFF');
                                        TRIGon = gettrigger(data);
                                        %TRIGoff = gettrigger(-double(H.data(:,H.AS.TRIGCHAN)));
                                        if isfield(H,'TriggerOffset')
                                                TRIGon  = TRIGon - round(H.TriggerOffset/1000*H.SampleRate);
                                        %        TRIGoff = TRIGoff - round(H.TriggerOffset/1000*H.SampleRate);
                                        end;
                                end;
                                H.TRIG = TRIGon(:);
                                H.EVENT.POS = TRIGon(:); %[TRIGon(:); TRIGoff(:)]; 
                                H.EVENT.TYP = repmat(hex2dec('0300'),numel(TRIGon),1); %repmat(hex2dec('8300'),numel(TRIGoff),1)];
                        end;
                        if length(H.TRIG)~=length(H.Classlabel),
                                % hack to deal with BCI22 data
                                fprintf(2,'Warning BKROPEN: Number of triggers (%i) and number of Classlabels (%i) do not fit\n',length(H.TRIG),length(H.Classlabel));
                                H.TRIG = [];
                                H.Classlabel = [];
                                H.ArtifactSelection = [];
                        end;
			%% end of BKR

		elseif strcmp(H.TYPE,'BrainVision');
        	        try 
        	        	H = bv2biosig_events(H); 
       	        	catch
       	        		%warning('bv2biosig_events not executed');
       	        	end; 

		elseif strcmp(H.TYPE,'EDF');
			if length(HDR.EVENT.TYP)==length(HDR.EVENT.Desc)
	                        [HDR.EVENT.CodeDesc, CodeIndex, HDR.EVENT.TYP] = unique(HDR.EVENT.Desc);
	                end;

			%% end of BKR

		elseif strcmp(H.TYPE,'BrainVision');
			HDR = sopen(fullfile(H.FILE.Path,[H.FILE.Name,'.vmrk']));
			HDR = sclose(HDR);
			H.EVENT = HDR.EVENT;

                end;        

	        H.CHANTYP = repmat(' ',1,H.NS);
		for k=1:H.NS,
			if     ~isempty(strfind(lower(H.Label{k}),'eeg')) 	H.CHANTYP(k) = 'E';
			elseif ~isempty(strfind(lower(H.Label{k}),'meg')) 	H.CHANTYP(k) = 'E';
			elseif ~isempty(strfind(lower(H.Label{k}),'emg')) 	H.CHANTYP(k) = 'M';
			elseif ~isempty(strfind(lower(H.Label{k}),'eog')) 	H.CHANTYP(k) = 'O';
			elseif ~isempty(strfind(lower(H.Label{k}),'ecg')) 	H.CHANTYP(k) = 'C';
			elseif ~isempty(strfind(lower(H.Label{k}),'air')) 	H.CHANTYP(k) = 'R';
			elseif ~isempty(strfind(lower(H.Label{k}),'trig')) 	H.CHANTYP(k) = 'T';
			end; 
		end;

	catch
		%fprintf(1,lasterr);
		fprintf(1, 'SLOAD: mexSLOAD failed - the slower M-function is used.\n');
	end;
else 
	global FLAG_HINT_mexSLOAD;
	if isempty(FLAG_HINT_mexSLOAD) 
		fprintf(1, 'Hint: the performance of SLOAD can be improved with mexSLOAD.mex which is part of biosig4c++.\n');
		FLAG_HINT_mexSLOAD = 1;	 	% turn off hint 
	end;
end;

if ~FlagLoaded,
	H = getfiletype(FILENAME);

	if isempty(H)
		fprintf(2,'Warning SLOAD: no file found\n');
		return;
	else	
		% FILENAME can be fn.name struct, or HDR struct. 
		FILENAME = H.FileName; 
	end;
	H.FLAG.UCAL = STATE.UCAL;
	H.FLAG.OVERFLOWDETECTION = STATE.OVERFLOWDETECTION;
	H.FLAG.EOG_CORRECTION = STATE.EOG_CORRECTION; 
	H.FLAG.OUTPUT = STATE.OUTPUT; 


if strncmp(H.TYPE,'IMAGE:',5)
	[H] = iopen(H);
	if H.FILE.OPEN,
		signal = iread(H);
		H.FILE.OPEN = 0; 
		fclose(H.FILE.FID);
	end;
	return;
end;

%%%%%%%%%%%%%%% --------- Load single file ------------%%%%%%%%%%%

H = sopen(H,'r',CHAN,MODE);
if 0, ~isnan(H.NS),
%------ ignore 'NaC'-channels
	NS = size(H.Calib,2);
	SelMx = speye(NS); 
	ch = 1:NS; %ch(strmatch('NaC',H.Label))=[];
	if length(ch)<NS,
		fprintf(2,'Warning SLOAD: Some NaC channels have been removed %s\n',H.FileName); 
	end; 
	SelMx = SelMx(:,ch); 
	H.Calib = H.Calib*SelMx; 

	% generate correction matrix 
	if H.FLAG.EOG_CORRECTION, 
		H.Calib = H.Calib*h.REGRESS.r0; 
		if ~isequal(H.Label(ch),h.Label)
			fprintf(2, 'Warning SLOAD (EOG correction): Channel labels in %s do not fit with arti* recording!\n',H.FileName); 
		end; 
	elseif STATE.EOG_CORRECTION, 
		fprintf(2, 'Warning: EOG correction not supported for this file (%s)!\n',H.FileName); 
	end; 

	%----- generate HDR.Calib -----
	if all(size(CHAN)>1) || any(floor(CHAN)~=CHAN) || any(CHAN<0) || (any(CHAN==0) && (numel(CHAN)>1));
        	ReRefMx = CHAN; 
        	CHAN = find(any(CHAN,2));
	elseif all(CHAN>0) && all(floor(CHAN)==CHAN), 
		if any(diff(CHAN)<=0),
		%	fprintf(HDR.FILE.FID,'Warning SOPEN: CHAN-argument not sorted - header information like Labels might not correspond to data.\n');
		end;	
        	ReRefMx = sparse(CHAN,1:length(CHAN),1,H.NS,length(CHAN));
	else %if (CHAN==0)    
		ReRefMx = []; 
	end
	if ~isempty(ReRefMx),
		%[size(H.Calib),size(SelMx),size(h.REGRESS.r0),size(ReRefMx)] 
		ReRefMx = [ReRefMx(1:min(size(ReRefMx,1),size(H.Calib,2)),:); zeros(max(0,size(H.Calib,2)-size(ReRefMx,1)),size(ReRefMx,2))]; 
		H.Calib = H.Calib*ReRefMx; 
	end; 	

	H.InChanSelect = find(any(H.Calib(2:end,:),2));
	H.Calib = H.Calib([1;1+H.InChanSelect(:)],:);
end
	

if 0,
        
elseif isfield(H,'data') && all(diag(H.Calib(2:end,1:end))==1)
	% only a single copy of the data
	% important for large data sets, close to the available memory 
	signal = H.data; 
	H = rmfield(H,'data'); 
	if (CHAN>0)
		signal = signal(:,CHAN);
	end; 	 


elseif any(strmatch(H.TYPE,{'native','TFM_EXCEL_Beat_to_Beat','EEProbe-CNT','EEProbe-AVR'})); 
	[signal,H] = sread(H);  
        H = sclose(H);

elseif (H.FILE.OPEN > 0)
	[signal,H] = sread(H);  
        H = sclose(H);

elseif (H.FILE.OPEN > 0)
        signal = repmat(NaN,H.SPR*H.NRec,size(H.Calib,2));
	k1 = 0;
        while ~seof(H),
	        [s,H] = sread(H,100);
        	if 1,
			k2 = size(s,1);
			signal(k1+1:k1+k2,:)=s;
			k1 = k1+k2;
		else
	        	 %fprintf('%i/%i\n',stell(H),H.SPR*H.NRec);
		        signal=[signal;s];
		end;
	end;      
        H = sclose(H);

        
elseif strncmp(H.TYPE,'EVENT',5)
        signal = H.EVENT;
        return; 
        

elseif strncmp(H.TYPE,'ELPOS',5)
        signal = H.ELEC.XYZ;
	return; 

elseif strcmp(H.TYPE,'DAQ')
	fprintf(1,'Loading a matlab DAQ data file - this can take a while.\n');
	tic;
        [signal, tmp, H.DAQ.T0, H.DAQ.events, DAQ.info] = daqread(H.FileName);
        fprintf(1,'Loading DAQ file finished after %.0f s.\n',toc);
        [H.SPR,H.NS] = size(signal);
        
        H.SampleRate = DAQ.info.ObjInfo.SampleRate;
        sz     = size(signal);
        if length(sz)==2, sz=[1,sz]; end;
        H.NRec = sz(1);
        H.Dur  = sz(2)/H.SampleRate;
        H.NS   = sz(3);
        H.FLAG.TRIGGERED = H.NRec>1;
        H.FLAG.UCAL = 1;
        
        H.PhysDim = {DAQ.info.ObjInfo.Channel.Units};
        H.DAQ   = DAQ.info.ObjInfo.Channel;
        
        H.Cal   = diff(cat(1,DAQ.info.ObjInfo.Channel.InputRange),[],2).*(2.^(-DAQ.info.HwInfo.Bits));
        H.Off   = cat(1,DAQ.info.ObjInfo.Channel.NativeOffset); 
        H.Calib = sparse([H.Off';eye(H.NS)]*diag(H.Cal));
        
        if CHAN<1,
                CHAN = 1:H.NS; 
        end;
        if ~H.FLAG.UCAL,
                Calib = H.Calib;	% Octave can not index sparse matrices within a struct
                signal = [ones(size(signal,1),1),signal]*Calib(:,CHAN);
        end;
        
        
elseif 0, strcmp(H.TYPE,'BIFF'),
	try, 
                [H.TFM.S,H.TFM.E] = xlsread(H.FileName,'Beat-To-Beat');
	        if size(H.TFM.S,1)+1==size(H.TFM.E,1),
                        H.TFM.S = [repmat(NaN,1,size(H.TFM.S,2));H.TFM.S];
                end;
                H.TYPE = 'TFM_EXCEL_Beat_to_Beat'; 
	catch
		try
	                [H.TFM.S,H.TFM.E] = xlsread(H.FileName,'Beat-to-Beat');
        	        H.TYPE = 'TFM_EXCEL_Beat_to_Beat'; 
        	catch 
			0; 
        	end;
	end; 	

	if strcmp(H.TYPE, 'TFM_EXCEL_Beat_to_Beat');
                if ~isempty(strfind(H.TFM.E{3,1},'---'))
                        H.TFM.S(3,:) = [];    
                        H.TFM.E(3,:) = [];    
                end;
                
                H.Label   = strvcat(H.TFM.E{4,:});
                H.PhysDim = strvcat(H.TFM.E{5,:});
           
                H.TFM.S = H.TFM.S(6:end,:);
                H.TFM.E = H.TFM.E(6:end,:);
                
                ix = find(isnan(H.TFM.S(:,2)) & ~isnan(H.TFM.S(:,1)));
                H.EVENT.Desc = H.TFM.E(ix,2);
                H.EVENT.POS  = ix;
                
                if any(CHAN),
			H.TFM.S = H.TFM.S(:,CHAN);
			H.TFM.E = H.TFM.E(:,CHAN);
		end;
		[H.SPR,H.NS] = size(H.TFM.S);
		H.NRec = 1; 
		H.THRESHOLD  = repmat([0,NaN],H.NS,1);

		signal  = H.TFM.S;
		signal(signal==0) = NaN;
        end;


elseif strcmp(H.TYPE,'MatrixMarket'),
        H.FILE.FID = fopen(H.FileName,'rt','ieee-le');

    	line = fgetl(H.FILE.FID);
	
	H.FLAG.Coordinate = ~isempty(strfind(line,'coordinate'));
	H.FLAG.Array 	  = ~isempty(strfind(line,'array'));

	H.FLAG.Complex = ~isempty(strfind(line,'complex'));
	H.FLAG.Real = ~isempty(strfind(line,'real'));
	H.FLAG.Integer = ~isempty(strfind(line,'integer'));
	H.FLAG.Pattern = ~isempty(strfind(line,'pattern'));
	
	H.FLAG.General = ~isempty(strfind(line,'general'));
	H.FLAG.Symmetric = ~isempty(strfind(line,' symmetric'));
	H.FLAG.SkewSymmetric = ~isempty(strfind(line,'skew-symmetric'));
	H.FLAG.Hermitian = ~isempty(strfind(lower(line),'hermitian'));

	while strncmp(line,'%',1)
        	line = fgetl(H.FILE.FID);
	end;

	[tmp,status] = str2double(line);
	if any(status)
		fprintf(H.FILE.stderr,'SLOAD (MM): invalid size %s\n',line);
	else
		H.MATRIX.Size = tmp;
	end;	

	if length(H.MATRIX.Size)==3,
		H.Length = tmp(3);
		signal = sparse([],[],[],tmp(1),tmp(2),tmp(3));
		for k = 1:H.Length,
	        	line = fgetl(H.FILE.FID);
			[tmp,status] = str2double(line);
			if any(status)
				fprintf(H.FILE.stderr,'SLOAD (MM): invalid size %s\n',line);
			elseif length(tmp)==4,	
		    		val = tmp(3) + i*tmp(4);
			elseif length(tmp)==3,	
				val = tmp(3);
			elseif length(tmp)==2,	
				val = 1;
			else
				fprintf(H.FILE.stderr,'SLOAD (MM): invalid size %s\n',line);
			end;

			if H.FLAG.General,
				signal(tmp(1),tmp(2)) = val;
			elseif H.FLAG.Symmetric,
				signal(tmp(1),tmp(2)) = val;
				signal(tmp(2),tmp(1)) = val;
			elseif H.FLAG.SkewSymmetric,
				signal(tmp(1),tmp(2)) = val;
				signal(tmp(2),tmp(1)) =-val;
			elseif H.FLAG.Hermitian,
				signal(tmp(1),tmp(2)) = val;
				signal(tmp(2),tmp(1)) = conj(val);
			else	
				fprintf(H.FILE.stderr,'SLOAD (MM): invalid size %s\n',line);
			end;	
		end;
					
	elseif length(H.MATRIX.Size)==2
		H.Length = prod(tmp);
		signal = zeros(H.MATRIX.Size);
		if H.FLAG.General==1,
			[IX,IY]=find(ones(H.MATRIX.Size));
		else
			[IX,IY]=find(cumsum(eye(H.MATRIX.Size)));
		end;
				
		for k = 1:H.Length,
	        	line = fgetl(H.FILE.FID);
			[tmp,status] = str2double(line);
			if any(status)
				error('SLOAD (MM)');
			elseif length(tmp)==2,	
				val=tmp(1) + i*tmp(2);
			elseif length(tmp)==1,	
				val=tmp(1);
			else
				fprintf(H.FILE.stderr,'SLOAD (MM): invalid size %s\n',line);
			end;

			signal(IX(k),IY(k)) = val;
			if H.FLAG.Symmetric,
				signal(IY(k),IX(k)) = val;
			elseif H.FLAG.SkewSymmetric,
				signal(IY(k),IX(k)) =-val;
			elseif H.FLAG.Hermitian,
				signal(IY(k),IX(k)) = conj(val);
			else	
				fprintf(H.FILE.stderr,'SLOAD (MM): invalid size %s\n',line);
			end;	
		end;
        end;
        fclose(H.FILE.FID);

        
elseif strcmp(H.TYPE,'OFF'),
	        H.FILE.FID = fopen(H.FileName,'rt','ieee-le');
		
                line1 = fgetl(H.FILE.FID);
                line2 = fgetl(H.FILE.FID);
		while ~feof(H.FILE.FID) & (line2(1)=='#'),
	                line2 = fgetl(H.FILE.FID);
		end;
		[tmp,status] = str2double(line2);
		if status | (size(tmp,2)~=3), 
			fclose(H.FILE.FID);
			error('SOPEN (OFF)');
		else
			H.VertexCount = tmp(1);
			H.FaceCount = tmp(2);
			H.EdgeCount = tmp(3);
		end	
		
		H.Vertex = repmat(NaN,H.VertexCount,3);
		for k = 1:H.VertexCount,
			line = '';
			while isempty(line) | strncmp(line,'#',1)
				line = fgetl(H.FILE.FID);
			end;
			len = min(length(line),min(find(line=='#')));
			tmp = str2double(line(1:len));
			H.Vertex(k,:) = tmp(1:H.ND);
		end;	
		
%		H.Face = repmat(NaN,H.FaceCount,3);
		for k = 1:H.FaceCount,
			line = '';
			while isempty(line) | strncmp(line,'#',1)
				line = fgetl(H.FILE.FID);
			end;
			len = min(length(line),min(find(line=='#')));
			tmp = str2double(line(1:len));
			H.Ngon(k) = tmp(1);
			H.Face{k} = tmp(2:tmp(1)+1) + 1;
		end;	
		if all(H.Ngon(1)==H.Ngon),
			H.Face = cat(1,H.Face{:});
		end;	
                fclose(H.FILE.FID);
        
        
elseif strcmp(H.TYPE,'POLY'),
        H.FILE.FID = fopen(H.FileName,'rt','ieee-le');

	K = 0;
	while ~feof(H.FILE.FID)
        	line = fgetl(H.FILE.FID);
		if isempty(line),
		elseif line(1)=='#',
		else
			K = K + 1;
			
		end;
        end;
        fclose(H.FILE.FID);
        
        
elseif strcmp(H.TYPE,'SMF'),
	        H.FILE.FID = fopen(H.FileName,'rt','ieee-le');
		
		VertexCount = 0;
		FaceCount = 0;
		PalLen = 0; 
		K = 1;
		while ~feof(H.FILE.FID)
	                line = fgetl(H.FILE.FID);
			if isempty(line)
			elseif line(1)=='#';

			elseif line(1)=='v';
				[tmp,status] = str2double(line(3:end));
				if ~any(status)
					VertexCount = VertexCount + 1 ;
					H.Vertex(VertexCount,:) = tmp;
				else
					fprintf(H.FILE.stderr,'Warning SLOAD: could not read line %i in file %s\n',K,H.FileName); 	
				end;	

			elseif line(1)=='f';
				[tmp,status] = str2double(line(3:end));
				if ~any(status)
					FaceCount  = FaceCount + 1; 
					H.Ngon(FaceCount) = length(tmp);
					H.Face{FaceCount} = tmp;
				else
					fprintf(H.FILE.stderr,'Warning SLOAD: could not read line %i in file %s\n',K,H.FileName); 	
				end;	

			elseif line(1)=='n';
				[tmp,status] = str2double(line(3:end));
				if ~any(status)
					H.NormalVector = tmp;
				else
					fprintf(H.FILE.stderr,'Warning SLOAD: could not read line %i in file %s\n',K,H.FileName); 	
				end;	

			elseif line(1)=='c';
				[tmp,status] = str2double(line(3:end));
				if ~any(status)
					PalLen = PalLen +1; 
					H.Palette(PalLen,:)= tmp;
				else
					fprintf(H.FILE.stderr,'Warning SLOAD: could not read line %i in file %s\n',K,H.FileName); 	
				end;	
			else

			end;
			K = K+1;
		end;
		fclose(H.FILE.FID);
		if all(H.Ngon(1)==H.Ngon),
			H.Face = cat(1,H.Face{:});
		end;	
	
        
elseif strcmp(H.TYPE,'TVF 1.1A'),
        H.FILE.FID = fopen(H.FileName,'rt');

	tmp = fgetl(H.FILE.FID);
	H.TVF.Name = fgetl(H.FILE.FID);
	H.TVF.Desc = fgetl(H.FILE.FID);
	[tmp,status] = str2double(fgetl(H.FILE.FID));
	H.TVF.NTR = tmp(1);
	H.TVF.NV  = tmp(2);
	H.FLAG.CULL = ~~tmp(3);
	[tmp,status] = str2double(fgetl(H.FILE.FID));
	H.TVF.NTRC  = tmp(1);
	H.TVF.NTRM  = tmp(2);
	[tmp,status] = str2double(fgetl(H.FILE.FID));
	H.TVF.NVC   = tmp(1);
	H.TVF.NVN   = tmp(2);
	[tmp,status] = str2double(fgetl(H.FILE.FID));
	H.TVF.GlobalColor = tmp;
	[tmp,status] = str2double(fgetl(H.FILE.FID));
	H.TVF.GlobalMtrlProps = tmp;
	
	H.TVF.Triangles = repmat(NaN,[H.TVF.NTR,3]);
	for k = 1:H.TVF.NTR,
		[tmp, status] = str2double(fgetl(H.FILE.FID));
		H.TVF.Triangles(k,:) = tmp;
	end;	
	H.TVF.TrColorSets = reshape(NaN,[H.TVF.NTRC,H.TVF.NTR]);
	for k = 1:H.TVF.NTRC,
		[tmp, status] = str2double(fgetl(H.FILE.FID));
		H.TVF.TrColorSets(k,:) = tmp;
	end;
	H.TVF.TrMtrlSets = repmat(NaN, [H.TVF.NTRM,H.TVF.NTR]);
	for k = 1:H.TVF.NTRM,
		[tmp, status] = str2double(fgetl(H.FILE.FID));
		H.TVF.TrMtrlSets(k,:) = tmp;
	end;

	H.TVF.Vertices   = repmat(NaN, [H.TVF.NV,3]);
	for k = 1:H.TVF.NV,
		[tmp, status] = str2double(fgetl(H.FILE.FID));
		H.TVF.Vertices(k,:) = tmp;
	end;
	H.TVF.VColorSets = repmat(NaN, [H.TVF.NVC,H.TVF.NV]);
	for k = 1:H.TVF.NVC,
		[tmp, status] = str2double(fgetl(H.FILE.FID));
		H.TVF.VColorSets(k,:) = tmp;
	end;
	H.TVF.VNrmlSets  = repmat(NaN, [H.TVF.NVN,3]);
	for k = 1:H.TVF.NVN,
		[tmp, status] = str2double(fgetl(H.FILE.FID));
		H.TVF.VNrmlSets(k,:) = tmp;
	end;

	fclose(H.FILE.FID);	

        
elseif strcmp(H.TYPE,'TVF 1.1B'),
        H.FILE.FID = fopen(H.FileName,'rb',H.Endianity);

	tmp = fread(H.FILE.FID,12,'uchar');
	H.TVF.Name = char(fread(H.FILE.FID,32,'uchar'));
	H.TVF.Desc = char(fread(H.FILE.FID,80,'uchar'));
	tmp = fread(H.FILE.FID,7,'uint32');
	H.TVF.NTR = tmp(1);
	H.TVF.NV  = tmp(2);
	H.FLAG.CULL = ~~tmp(3);
	H.TVF.NTRC  = tmp(4);
	H.TVF.NTRM  = tmp(5);
	H.TVF.NVC   = tmp(6);
	H.TVF.NVN   = tmp(7);
	tmp = fread(H.FILE.FID,[4,2],'float32');
	H.TVF.GlobalColor = tmp(:,1)';
	H.TVF.GlobalMtrlProps = tmp(:,2)';

	H.TVF.Triangles   = fread(H.FILE.FID, [3,H.TVF.NTR],'uint32')';
	H.TVF.TrColorSets = fread(H.FILE.FID, [H.TVF.NTR,H.TVF.NTRC],'uint32')';
	H.TVF.TrMtrlSets  = fread(H.FILE.FID, [H.TVF.NTR,H.TVF.NTRM],'uint32')';

	H.TVF.Vertices   = fread(H.FILE.FID, [3,H.TVF.NV],'float32')';
	H.TVF.VColorSets = fread(H.FILE.FID, [H.TVF.NV,H.TVF.NVC],'uint32')';
	H.TVF.VNrmlSets  = fread(H.FILE.FID, [3,H.TVF.NVN],'float32')';

	fclose(H.FILE.FID);	

        
elseif strcmp(H.TYPE,'VTK'),
                H.FILE.FID = fopen(H.FileName,'rt','ieee-le');
                
                H.VTK.version = fgetl(H.FILE.FID);
                H.VTK.Title   = fgetl(H.FILE.FID);
                H.VTK.type    = fgetl(H.FILE.FID);
		
		while ~feof(H.FILE.FID),
			tline = fgetl(H.FILE.FID);
			
			if 0, 
			elseif strncmp(tline,'CELLS',5);
			elseif strncmp(tline,'CELL_DATA',9);
			elseif strncmp(tline,'DATASET',7);
		                H.VTK.DATASET = tline(8:end);
			elseif strncmp(tline,'LINES',6);
			elseif strncmp(tline,'POINTS',6);
			elseif strncmp(tline,'POINT_DATA',10);
			elseif strncmp(tline,'POLYGON',7);
			elseif strncmp(tline,'SCALARS',7);
				[t1,r]=strtok(tline);
				[dataName,r]=strtok(r);
				[dataType,r]=strtok(r);
				[numComp ,r]=strtok(r);
				if isempty(numComp), numComp=1;
				else numComp = str2double(numComp); end;
				tline = fgetl(fid);
				if strcmp(tline,'LOOKUP_TABLE');
	    				[t1,r]=strtok(tline);
	    				[tableName,r]=strtok(tline);
				else	
					tline = fgetl(fid);
				end;
			end;

		
		end;
                fclose(H.FILE.FID);

                fprintf(H.FILE.stderr,'Warning SOPEN: VTK-format not supported, yet.\n');
	
        
elseif strcmp(H.TYPE,'unknown')
        TYPE = upper(H.FILE.Ext);
        if 0, 
        elseif strcmp(TYPE,'RAW')
                loadraw;
        elseif strcmp(TYPE,'RDT')
                [signal] = loadrdt(FILENAME,CHAN);
                fs = 128;
        elseif strcmp(TYPE,'XLS')
                loadxls;
        elseif strcmp(TYPE,'DA_')
                fprintf('Warning SLOAD: Format DA# in testing state and is not supported\n');
                loadda_;
        elseif strcmp(TYPE,'RG64')
                [signal,H.SampleRate,H.Label,H.PhysDim,H.NS]=loadrg64(FILENAME,CHAN);
                %loadrg64;
        else
                fprintf('Error SLOAD: Unknown Data Format\n');
                signal = [];
        end;
end;

end; 

%%%%%%%%%% Post-Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(H.TYPE,'CNT');    
        f = fullfile(H.FILE.Path, [H.FILE.Name,'.txt']); 
        if exist(f,'file'),
                fid = fopen(f,'r');
		tmp = fread(fid,inf,'uint8');
		fclose(fid);
		[tmp,v] = str2double(char(tmp'));
		if ~any(v), 
            		H.Classlabel=tmp(:);                        
	        end;
        end
        f = fullfile(H.FILE.Path, [H.FILE.Name,'.par']); 
        if exist(f,'file'),
                fid = fopen(f,'r');
		tmp = fread(fid,inf,'uint8');
		fclose(fid);
		[tmp,v] = str2double(char(tmp'));
		if ~any(v), 
            		H.Classlabel=tmp(:);                        
	        end;
        end
        f = fullfile(H.FILE.Path, [H.FILE.Name,'.mat']);
        if exist(f,'file'),
                try,
                        tmp = load(f);
                catch
                        tmp = []; 
                end
                if isfield(tmp,'classlabel') & ~isfield(H,'Classlabel')
                        H.Classlabel=tmp.classlabel(:);                        
                elseif isfield(tmp,'classlabel') & isfield(tmp,'header') & isfield(tmp.header,'iniFile') & strcmp(tmp.header.iniFile,'oom.ini'), %%% for OOM study only. 
                        H.Classlabel=tmp.classlabel(:);                        
                end;
        end;
        f = fullfile(H.FILE.Path, [H.FILE.Name,'c.mat']);
        if exist(f,'file'),
                try,
                        tmp = load(f);
                catch
                        tmp = []; 
                end
                if isfield(tmp,'classlabel') & ~isfield(H,'Classlabel')
                        H.Classlabel=tmp.classlabel(:);                        
                end;
        end;
        f = fullfile(H.FILE.Path, [H.FILE.Name,'_classlabel.mat']);
        if exist(f,'file'),
                try,
                        tmp = load(f);
                catch
                        tmp = []; 
                end
                if isfield(tmp,'Classlabel') & (size(tmp.Classlabel,2)==4)
                        [x,H.Classlabel] = max(tmp.Classlabel,[],2);                        
                end;
                if isfield(tmp,'classlabel') & (size(tmp.classlabel,2)==4)
                        [x,H.Classlabel] = max(tmp.classlabel,[],2);                        
                end;
        end;
        
        f=fullfile(H.FILE.Path,[H.FILE.Name,'.sel']);
        if ~exist(f,'file'),
                f=fullfile(H.FILE.Path,[H.FILE.Name,'.SEL']);
        end
        if exist(f,'file'),
                fid = fopen(f,'r');
		tmp = fread(fid,inf,'uint8');
		fclose(fid);
		[tmp,v] = str2double(char(tmp'));
		if any(isnan(tmp)) |any(tmp~=ceil(tmp)) | any(tmp<0) | (any(tmp==0) & any(tmp>1))
                        fprintf(2,'Warning SLOAD(CNT): corrupted SEL-file %s\n',f);
                else
                        if ~isfield(H,'Classlabel'), 
                                H.Classlabel = H.EVENT.TYP;
                        end;        
                        n = length(H.Classlabel);
                        
                        H.ArtifactSelection = zeros(n,1);
                        if all((tmp==0) | (tmp==1)) & (length(tmp)>1) & (sum(diff(sort(tmp))~=0) ~= length(tmp)-1)
                                H.ArtifactSelection = logical(tmp);         
                        else
                                H.ArtifactSelection(tmp) = 1;         
                        end;
                end;
        end;
end;

if ~isempty(strfind(upper(MODE),'TSD'));
        f = fullfile(H.FILE.Path, [H.FILE.Name,'.tsd']);
        if ~exist(f,'file'),
                        fprintf(2,'Warning SLOAD-TSD: file %s.tsd found\n',H.FILE(1).Name,H.FILE(1).Name);
        else
                fid = fopen(f,'rb');
                tsd = fread(fid,inf,'float');
                fclose(fid);
                nc = size(signal,1)\size(tsd,1);
                if (nc == round(nc)),
                        signal = [signal, reshape(tsd,nc,size(tsd,1)/nc)'];
                else
                        fprintf(2,'Warning SLOAD: size of %s.tsd does not fit to size of %s.bkr\n',H.FILE(1).Name,H.FILE(1).Name);
                end;
        end;
end;


if strcmp(H.TYPE,'GDF')
        if isfield(H.EVENT,'TYP')
                if isempty(H.EVENT.TYP)
                        %%%%% if possible, load Reinhold's configuration files
                        f = fullfile(H.FILE.Path, [H.FILE.Name,'.mat']);
                        if exist(f,'file'),
                                try
                                        x = load(f,'header');
                                catch;
                                        x=[];
                                end;
                                if isfield(x,'header'),
                                        if isfield(x.header,'Paradigm')
                                                fprintf(2,'Warning SLOAD (GDF): Obsolete feature is used (header information loaded from MAT-file %s).\n',H.FileName);
                                                fprintf(2,'This feature will be removed in future. Contact <a.schloegl@ieee.org> if you need this.\n ');
                                                H.BCI.Paradigm = x.header.Paradigm;
                                                if isfield(H.BCI.Paradigm,'TriggerTiming');
                                                        H.TriggerOffset = H.BCI.Paradigm.TriggerTiming;
                                                elseif isfield(H.BCI.Paradigm,'TriggerOnset');
                                                        H.TriggerOffset = H.BCI.Paradigm.TriggerOnset;
                                                end;

                                                if isfield(H,'Classlabel') & isempty(H.Classlabel),
                                                        H.Classlabel = x.header.Paradigm.Classlabel;
                                                end;
                                        end;
                                end;
                        end;
                end;
        end;

        fid=fopen(fullfile(H.FILE.Path,[H.FILE.Name,'.sel']),'r');
        if fid<0,
	        fid=fopen(fullfile(H.FILE.Path,[H.FILE.Name,'.SEL']),'r');
        end
        if fid>0,
                [tmp,c] = fread(fid,[1,inf],'uint8');
                fclose(fid);
                [tmp,v,sa] = str2double(tmp);
                if isempty(sa) || isempty(sa{1})
                        H.ArtifactSelection = repmat(0,length(H.TRIG),1);
                elseif any(isnan(tmp)) |any(tmp~=ceil(tmp)) | any(tmp<0) | (any(tmp==0) & any(tmp>1))
                        fprintf(2,'Warning SLOAD(GDF): corrupted SEL-file %s\n',fullfile(H.FILE.Path,[H.FILE.Name,'.sel']));
		else                        
                        H.ArtifactSelection = zeros(length(H.TRIG),1);
                        if all((tmp==0) | (tmp==1)) & (length(tmp)>1) & (sum(diff(sort(tmp))~=0) ~= length(tmp)-1)
                                H.ArtifactSelection = logical(tmp);
                        else
                                H.ArtifactSelection(tmp) = 1;
                        end;
                end;
        end;
end;



% resampling 
if ~isnan(Fs) & (H.SampleRate~=Fs);
        tmp = ~mod(H.SampleRate,Fs) | ~mod(Fs,H.SampleRate);
        tmp2= ~mod(H.SampleRate,Fs*2.56);
        if tmp,
                signal = rs(signal,H.SampleRate,Fs);
                H.EVENT.POS = H.EVENT.POS/H.SampleRate*Fs;
                if isfield(H.EVENT,'DUR');
                        H.EVENT.DUR = H.EVENT.DUR/H.SampleRate*Fs;
                end;
                if isfield(H,'TRIG');
                        H.TRIG = H.TRIG/H.SampleRate*Fs;
                end;
                H.SampleRate = Fs;
        elseif tmp2,
                x = load('resample_matrix.mat');
                signal = rs(signal,x.T256100);
                if H.SampleRate*100~=Fs*256,
                        signal = rs(signal,H.SampleRate/(Fs*2.56),1);
                end;
                H.EVENT.POS = H.EVENT.POS/H.SampleRate*Fs;
                if isfield(H.EVENT,'DUR');
                        H.EVENT.DUR = H.EVENT.DUR/H.SampleRate*Fs;
                end;
                if isfield(H,'TRIG');
                        H.TRIG = H.TRIG/H.SampleRate*Fs;
                end;
                H.SampleRate = Fs;
        else 
                fprintf(2,'Warning SLOAD: resampling %f Hz to %f Hz not implemented.\n',H.SampleRate,Fs);
        end;                
end;

return; 
ratiomissing = mean(isnan(signal)); 
if any(ratiomissing>.1)
	fprintf(2,'Warning SLOAD: ratio of missing samples exceeds 10%% in file %s.\n',H.FileName);
	ix = find(ratiomissing); 
	fprintf(1,'#%3i:  %4.1f%%\n',[ix;ratiomissing(ix)*100])
end;	
