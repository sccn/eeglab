% sdfopen() - Opens EDF/GDF/SDF files for reading and writing. The EDF format is specified 
%             in [1], GDF in [2]. SDF is the SIESTA convention (see [3], pp.8-9); SDF uses 
%             the EDF format.
% Usage:
%       >> EDF     = sdfopen(Filename,'r' [,CHAN]);  % prepare to read CHAN channels
%       >> EDF     = sdfopen(Filename,'r' [,CHAN [,MODE [,TSR]]]); % all arguments
% % then
%       >> [S,EDF] = sdfread(EDF, SecRead, SecStart); % reads the opened EDF data file.
%       >> [S,EDF] = sdfread(EDF, Inf); % reads the whole opened EDF data file.
%
% Inputs: 
%      CHAN - [int vector] Specifies the channel(s) to read. Else, a re-referencing 
%                  matrix. [Default|0 -> read all channels].
%      MODE  -
%
%   'UCAL' [Default mode] Indicates that no calibration (re-scaling) to physical dim. 
%                  is performed. Outputs are 16-bit integers. 
% 
%   'SIESTA' - Indicates that channels #1-#6 are re-referenced to (M1+M2)/2          
%                 (Note: 'UCAL' overrides 'SIESTA')
%   'AFIR'   - Indicates that Adaptive FIR filtering is used for ECG removal.
%		       Implements Adaptive FIR filtering for ECG removal in EDF/GDF-tb.
% 		       based on the algorithm of Mikko Koivuluoma <k7320@cs.tut.fi>
%                  A delay of EDF.AFIR.delay number of samples has to be considered. 
%   'SIESTA+AFIR' - Does both.
%   'RAW'       - One column represents one EDF-block
%   'Notch50Hz' - Implements a simple FIR-notch filter at 50 Hz
%   'RECG'      - Implements ECG minimization with regression analysis
%   'TECG'      - Implements ECG minimization with template removal (test status)
%   'HPF0.0Hz'  - Implements a high-pass filter (with zero at z=+1, i.e. a differentiator).
%                 In this case, a notch-filter and/or sub-sampling is recommended. 
%   'TAUx.yS'   - Compensates time-constant of x.y seconds
%   'EOG[hvr]'  - Produces HEOG, VEOG and/or REOG output (CHAN not considered)
%   'OVERFLOW'  - Performs overflow detection
%   'Units_Blocks' - Requests the EDF-field arguments to SDFREAD in blocks 
%                    [default is seconds]
%      TSR - [optional] The target (re)sampling rate. Currently, only downsampling 
%            from 256 Hz or 200 Hz to 100 Hz is supported.  The details are described 
%            in the appendix of [4].
%
% Outputs: 
%         EDF - data structure read from the input file header.
%           EDF.ErrNo   ~= 0  Indicates that an error occurred 
%              1: First 8 bytes are not '0       ', violating the EDF spec.
%              2: Invalid date (unable to guess correct date)
%              4: Incorrect date information (later than current date) 
%             16: Incorrect filesize: Header information does not match actual size
%           EDF.FILE.FID = -1 indicates that file has not been opened
%                           (for compatibility with former versions).
%           EDF.WarnNo  ~=0 Indicates EDF structure problems
%              1: ascii(0) in 1st header
%              2: ascii(0) in 2nd header
%              4: invalid SPR
%              4: invalid samples_per_record-values
%              8: date not in EDF-format (tries to guess correct date, see also E2)
%             16: invalid number_of-channels-value
%             32: invalid value of the EDF header length
%             64: invalid value of block_duration
%            128: Polarity of #7 probably inverted  
%
% Example: 
%      To open an EDF/SDF file for writing:
%      >> [EDF] = sdfopen(EDF,'w') % or equivalently
%      >> [EDF] = sdfopen(EDF.FileName,'w',EDF.Dur,EDF.SampleRate);
%
% Note: Fields EDF.FileName, EDF.NS, EDF.Dur and EDF.EDF.SampleRate must be defined.
% 
% Author: (C) 1997-2002 by Alois Schloegl, 15 Jun 2002 #0.85, (Header reworked for 
%         EEGLAB format, Arnaud Delorme and Scott Makeig, 27 Dec 2002)
%
% See also: fopen, SDFREAD, SDFWRITE, SDFCLOSE, SDFSEEK, SDFREWIND, SDFTELL, SDFEOF

% References: 
% [1] Bob Kemp, Alpo Värri, Agostinho C. Rosa, Kim D. Nielsen and John Gade.
%     A simple format for exchange of digitized polygraphic recordings.
%     Electroencephalography and Clinical Neurophysiology, 82 (1992) 391-393.
% See also: http://www.medfac.leidenuniv.nl/neurology/knf/kemp/edf/edf_spec.htm
%
% [2] Alois Schlögl, Oliver Filz, Herbert Ramoser, Gert Pfurtscheller.
%     GDF - A GENERAL DATAFORMAT FOR BIOSIGNALS
%     Technical Report, Department for Medical Informatics, 
%     Universtity of Technology, Graz (1999)
% See also: http://www-dpmi.tu-graz.ac.at/~schloegl/matlab/eeg/gdf4/tr_gdf.ps
%
% [3] The SIESTA recording protocol. 
% See http://www.ai.univie.ac.at/siesta/protocol.html
% and http://www.ai.univie.ac.at/siesta/protocol.rtf 
%
% [4] Alois Schlögl
%     The electroencephalogram and the adaptive autoregressive model: theory and applications. 
%     (ISBN 3-8265-7640-3) Shaker Verlag, Aachen, Germany.
% See also: "http://www.shaker.de/Online-Gesamtkatalog/Details.idc?ISBN=3-8265-7640-3"

% Testing state
%
% (1) reads header information and writes the header; can be used to check SDFOPEN or for correcting headerinformation
% EDF=sdfopen(EDF,'r+'); EDF=sdfclose(EDF); 
% 
% (2a) Minimal requirements for opening an EDF-File
%       EDF.FileName='file.edf'; % define filename
%       EDF.NS = 5; % fix number of channels
%       EDF=sdfopen(EDF,'w');
%           write file
%           define header somewhen before 
%       EDF=sdfclose(EDF); % writes the corrected header information
% 
% (2b) Minimal requirements for opening an EDF-File
%       EDF=sdfopen('file.edf','w',N); % N number of channels
%            .. do anything, e.g define header
%       EDF=sdfopen(EDF,'w+'); % writes Header information

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

function [EDF,H1,h2]=sdfopen(arg1,arg2,arg3,arg4,arg5,arg6)

    INFO='(C) 1997-2002 by Alois Schloegl, 15 Jun 2002 #0.85';
%	a.schloegl@ieee.org
%	Version 0.85

if nargin<2, arg2='r';end;

if 1,exist('OCTAVE_VERSION');
        EDF.AS.Method='sdfopen';
        mfilename=EDF.AS.Method;
else 
        EDF.AS.Method=mfilename;
end;

EDF.AS.Method=[EDF.AS.Method '-' arg2];
EDF.AS.Date=fix(clock);
EDF.AS.Info=INFO;
EDF.AS.Ver = 0.82;

if isstruct(arg1) 
        EDF=arg1; 
        FILENAME=EDF.FileName;
else
        FILENAME=arg1;
end;

H1idx=[8 80 80 8 8 8 44 8 8 4];
H2idx=[16 80 8 8 8 8 8 80 8 32];

%%%%% Define Valid Data types %%%%%%
%GDFTYPES=[0 1 2 3 4 5 6 7 16 17 255+(1:64) 511+(1:64)];
GDFTYPES=[0 1 2 3 4 5 6 7 16 17 255+[1 12 22 24] 511+[1 12 22 24]];

%%%%% Define Size for each data type %%%%%
GDFTYP_BYTE=zeros(1,512+64);
GDFTYP_BYTE(256+(1:64))=(1:64)/8;
GDFTYP_BYTE(512+(1:64))=(1:64)/8;
GDFTYP_BYTE(1:18)=[1 1 1 2 2 4 4 8 8 4 8 0 0 0 0 0 4 8]';

%EDF.GDFTYP.TEXT={'char','int8','uint8','int16','uint16','int32','uint32','int64','uint64','float32','float64'};
%GDFTYP_BYTE=[1 1 1 2 2 4 4 8 8 4 8 0 0 0 0 0 4 8]';
%GDFTYPES=[0 1 2 3 4 5 6 7 16 17];

EDF.ErrNo = 0;

%%%%%%% ============= READ ===========%%%%%%%%%%%%
if (strcmp(arg2,'r') | strcmp(arg2,'r+')) 

[EDF.FILE.FID,MESSAGE]=fopen(FILENAME,arg2,'ieee-le');          
%EDF.FILE.FID=fid;

EDF.FILE.stderr=2;
if EDF.FILE.FID<0 
        %fprintf(EDF.FILE.stderr,'Error GDFOPEN: %s %s\n',MESSAGE,FILENAME);  
        H1=MESSAGE; H2=FILENAME;
        EDF.ErrNo = [32,EDF.ErrNo];
	return;
end;
if (arg2=='r') 
        EDF.FILE.OPEN = 1;
elseif (arg2=='r+') 
        EDF.FILE.OPEN = 2;
end;
EDF.FileName = FILENAME;

PPos=min([max(find(FILENAME=='.')) length(FILENAME)+1]);
SPos=max([0 find(FILENAME==filesep)]);
EDF.FILE.Ext = FILENAME(PPos+1:length(FILENAME));
EDF.FILE.Name = FILENAME(SPos+1:PPos-1);
if SPos==0
	EDF.FILE.Path = pwd;
else
	EDF.FILE.Path = FILENAME(1:SPos-1);
end;
EDF.FileName = [EDF.FILE.Path filesep EDF.FILE.Name '.' EDF.FILE.Ext];

%%% Read Fixed Header %%%
[tmp,count]=fread(EDF.FILE.FID,184,'uchar');     %
if count<184
        EDF.ErrNo = [64,EDF.ErrNo];
        return;
end;
H1=setstr(tmp');     %

EDF.VERSION=H1(1:8);                     % 8 Byte  Versionsnummer 
if ~(strcmp(EDF.VERSION,'0       ') | strcmp(EDF.VERSION(1:3),'GDF'))
        EDF.ErrNo = [1,EDF.ErrNo];
	if ~strcmp(EDF.VERSION(1:3),'   '); % if not a scoring file, 
%	    return; 
	end;
end;
EDF.PID = deblank(H1(9:88));                  % 80 Byte local patient identification
EDF.RID = deblank(H1(89:168));                % 80 Byte local recording identification

IsGDF=strcmp(EDF.VERSION(1:3),'GDF');

if IsGDF
        %EDF.T0=[str2num(H1(168+[1:4])) str2num(H1(168+[5 6])) str2num(H1(168+[7 8])) str2num(H1(168+[9 10])) str2num(H1(168+[11 12])) str2num(H1(168+[13:16]))/100 ];
        if 1, % if strcmp(EDF.VERSION(4:8),' 0.12'); 
	% Note: in future versions the date format might change. 
      		EDF.T0(1) = str2num( H1(168 + [ 1:4]));
		EDF.T0(2) = str2num( H1(168 + [ 5 6]));
        	EDF.T0(3) = str2num( H1(168 + [ 7 8]));
        	EDF.T0(4) = str2num( H1(168 + [ 9 10]));
        	EDF.T0(5) = str2num( H1(168 + [11 12]));
        	EDF.T0(6) = str2num( H1(168 + [13:16]))/100;
     	end; 
     
	if str2num(EDF.VERSION(4:8))<0.12
                tmp = setstr(fread(EDF.FILE.FID,8,'uchar')');    % 8 Byte  Length of Header
                EDF.HeadLen = str2num(tmp);    % 8 Byte  Length of Header
        else
                EDF.HeadLen = fread(EDF.FILE.FID,1,'int64');    % 8 Byte  Length of Header
        end;
        EDF.reserved1 = fread(EDF.FILE.FID,8+8+8+20,'uchar');     % 44 Byte reserved
        
        EDF.NRec = fread(EDF.FILE.FID,1,'int64');     % 8 Byte  # of data records
        if strcmp(EDF.VERSION(4:8),' 0.10')
                EDF.Dur =  fread(EDF.FILE.FID,1,'float64');    % 8 Byte  # duration of data record in sec
        else
                tmp =  fread(EDF.FILE.FID,2,'uint32');    % 8 Byte  # duration of data record in sec
                EDF.Dur =  tmp(1)./tmp(2);
        end;
        EDF.NS =   fread(EDF.FILE.FID,1,'uint32');     % 4 Byte  # of signals
else 
        if exist('OCTAVE_VERSION')>=5
                tmp=(find((toascii(H1)<32) | (toascii(H1)>126))); 	%%% snytax for OCTAVE
        else
                tmp=(find((H1<32) | (H1>126))); 		%%% syntax for Matlab
        end;
        if ~isempty(tmp) %%%%% not EDF because filled out with ASCII(0) - should be spaces
                %H1(tmp)=32; 
                EDF.ErrNo=[1025,EDF.ErrNo];
        end;
        %EDF.T0=[str2num(H1(168+[7 8])) str2num(H1(168+[4 5])) str2num(H1(168+[1 2])) str2num(H1(168+[9 10])) str2num(H1(168+[12 13])) str2num(H1(168+[15 16])) ];
        
        EDF.T0 = zeros(1,6);
        ErrT0=0;
        tmp = str2num( H1(168 + [ 7  8]));
        if ~isempty(tmp), EDF.T0(1) = tmp; else ErrT0 = 1; end;
        tmp = str2num( H1(168 + [ 4  5]));
        if ~isempty(tmp), EDF.T0(2) = tmp; else ErrT0 = 1; end;
        tmp = str2num( H1(168 + [ 1  2]));
        if ~isempty(tmp), EDF.T0(3) = tmp; else ErrT0 = 1; end;
        tmp = str2num( H1(168 + [ 9 10]));
        if ~isempty(tmp), EDF.T0(4) = tmp; else ErrT0 = 1; end;
        tmp = str2num( H1(168 + [12 13]));
        if ~isempty(tmp), EDF.T0(5) = tmp; else ErrT0 = 1; end;
        tmp = str2num( H1(168 + [15 16]));
        if ~isempty(tmp), EDF.T0(6) = tmp; else ErrT0 = 1; end;
        
	if any(EDF.T0~=fix(EDF.T0)); ErrT0=1; end;

        if ErrT0,
                ErrT0=0;
                EDF.ErrNo = [1032,EDF.ErrNo];
                
                tmp = H1(168 + [1:8]);
                for k = [3 2 1],
                        %fprintf(1,'\n zz%szz \n',tmp);
                        [tmp1,tmp] = strtok(tmp,' :./-');
			tmp1 = str2num([tmp1,' ']);
			
                        if isempty(tmp1)
                                ErrT0 = ErrT0 | 1;
                        else
                                EDF.T0(k)  = tmp1;
                        end;
                end;
                tmp = H1(168 + [9:16]);
                for k = [4 5 6],
                        [tmp1,tmp] = strtok(tmp,' :./-');
                        tmp1=str2num([tmp1,' ']);
                        if isempty(tmp1)
                                ErrT0 = ErrT0 | 1;
                        else
                                EDF.T0(k)  = tmp1;
                        end;
                end;
                if ErrT0
                        EDF.ErrNo = [2,EDF.ErrNo];
                end;
        else
                % Y2K compatibility until year 2084
                if EDF.T0(1) < 85    % for biomedical data recorded in the 1950's and converted to EDF
                        EDF.T0(1) = 2000+EDF.T0(1);
                elseif EDF.T0(1) < 100
                        EDF.T0(1) = 1900+EDF.T0(1);
                %else % already corrected, do not change
                end;
        end;     
        H1(185:256)=setstr(fread(EDF.FILE.FID,256-184,'uchar')');     %
        EDF.HeadLen = str2num(H1(185:192));           % 8 Bytes  Length of Header
        EDF.reserved1=H1(193:136);              % 44 Bytes reserved   
        EDF.NRec    = str2num(H1(237:244));     % 8 Bytes  # of data records
        EDF.Dur     = str2num(H1(245:252));     % 8 Bytes  # duration of data record in sec
        EDF.NS      = str2num(H1(253:256));     % 4 Bytes  # of signals
	EDF.AS.H1   = H1;	                     % for debugging the EDF Header
end;

if isempty(EDF.NS) %%%%% not EDF because filled out with ASCII(0) - should be spaces
        fprintf(EDF.FILE.stderr, 'Warning SDFOPEN: invalid NS-value in header of %s\n',EDF.FileName);
        EDF.ErrNo=[1040,EDF.ErrNo];
        EDF.NS=1;
end;

if isempty(EDF.HeadLen) %%%%% not EDF because filled out with ASCII(0) - should be spaces
        EDF.ErrNo=[1056,EDF.ErrNo];
        EDF.HeadLen=256*(1+EDF.NS);
end;


if any(~isempty(EDF.reserved1)) %%%%% not EDF because filled out with ASCII(0) - should be spaces
        fprintf(2, 'Warning SDFOPEN: 44bytes-reserved-field of fixed header is not empty in %s\n',EDF.FILE.Name);
        %EDF.ErrNo=[1057,EDF.ErrNo];
end;


if isempty(EDF.NRec) %%%%% not EDF because filled out with ASCII(0) - should be spaces
        EDF.ErrNo=[1027,EDF.ErrNo];
        EDF.NRec = -1;
end;

if isempty(EDF.Dur) %%%%% not EDF because filled out with ASCII(0) - should be spaces
        EDF.ErrNo=[1088,EDF.ErrNo];
        EDF.Dur=30;
else 
	if ~IsGDF,
		if (EDF.Dur>1) & (EDF.Dur~=fix(EDF.Dur)), % Blockduration must be in seconds or subsecond range 
		        EDF.ErrNo=[1088,EDF.ErrNo];
		end;
	end;	
end;

if  any(EDF.T0>[2084 12 31 24 59 59]) | any(EDF.T0<[1985 1 1 0 0 0])
        EDF.ErrNo = [4, EDF.ErrNo];
end;

%%% Read variable Header %%%
if ~IsGDF
        idx1=cumsum([0 H2idx]);
        idx2=EDF.NS*idx1;

        h2=zeros(EDF.NS,256);
        [H2,count]=fread(EDF.FILE.FID,EDF.NS*256,'uchar');
        if count < EDF.NS*256 
	        EDF.ErrNo=[8,EDF.ErrNo];
                return; 
        end;
                
        %tmp=find((H2<32) | (H2>126)); % would confirm 
        tmp = find((H2<32) | ((H2>126) & (H2~=255) & (H2~=181)& (H2~=230))); 
        if ~isempty(tmp) %%%%% not EDF because filled out with ASCII(0) - should be spaces
                H2(tmp) = 32; 
	        EDF.ErrNo = [1026,EDF.ErrNo];
        end;
        
        for k=1:length(H2idx);
                %disp([k size(H2) idx2(k) idx2(k+1) H2idx(k)]);
                h2(:,idx1(k)+1:idx1(k+1))=reshape(H2(idx2(k)+1:idx2(k+1)),H2idx(k),EDF.NS)';
        end;
        %size(h2),
        h2=setstr(h2);
        %(h2(:,idx1(9)+1:idx1(10))),
        %abs(h2(:,idx1(9)+1:idx1(10))),
        
        EDF.Label      =         h2(:,idx1(1)+1:idx1(2));
        EDF.Transducer =         h2(:,idx1(2)+1:idx1(3));
        EDF.PhysDim    =         h2(:,idx1(3)+1:idx1(4));
        EDF.PhysMin    = str2num(h2(:,idx1(4)+1:idx1(5)));
        EDF.PhysMax    = str2num(h2(:,idx1(5)+1:idx1(6)));
        EDF.DigMin     = str2num(h2(:,idx1(6)+1:idx1(7)));
        EDF.DigMax     = str2num(h2(:,idx1(7)+1:idx1(8)));
        EDF.PreFilt    =         h2(:,idx1(8)+1:idx1(9));
        EDF.SPR        = str2num(h2(:,idx1(9)+1:idx1(10)));
        %EDF.reserved  =       h2(:,idx1(10)+1:idx1(11));
        EDF.GDFTYP     = 3*ones(1,EDF.NS);	%	datatype

        if isempty(EDF.SPR), 
                fprintf(EDF.FILE.stderr, 'Warning SDFOPEN: invalid SPR-value in header of %s\n',EDF.FileName);
                EDF.SPR=ones(EDF.NS,1);
	        EDF.ErrNo=[1028,EDF.ErrNo];
        end;
else
        fseek(EDF.FILE.FID,256,'bof');
        EDF.Label      =  setstr(fread(EDF.FILE.FID,[16,EDF.NS],'uchar')');		
        EDF.Transducer =  setstr(fread(EDF.FILE.FID,[80,EDF.NS],'uchar')');	
        EDF.PhysDim    =  setstr(fread(EDF.FILE.FID,[ 8,EDF.NS],'uchar')');
%       EDF.AS.GDF.TEXT = EDF.GDFTYP.TEXT;
        EDF.PhysMin    =         fread(EDF.FILE.FID,[EDF.NS,1],'float64');	
        EDF.PhysMax    =         fread(EDF.FILE.FID,[EDF.NS,1],'float64');	
        EDF.DigMin     =         fread(EDF.FILE.FID,[EDF.NS,1],'int64');	
        EDF.DigMax     =         fread(EDF.FILE.FID,[EDF.NS,1],'int64');	
        
        EDF.PreFilt    =  setstr(fread(EDF.FILE.FID,[80,EDF.NS],'uchar')');	%	
        EDF.SPR        =         fread(EDF.FILE.FID,[ 1,EDF.NS],'uint32')';	%	samples per data record
        EDF.GDFTYP     =         fread(EDF.FILE.FID,[ 1,EDF.NS],'uint32');	%	datatype
        %                        fread(EDF.FILE.FID,[32,EDF.NS],'uchar')';	%	datatype
end;

%		EDF=gdfcheck(EDF,1);
if any(EDF.PhysMax==EDF.PhysMin), EDF.ErrNo=[1029,EDF.ErrNo]; end;	
if any(EDF.DigMax ==EDF.DigMin ), EDF.ErrNo=[1030,EDF.ErrNo]; end;	
EDF.Cal = (EDF.PhysMax-EDF.PhysMin)./(EDF.DigMax-EDF.DigMin);
EDF.Off = EDF.PhysMin - EDF.Cal .* EDF.DigMin;
EDF.Calib=[EDF.Off';(diag(EDF.Cal))];
EDF.SampleRate = EDF.SPR / EDF.Dur;

EDF.AS.spb = sum(EDF.SPR);	% Samples per Block
EDF.AS.bi = [0;cumsum(EDF.SPR)]; 
EDF.AS.BPR  = ceil(EDF.SPR.*GDFTYP_BYTE(EDF.GDFTYP+1)'); 
EDF.AS.SAMECHANTYP = all(EDF.AS.BPR == (EDF.SPR.*GDFTYP_BYTE(EDF.GDFTYP+1)')) & all(EDF.GDFTYP(:)==EDF.GDFTYP(1));
EDF.AS.GDFbi = [0;cumsum(ceil(EDF.SPR.*GDFTYP_BYTE(EDF.GDFTYP+1)'))]; 
EDF.AS.bpb = sum(ceil(EDF.SPR.*GDFTYP_BYTE(EDF.GDFTYP+1)'));	% Bytes per Block
EDF.AS.startrec = 0;
EDF.AS.numrec = 0;
EDF.FILE.POS = 0;

status = fseek(EDF.FILE.FID, 0, 'eof');
EDF.AS.endpos = ftell(EDF.FILE.FID);
fseek(EDF.FILE.FID, EDF.HeadLen, 'bof');

%[status EDF.AS.endpos EDF.HeadLen EDF.AS.bpb EDF.NRec]
if EDF.NRec == -1   % unknown record size, determine correct NRec
        EDF.NRec = floor((EDF.AS.endpos - EDF.HeadLen) / EDF.AS.bpb);
elseif  EDF.NRec ~= ((EDF.AS.endpos - EDF.HeadLen) / EDF.AS.bpb);
        EDF.ErrNo=[16,EDF.ErrNo];
%        EDF.NRec = floor((EDF.AS.endpos - EDF.HeadLen) / EDF.AS.bpb);
end; 

% if Channelselect, ReReferenzing and Resampling
% Overflowcheck, Adaptive FIR
% Layer 4 

%if nargin<3 %%%%%          Channel Selection 
if nargin <3 %else
        arg3=0;
end;

EDF.SIE.ChanSelect = 1:EDF.NS;
EDF.SIE.InChanSelect = 1:EDF.NS;

EDF.SIE.RR=1;
EDF.SIE.RS=0; %exist('arg5')==1; if EDF.SIE.RS, EDF.SIE.RS==(arg5>0); end;
EDF.SIE.TH=0; %(exist('arg6')==1);
EDF.SIE.RAW=0;
EDF.SIE.REGC=0;
EDF.SIE.TECG=0;
EDF.SIE.AFIR=0;
EDF.SIE.FILT=0;
EDF.SIE.TimeUnits_Seconds=1;
EDF.AS.MAXSPR=max(EDF.SPR(EDF.SIE.ChanSelect)); % Layer 3 defines EDF.AS.MAXSPR in GDFREAD

EDF.SIE.ReRefMx=eye(EDF.NS);
EDF.SIE.REG=eye(EDF.NS);
EDF.SIE.ReRefMx = EDF.SIE.ReRefMx(:,EDF.SIE.ChanSelect);
EDF.Calib=EDF.Calib*EDF.SIE.REG*EDF.SIE.ReRefMx; 

%if nargin>2 %%%%%          Channel Selection 
        EDF.SIE.REG=eye(EDF.NS);
        if arg3==0
                EDF.SIE.ChanSelect = 1:EDF.NS;
                EDF.SIE.InChanSelect = 1:EDF.NS;
                EDF.SIE.ReRefMx = eye(EDF.NS);
        else
                [nr,nc]=size(arg3);
                if all([nr,nc]>1), % Re-referencing
                        EDF.SIE.ReRefMx = [[arg3 zeros(size(arg3,1),EDF.NS-size(arg3,2))]; zeros(EDF.NS-size(arg3,1),EDF.NS)];
                        EDF.SIE.InChanSelect = find(any(EDF.SIE.ReRefMx'));
                        EDF.SIE.ChanSelect = find(any(EDF.SIE.ReRefMx));
                        if nargin>3
                        %        fprintf(EDF.FILE.stderr,'Error SDFOPEN: Rereferenzing does not work correctly with %s (more than 3 input arguments)\n',arg4);
                        end;
                else
                        EDF.SIE.ChanSelect = arg3; %full(sparse(1,arg3(:),1,1,EDF.NS));
                        EDF.SIE.InChanSelect = arg3;
                        EDF.SIE.ReRefMx = eye(EDF.NS);
                end;
        end;
        
        if (exist('sedfchk')==2),  
                EDF=sedfchk(EDF); % corrects incorrect Header information. 
        end;
        
	%EDF.SIE.CS=1;
	EDF.SIE.RR=1;
        EDF.SIE.RS=exist('arg5')==1; if EDF.SIE.RS, EDF.SIE.RS==(arg5>0); end;
        EDF.SIE.TH=0; %(exist('arg6')==1);
        EDF.SIE.RAW=0;
        EDF.SIE.REGC=0;
        EDF.SIE.TECG=0;
        EDF.SIE.AFIR=0;
        EDF.SIE.FILT=0;
        EDF.AS.MAXSPR=max(EDF.SPR(EDF.SIE.ChanSelect)); % Layer 3 defines EDF.AS.MAXSPR in GDFREAD
        
        EDF.SIE.REG=eye(EDF.NS);

        %elseif nargin>3   %%%%% RE-REFERENCING
        if nargin>3   
                if ~isempty(findstr(upper(arg4),'SIESTA'))
                        EDF.SIE.ReRefMx = eye(EDF.NS);% speye(EDF.NS);
                        EDF.SIE.ReRefMx(7,1:6)=[1 1 1 -1 -1 -1]/2;
                        EDF.SIE.TH=1;
	                % calculate (In)ChanSelect based on ReRefMatrix
    		        EDF.SIE.InChanSelect = find(any(EDF.SIE.ReRefMx'));
            		EDF.SIE.ChanSelect = find(any(EDF.SIE.ReRefMx));
                end;
                if ~isempty(findstr(upper(arg4),'EOG'))
                        tmp=findstr(upper(arg4),'EOG');
                        tmp=lower(strtok(arg4(tmp:length(arg4)),' +'));
                        EDF.SIE.ReRefMx = sparse(EDF.NS,EDF.NS);% speye(EDF.NS);
                        if any(tmp=='h'); EDF.SIE.ReRefMx(8:9,1)=[ 1 -1 ]';  end;
                        if any(tmp=='v'); EDF.SIE.ReRefMx(1:9,2)=[ 1 0 0 1 0 0 1 -1 -1]';end;
                        if any(tmp=='r'); EDF.SIE.ReRefMx(3:9,3)=[ 1 0 0 1 2 -2 -2]'/2;  end;
                        %EDF.SIE.TH=1;
	                % calculate (In)ChanSelect based on ReRefMatrix
    		        EDF.SIE.InChanSelect = find(any(EDF.SIE.ReRefMx'));
            		EDF.SIE.ChanSelect = find(any(EDF.SIE.ReRefMx));
                end;
                
                if ~isempty(findstr(upper(arg4),'UCAL'))
                        %EDF.SIE.ReRefMx = speye(EDF.NS); % OVERRIDES 'SIESTA' and 'REGRESS_ECG'
                        EDF.SIE.RR = 0;
                end;
                if ~isempty(findstr(upper(arg4),'RAW'))
                        EDF.SIE.RAW = 1;
                end;
                if ~isempty(findstr(upper(arg4),'OVERFLOW'))
                        EDF.SIE.TH = 1;
                end;
                if ~isempty(findstr(upper(arg4),'FailingElectrodeDetector'))
                        EDF.SIE.FED = 1;
			EDF.SIE.TH = 2; 
                end;
		
                if ~isempty(findstr(upper(arg4),'ECG')), % identify ECG channel for some ECG artifact processing method   
                        if ~isempty(findstr(upper(arg4),'SIESTA'))
                                channel1=12;
                                %channel2=1:9;
                                M=zeros(EDF.NS,1);M(channel1)=1;
                        end
                        if isfield(EDF,'ChanTyp')
                                M=upper(char(EDF.ChanTyp(channel1)))=='C';
                                channel1=find(M);
                                M=(upper(char(EDF.ChanTyp))=='E' | upper(char(EDF.ChanTyp))=='O' );
                                %channel2=find(M);
                        else
                                channel1=12;
                                %channel2=1:9;
                                M=zeros(EDF.NS,1);M(channel1)=1;
                        end;
                        
                end;
                if ~isempty(findstr(upper(arg4),'RECG'))
                        EDF.SIE.REGC = 1;
                        if all(EDF.SIE.InChanSelect~=channel1)
                                EDF.SIE.InChanSelect=[EDF.SIE.InChanSelect channel1];        
                        end;
                end;
                
                if ~isempty(findstr(upper(arg4),'TECG'))
                        fprintf(EDF.FILE.stderr,'SDFOPEN: status of TECG Mode: alpha test passed\n');    
                        %%%% TECG - ToDo
                        % - optimize window
                        if exist([lower(EDF.FILE.Name) 'ECG.mat']);
                                if exist('OCTAVE_VERSION')==5
                                        load(file_in_loadpath([lower(EDF.FILE.Name) 'ECG.mat']));
                                else
                                        load([lower(EDF.FILE.Name) 'ECG.mat']);
                                end;
                                if isstruct(QRS)
                                        EDF.SIE.TECG = 1;
					%%%  EDF.AS.MAXSPR=size(QRS.Templates,1)/3;
                                        EDF.AS.MAXSPR=max(EDF.AS.MAXSPR,EDF.SPR(channel1));
				else
                                        fprintf(EDF.FILE.stderr,'WARNING SDFOPEN: %s invalid for TECG\n',[ lower(EDF.FILE.Name) 'ECG.mat']);
                                end;
                        else
                                fprintf(EDF.FILE.stderr,'WARNING SDFOPEN: %s not found\t (needed for Mode=TECG)\n',[lower(EDF.FILE.Name) 'ECG.mat']);
                        end;
                end;
                
                if ~isempty(findstr(upper(arg4),'AFIR')) 
                        % Implements Adaptive FIR filtering for ECG removal in EDF/GDF-tb.
                        % based on the Algorithm of Mikko Koivuluoma <k7320@cs.tut.fi>
                                
                        % channel2 determines which channels should be corrected with AFIR 
                        
                        %EDF = sdf_afir(EDF,12,channel2);
                        channel2=EDF.SIE.ChanSelect;	
                        fprintf(EDF.FILE.stderr,'Warning SDFOPEN: option AFIR still buggy\n');    
                        if isempty(find(channel2))
                                EDF.SIE.AFIR=0;
                        else
                                EDF.SIE.AFIR=1;
                                EDF.AFIR.channel2=channel2;
                                
                                EDF.AFIR.alfa=0.01;
                                EDF.AFIR.gamma=1e-32;
                                
                                EDF.AFIR.delay  = ceil(0.05*EDF.AS.MAXSPR/EDF.Dur); 
                                EDF.AFIR.nord = EDF.AFIR.delay+EDF.AS.MAXSPR/EDF.Dur; 
                                
                                EDF.AFIR.nC = length(EDF.AFIR.channel2);
                                EDF.AFIR.w = zeros(EDF.AFIR.nC, max(EDF.AFIR.nord));
                                EDF.AFIR.x = zeros(1, EDF.AFIR.nord);
                                EDF.AFIR.d = zeros(EDF.AFIR.delay, EDF.AFIR.nC);
                                
                                channel1=12;
                                
                                if isfield(EDF,'ChanTyp')
                                        if upper(char(EDF.ChanTyp(channel1)))=='C';
                                                EDF.AFIR.channel1 = channel1;
                                        else
                                                EDF.AFIR.channel1 = find(EDF.ChanTyp=='C');
                                                fprintf(EDF.FILE.stderr,'Warning %s: #%i is not an ECG channel, %i used instead\n' ,filename,channel1,EDF.AFIR.channel1);
                                        end;
                                else
                                        EDF.AFIR.channel1 = channel1;
                                end;
                                
                                if all(EDF.SIE.InChanSelect~=channel1)
                                        EDF.SIE.InChanSelect=[EDF.SIE.InChanSelect channel1];        
                                end;
                        end;
                end;


                if isempty(findstr(upper(arg4),'NOTCH50')) 
                        EDF.Filter.A=1;
                        EDF.Filter.B=1;
                else
                        EDF.SIE.FILT=1;
                        EDF.Filter.A=1;
                        EDF.Filter.B=1;
                        %if all(EDF.SampleRate(EDF.SIE.ChanSelect)==100)
                        if EDF.AS.MAXSPR/EDF.Dur==100
                                EDF.Filter.B=[1 1]/2;
                                %elseif all(EDF.SampleRate(EDF.SIE.ChanSelect)==200)
                        elseif EDF.AS.MAXSPR/EDF.Dur==200
                                EDF.Filter.B=[1 1 1 1]/4;
                                %elseif all(EDF.SampleRate(EDF.SIE.ChanSelect)==256)
                        elseif EDF.AS.MAXSPR/EDF.Dur==400
                                EDF.Filter.B=ones(1,8)/8;
                                %elseif all(EDF.SampleRate(EDF.SIE.ChanSelect)==256)
                        elseif EDF.AS.MAXSPR/EDF.Dur==256
                                EDF.Filter.B=poly([exp([-1 1 -2 2]*2*pi*i*50/256)]); %max(EDF.SampleRate(EDF.SIE.ChanSelect)))]);
                                EDF.Filter.B=EDF.Filter.B/sum(EDF.Filter.B);
                        else
                                fprintf(EDF.FILE.stderr,'Warning SDFOPEN: 50Hz Notch does not fit\n');
                        end;
                end;



                if ~isempty(findstr(upper(arg4),'NOTCH60')) 
                        fprintf(EDF.FILE.stderr,'Warning SDFOPEN: option NOTCH60 not implemented yet.\n');    
                end;


                if ~isempty(findstr(upper(arg4),'HPF')),  % high pass filtering
                        if EDF.SIE.FILT==0; EDF.Filter.B=1; end;
                        EDF.SIE.FILT=1;
                        EDF.Filter.A=1;
		end;
                if ~isempty(findstr(upper(arg4),'HPF0.0Hz')),  % high pass filtering
                        EDF.Filter.B=conv([1 -1], EDF.Filter.B);
                elseif ~isempty(findstr(upper(arg4),'TAU')),  % high pass filtering / compensate time constant
                        tmp=findstr(upper(arg4),'TAU');
                        TAU=strtok(upper(arg4(tmp:length(arg4))),'S');
                        tau=str2num(TAU);
                        if isempty(tau)
                                fprintf(EDF.FILE.stderr,'Warning SDFOPEN: invalid tau-value.\n');
                        else
                                EDF.Filter.B=conv([1 (EDF.Dur/EDF.AS.MAXSPR/tau-1)], EDF.Filter.B);
                        end;
			
                %%%% example 'HPF_1.0Hz_Hamming',  % high pass filtering
                elseif ~isempty(findstr(upper(arg4),'HPF')),  % high pass filtering
			    tmp=findstr(upper(arg4),'HPF');
			    FilterArg0=arg4(tmp+4:length(arg4));
			    %[tmp,FilterArg0]=strtok(arg4,'_');
			    [FilterArg1,FilterArg2]=strtok(FilterArg0,'_');
			    [FilterArg2,FilterArg3]=strtok(FilterArg2,'_');
			    tmp=findstr(FilterArg1,'Hz');
			    F0=str2num(FilterArg1(1:tmp-1));				    
			    B=feval(FilterArg2,F0*EDF.AS.MAXSPR/EDF.Dur);
			    B=B/sum(B);
			    B(ceil(length(B)/2))=(B(ceil(length(B)/2)))-1;
			    
                        EDF.Filter.B=conv(-B, EDF.Filter.B);
                end;


                if ~isempty(findstr(upper(arg4),'UNITS_BLOCK'))
			EDF.SIE.TimeUnits_Seconds=0; 
                end;

        end; % end nargin >3
        
        if EDF.SIE.FILT==1;
		EDF.Filter.Z=[];
                for k=1:length(EDF.SIE.ChanSelect),
                        [tmp,EDF.Filter.Z(:,k)]=filter(EDF.Filter.B,EDF.Filter.A,zeros(length(EDF.Filter.B+1),1));
                end;
    		EDF.FilterOVG.Z=EDF.Filter.Z;
        end;
        
        if EDF.SIE.REGC
                FN=[lower(EDF.FILE.Name) 'cov.mat'];
                if exist(FN)~=2
                        fprintf(EDF.FILE.stderr,'Warning %s: Covariance-file %s not found.\n',EDF.AS.Method,FN);
                        EDF.SIE.REGC=0;   
                else
                        if exist('OCTAVE_VERSION')==5
                                load(file_in_loadpath(FN));
                        else
                                load(FN);
                        end;
                        if exist('XC') == 1
                                %EDF.SIE.COV = tmp.XC;
                                %[N,MU,COV,Corr]=decovm(XC);
                                N=size(XC,2);
                                COV=(XC(2:N,2:N)/XC(1,1)-XC(2:N,1)*XC(1,2:N)/XC(1,1)^2);
                                
                                %clear tmp;
                                %cov = diag(EDF.Cal)*COV*diag(EDF.Cal);
                                mcov = M'*diag(EDF.Cal)*COV*diag(EDF.Cal);
                                %mcov(~())=0;
                                EDF.SIE.REG = eye(EDF.NS) - M*((mcov*M)\(mcov));
                                EDF.SIE.REG(channel1,channel1) = 1; % do not remove the regressed channels
                                %mcov, EDF.SIE.REG, 
                        else
                                fprintf(EDF.FILE.stderr,'Error %s: Regression Coefficients for ECG minimization not available.\n',EDF.AS.Method);
                        end;
                end;
                
        end;
        
        if EDF.SIE.TECG == 1; 
                % define channels that should be corrected
                if isfield(QRS,'Version')
                        OutChanSelect=[1:11 13:EDF.NS];
                        if EDF.SIE.REGC % correct templates
                                QRS.Templates=QRS.Templates*EDF.SIE.REG;
                                fprintf(EDF.FILE.stderr,'Warning SDFOPEN: Mode TECG+RECG not tested\n');
                        end;
                        if QRS.Version~=2
                                fprintf(EDF.FILE.stderr,'Warning SDFOPEN Mode TECG: undefined QRS-version\n');
                        end;
                else
                        %OutChanSelect=find(EDF.ChanTyp=='E' | EDF.ChanTyp=='O');
                        OutChanSelect=[1:9 ];
                        if any(EDF.SIE.ChanSelect>10)
                                fprintf(EDF.FILE.stderr,'Warning SDFOPEN: Mode TECG: Only #1-#9 are corrected\n');
                        end;
                        if EDF.SIE.REGC, % correct the templates
                                QRS.Templates=QRS.Templates*EDF.SIE.REG([1:9 12],[1:9 12]);
                                fprintf(EDF.FILE.stderr,'Warning SDFOPEN: Mode TECG+RECG not tested\n');
                        end;
                        
                end;
                fs = EDF.SPR(12)/EDF.Dur; 
                QRS.Templates=detrend(QRS.Templates,0); %remove mean
                EDF.TECG.idx = [(QRS.Index-fs/2-1) (EDF.NRec+1)*EDF.SPR(12)]; %include terminating element
                EDF.TECG.idxidx = 1; %pointer to next index

                % initialize if any spike is detected before first window    
	        pulse = zeros(length(QRS.Templates),1);
    		Index=[];
	        while EDF.TECG.idx(EDF.TECG.idxidx) < 1,
    		        Index=[Index EDF.TECG.idx(EDF.TECG.idxidx)-EDF.AS.startrec*EDF.SPR(12)];
            		EDF.TECG.idxidx=EDF.TECG.idxidx+1;
	        end;
	        if ~isempty(Index)
    		        pulse(Index+length(QRS.Templates)) = 1;  
	        end;
        
                for i=1:length(EDF.SIE.InChanSelect),
                        k=find(OutChanSelect==EDF.SIE.InChanSelect(i));
                        if isempty(k)
                                EDF.TECG.QRStemp(:,i) = zeros(fs,1);
                        else
                                EDF.TECG.QRStemp(:,i) = QRS.Templates(0.5*fs:1.5*fs-1,k).*hanning(fs);
                        end;
                        [tmp,EDF.TECG.Z(:,i)] = filter(EDF.TECG.QRStemp(:,i),1,pulse);
                end;
        end; % if EDF.SIE.TECG==1
        
        %syms Fp1 Fp2 M1 M2 O2 O1 A1 A2 C3 C4
        for k=EDF.SIE.ChanSelect,
                %fprintf(1,'#%i: ',k);
                tmp=find(EDF.SIE.ReRefMx(:,k))';
                
                if EDF.SIE.ReRefMx(tmp(1),k)==1,
                        x=EDF.Label(tmp(1),:);
                else
                        x=sprintf('%3.1f*%s',EDF.SIE.ReRefMx(tmp(1),k),deblank(EDF.Label(tmp(1),:)));
                end;
                for l=2:length(tmp), L=tmp(l);
                        if (EDF.SPR(tmp(l-1),:)~=EDF.SPR(tmp(l),:))  
                                fprintf(EDF.FILE.stderr,'Warning %s: SampleRate Mismatch in "%s", channel #%i and #%i\n',upper(EDF.AS.Method),EDF.FILE.Name,tmp(l-1),tmp(l));
                        end;
                        if ~strcmp(EDF.PhysDim(tmp(l-1),:),EDF.PhysDim(tmp(l),:))  
                                fprintf(EDF.FILE.stderr,'Warning %s: Dimension Mismatch in "%s", channel #%i and #%i\n',upper(EDF.AS.Method),EDF.FILE.Name,tmp(l-1),tmp(l));
                        end;
                        if ~strcmp(EDF.Transducer(tmp(l-1),:),EDF.Transducer(tmp(l),:))  
                                fprintf(EDF.FILE.stderr,'Warning %s: Transducer Mismatch in "%s", channel #%i and #%i\n',upper(EDF.AS.Method),EDF.FILE.Name,tmp(l-1),tmp(l));
                        end;
                        if ~strcmp(EDF.PreFilt(tmp(l-1),:),EDF.PreFilt(tmp(l),:))  
                                fprintf(EDF.FILE.stderr,'Warning %s: PreFiltering Mismatch in "%s", channel #%i and #%i\n',upper(EDF.AS.Method),EDF.FILE.Name,tmp(l-1),tmp(l));
                        end;
                        x=[x sprintf('+(%3.1f)*(%s)',EDF.SIE.ReRefMx(tmp(l),k),deblank(EDF.Label(tmp(l),:)))];
                end;
                EDF.Label(k,1:length(x))=x; %char(sym(x))
        end;
        %Label,
        EDF.SIE.ReRefMx = EDF.SIE.ReRefMx(:,EDF.SIE.ChanSelect);
	if exist('OCTAVE_VERSION')>=5,
            EDF.Calib = (EDF.Calib*EDF.SIE.REG*EDF.SIE.ReRefMx); % important to be sparse, otherwise overflow-check does not work correctly.
	    fprintf(EDF.FILE.stderr,'Warning SDFOPEN: Overflow Check does not work without SPARSE\n');
	else
            EDF.Calib = sparse(EDF.Calib*EDF.SIE.REG*EDF.SIE.ReRefMx); % important to be sparse, otherwise overflow-check does not work correctly.
	end;
	        
        if 1, % ??? not sure, whether it has any advantage
        EDF.PhysDim = EDF.PhysDim(EDF.SIE.ChanSelect,:);
        EDF.PreFilt = EDF.PreFilt(EDF.SIE.ChanSelect,:);
        EDF.Transducer = EDF.Transducer(EDF.SIE.ChanSelect,:);
        end;
        
        if EDF.SIE.RS,
		tmp = EDF.AS.MAXSPR/EDF.Dur;
                if arg5==0
                        %arg5 = max(EDF.SPR(EDF.SIE.ChanSelect))/EDF.Dur;
                        
                elseif ~rem(tmp,arg5); % The target sampling rate must divide the source sampling rate   
			EDF.SIE.RS = 1;
			tmp=tmp/arg5;
			EDF.SIE.T = ones(tmp,1)/tmp;
                elseif arg5==100; % Currently, only the target sampling rate of 100Hz are supported. 
                        EDF.SIE.RS=1;
                        tmp=EDF.AS.MAXSPR/EDF.Dur;
                        if exist('OCTAVE_VERSION')
                                load resample_matrix4octave.mat T256100 T200100;
                        else
                                load('resample_matrix');
                        end;
                        if 1,
                                if tmp==400,
                                        EDF.SIE.T=ones(4,1)/4;
                                elseif tmp==256,
                                        EDF.SIE.T=T256100;
                                elseif tmp==200,
                                        EDF.SIE.T=T200100; 
                                elseif tmp==100,
                                        EDF.SIE.T=1;
                                else
                                        fprintf('Warning %s-READ: sampling rates should be equal\n',upper(EDF.AS.Method));     
                                end;	
                        else
                                tmp=EDF.SPR(EDF.SIE.ChanSelect)/EDF.Dur;
                                if all((tmp==256) | (tmp<100)) 
                                        EDF.SIE.RS = 1;
                                        %tmp=load(RSMN,'T256100');	
                                        EDF.SIE.T = T256100;	
                                elseif all((tmp==400) | (tmp<100)) 
                                        EDF.SIE.RS = 1;
                                        EDF.SIE.T = ones(4,1)/4;
                                elseif all((tmp==200) | (tmp<100)) 
                                        EDF.SIE.RS = 1;
                                        %tmp=load(RSMN,'T200100');	
                                        EDF.SIE.T = T200100;	
                                elseif all(tmp==100) 
                                        %EDF.SIE.T=load('resample_matrix','T100100');	
                                        EDF.SIE.RS=0;
                                else
                                        EDF.SIE.RS=0;
                                        fprintf('Warning %s-READ: sampling rates should be equal\n',upper(EDF.AS.Method));     
                                end;
                        end;
                else
                        fprintf(EDF.FILE.stderr,'Error %s-READ: invalid target sampling rate of %i Hz\n',upper(EDF.AS.Method),arg5);
                        EDF.SIE.RS=0;
			EDF.ErrNo=[EDF.ErrNo,];
			%EDF=sdfclose(EDF);
			%return;
                end;
        end;
        
        FN=[lower(EDF.FILE.Name), 'th.mat'];
        if exist(FN)~=2,
	        if EDF.SIE.TH, % && ~exist('OCTAVE_VERSION'),
                        fprintf(EDF.FILE.stderr,'Warning %s: THRESHOLD-file %s not found.\n',EDF.AS.Method,FN);
                        EDF.SIE.TH=0;   
                end;
        else
                if exist('OCTAVE_VERSION')==5
                        tmp=load(file_in_loadpath(FN));
                else
                        tmp=load(FN);
                end;
    	        if isfield(tmp,'TRESHOLD') 
        	        EDF.SIE.THRESHOLD = tmp.TRESHOLD;
                %else
            	        %fprintf(EDF.FILE.stderr,'Error %s: TRESHOLD''s not found.\n',EDF.AS.Method);
                end;
        end;
    	if EDF.SIE.TH>1, % Failing electrode detector 
	        fprintf(2,'Warning SDFOPEN: FED not implemented yet\n');
                for k=1:length(InChanSelect),K=InChanSelect(k);
	        %for k=1:EDF.NS,
   % 		        [y1,EDF.Block.z1{k}] = filter([1 -1], 1, zeros(EDF.SPR(K)/EDF.Dur,1));
    %        		[y2,EDF.Block.z2{k}] = filter(ones(1,EDF.SPR(K)/EDF.Dur)/(EDF.SPR(K)/EDF.Dur),1,zeros(EDF.SPR(K)/EDF.Dur,1));
                
    %       		[y3,EDF.Block.z3{k}] = filter(ones(1,EDF.SPR(K)/EDF.Dur)/(EDF.SPR(K)/EDF.Dur),1,zeros(EDF.SPR(K)/EDF.Dur,1));
    		end;
        end;

% Initialization of Bufferblock for random access (without EDF-blocklimits) of data 
	if ~EDF.SIE.RAW & EDF.SIE.TimeUnits_Seconds
                EDF.Block.number=[0 0 0 0]; %Actual Blocknumber, start and end time of loaded block, diff(EDF.Block.number(1:2))==0 denotes no block is loaded;
                                            % EDF.Blcok.number(3:4) indicate start and end of the returned data, [units]=samples.
		EDF.Block.data=[];
		EDF.Block.dataOFCHK=[];
	end;
        
%end; % end of SDFOPEN-READ



%%%%%%% ============= WRITE ===========%%%%%%%%%%%%        

elseif (arg2=='w') | (arg2=='w+')
%        fprintf(EDF.FILE.stderr,'error EDFOPEN: write mode not possible.\n'); 
        H1=[]; H2=[];
%        return;
        
EDF.SIE.RAW = 0;
if ~isstruct(arg1)  % if arg1 is the filename 
        EDF.FileName=arg1;
        if nargin<3
                tmp=input('SDFOPEN: list of samplerates for each channel? '); 
                EDF.SampleRate = tmp(:);
        else
                EDF.SampleRate=arg3;
        end;
        EDF.NS=length(EDF.SampleRate);
        if nargin<4
                tmp=input('SDFOPEN: Duration of one block in seconds: '); 
                EDF.Dur = tmp;
                EDF.SPR=EDF.Dur*EDF.SampleRate;
        else
                if ~isempty(findstr(upper(arg4),'RAW'))
                        EDF.SIE.RAW = 1;
                else
                        EDF.Dur = arg4;
                        EDF.SPR=EDF.Dur*EDF.SampleRate;
                end;
        end;
end;

FILENAME=EDF.FileName;
if (arg2=='w') 
        [fid,MESSAGE]=fopen(FILENAME,'w','ieee-le');          
elseif (arg2=='w+')  % may be called only by SDFCLOSE
        if EDF.FILE.OPEN==2 
                [fid,MESSAGE]=fopen(FILENAME,'r+','ieee-le');          
        else
                fprintf(EDF.FILE.stderr,'Error SDFOPEN-W+: Cannot open %s for write access\n',FILENAME);
                return;
        end;
end;
if fid<0 
        %fprintf(EDF.FILE.stderr,'Error EDFOPEN: %s\n',MESSAGE);  
        H1=MESSAGE;H2=[];
        EDF.ErrNo = EDF.ErrNo + 32;
        fprintf(EDF.FILE.stderr,'Error SDFOPEN-W: Could not open %s \n',FILENAME);
        return;
end;
EDF.FILE.FID = fid;
EDF.FILE.OPEN = 2;

%%%% generate optional parameters

PPos=min([max(find(FILENAME=='.')) length(FILENAME)+1]);
SPos=max([0 find(FILENAME==filesep)]);
EDF.FILE.Ext = FILENAME(PPos+1:length(FILENAME));
EDF.FILE.Name = FILENAME(SPos+1:PPos-1);
if SPos==0
	EDF.FILE.Path = pwd;
else
	EDF.FILE.Path = FILENAME(1:SPos-1);
end;
EDF.FileName = [EDF.FILE.Path filesep EDF.FILE.Name '.' EDF.FILE.Ext];

% Check all fields of Header1
if ~isfield(EDF,'VERSION')
        fprintf('Warning SDFOPEN-W: EDF.VERSION not defined; default=EDF assumed\n');
        EDF.VERSION='0       '; % default EDF-format
        %EDF.ErrNo = EDF.ErrNo + 128;
        %fclose(EDF.FILE.FID); return;
end;

IsGDF=strcmp(upper(EDF.VERSION(1:3)),'GDF');
if ~IsGDF
        EDF.VERSION = '0       ';
        fprintf(EDF.FILE.stderr,'\nData are stored with integer16.\nMeasures for minimizing round-off errors have been taken.\nDespite, overflow and round off errors may occur.\n');  
        
        if sum(EDF.SPR)>61440/2;
                fprintf(EDF.FILE.stderr,'\nWarning SDFOPEN: One block exceeds 61440 bytes.\n')
        end;
else
        EDF.VERSION = 'GDF 0.12';
end;

if ~isfield(EDF,'PID')
        fprintf(EDF.FILE.stderr,'Warning SDFOPEN-W: EDF.PID not defined\n');
        EDF.PID=setstr(32+zeros(1,80));
end;
if ~isfield(EDF,'RID')
        fprintf(EDF.FILE.stderr,'Warning SDFOPEN-W: EDF.RID not defined\n');
        EDF.RID=setstr(32+zeros(1,80));
end;
if ~isfield(EDF,'T0')
        EDF.T0=zeros(1,6);
        fprintf(EDF.FILE.stderr,'Warning SDFOPEN-W: EDF.T0 not defined\n');
end;
if ~isfield(EDF,'reserved1')
        EDF.reserved1=char(ones(1,44)*32);
else
        tmp=min(8,size(EDF.reserved1,2));
        EDF.reserved1=[EDF.reserved1(1,1:tmp) 32+zeros(1,44-tmp)];
end;
if ~isfield(EDF,'NRec')
        EDF.NRec=-1;
end;
if ~isfield(EDF,'Dur')
        fprintf('Warning SDFOPEN-W: EDF.Dur not defined\n');
        EDF.Dur=NaN;
        EDF.ErrNo = EDF.ErrNo + 128;
        fclose(EDF.FILE.FID); return;
end;
if ~isfield(EDF,'NS')
        EDF.ERROR = sprintf('Error SDFOPEN-W: EDF.NS not defined\n');
        EDF.ErrNo = EDF.ErrNo + 128;
        fclose(EDF.FILE.FID); return;
end;

% Check all fields of Header2
if ~isfield(EDF,'Label')
        EDF.Label=setstr(32+zeros(EDF.NS,16));
        fprintf(EDF.FILE.stderr,'Warning SDFOPEN-W: EDF.Label not defined\n');
else
        tmp=min(16,size(EDF.Label,2));
        EDF.Label=[EDF.Label(1:EDF.NS,1:tmp) 32+zeros(EDF.NS,16-tmp)];
end;
if ~isfield(EDF,'Transducer')
        EDF.Transducer=setstr(32+zeros(EDF.NS,80));
else
        tmp=min(80,size(EDF.Transducer,2));
        EDF.Transducer=[EDF.Transducer(1:EDF.NS,1:tmp) 32+zeros(EDF.NS,80-tmp)];
end;
if ~isfield(EDF,'PreFilt')
        EDF.PreFilt=setstr(32+zeros(EDF.NS,80));
else
        tmp=min(80,size(EDF.PreFilt,2));
        EDF.PreFilt=[EDF.PreFilt(1:EDF.NS,1:tmp) 32+zeros(EDF.NS,80-tmp)];
end;
if ~isfield(EDF,'PhysDim')
        EDF.PhysDim=setstr(32+zeros(EDF.NS,8));
        fprintf(EDF.FILE.stderr,'Warning SDFOPEN-W: EDF.PhysDim not defined\n');
else
        tmp=min(8,size(EDF.PhysDim,2));
        EDF.PhysDim=[EDF.PhysDim(1:EDF.NS,1:tmp) 32+zeros(EDF.NS,8-tmp)];
end;

if ~isfield(EDF,'PhysMin')
        fprintf(EDF.FILE.stderr,'Warning SDFOPEN-W: EDF.PhysMin not defined\n');
        EDF.PhysMin=repmat(nan,EDF.NS,1);
        %EDF.ERROR = sprintf('Error SDFOPEN-W: EDF.PhysMax not defined\n');
        %EDF.ErrNo = EDF.ErrNo + 128;
        %fclose(EDF.FILE.FID); return;
else
        EDF.PhysMin=EDF.PhysMin(1:EDF.NS);
end;
if ~isfield(EDF,'PhysMax')
        fprintf('Warning SDFOPEN-W: EDF.PhysMax not defined\n');
        EDF.PhysMax=repmat(nan,EDF.NS,1);
        %EDF.ERROR = sprintf('Error SDFOPEN-W: EDF.PhysMax not defined\n');
        %EDF.ErrNo = EDF.ErrNo + 128;
        %fclose(EDF.FILE.FID); return;
else
        EDF.PhysMax=EDF.PhysMax(1:EDF.NS);
end;
if ~isfield(EDF,'DigMin')
        fprintf(EDF.FILE.stderr,'Warning SDFOPEN-W: EDF.DigMin not defined\n');
        EDF.DigMin=repmat(nan,EDF.NS,1);
        %EDF.ERROR = sprintf('Error SDFOPEN-W: EDF.DigMax not defined\n');
        %EDF.ErrNo = EDF.ErrNo + 128;
        %fclose(EDF.FILE.FID); return;
else
        EDF.DigMin=EDF.DigMin(1:EDF.NS);
end;
if ~isfield(EDF,'DigMax')
        fprintf('Warning SDFOPEN-W: EDF.DigMax not defined\n');
        EDF.DigMax=repmat(nan,EDF.NS,1);
        %EDF.ERROR = sprintf('Error SDFOPEN-W: EDF.DigMax not defined\n');
        %EDF.ErrNo = EDF.ErrNo + 128;
        %fclose(EDF.FILE.FID); return;
else
        EDF.DigMax=EDF.DigMax(1:EDF.NS);
end;
if ~isfield(EDF,'SPR')
        fprintf('Warning SDFOPEN-W: EDF.SPR not defined\n');
        EDF.SPR=repmat(nan,EDF.NS,1);
        EDF.ERROR = sprintf('Error SDFOPEN-W: EDF.SPR not defined\n');
        EDF.ErrNo = EDF.ErrNo + 128;
        fclose(EDF.FILE.FID); return;
else
        EDF.SPR=reshape(EDF.SPR(1:EDF.NS),EDF.NS,1);
end;

if IsGDF
        if ~isfield(EDF,'GDFTPY')
                EDF.ERROR = sprintf('Error SDFOPEN-W: EDF.GDFTYP not defined\n');
                EDF.ErrNo = EDF.ErrNo + 128;
                fclose(EDF.FILE.FID); return;
        else
                EDF.GDFTYP=EDF.GDFTYP(1:EDF.NS);
        end;
else
        EDF.GDFTYP=3+zeros(1,EDF.NS);
end;

%%%%%% generate Header 1, first 256 bytes 
EDF.HeadLen=(EDF.NS+1)*256;
H1=setstr(32*ones(1,256));
H1(1:8)=EDF.VERSION; %sprintf('%08i',EDF.VERSION);     % 8 Byte  Versionsnummer 
H1( 8+(1:length(EDF.PID)))=EDF.PID;       
H1(88+(1:length(EDF.RID)))=EDF.RID;
%H1(185:192)=sprintf('%-8i',EDF.HeadLen);

if IsGDF
        H1(168+(1:16))=sprintf('%04i%02i%02i%02i%02i%02i%02i',floor(EDF.T0),rem(EDF.T0(6),1));
        fwrite(fid,H1(1:184),'uchar');
        fwrite(fid,EDF.HeadLen,'int64');
        fwrite(fid,ones(8,1)*32,'byte'); % EP_ID=ones(8,1)*32;
        fwrite(fid,ones(8,1)*32,'byte'); % Lab_ID=ones(8,1)*32;
        fwrite(fid,ones(8,1)*32,'byte'); % T_ID=ones(8,1)*32;
        fwrite(fid,ones(20,1)*32,'byte'); % 
        fwrite(fid,EDF.NRec,'int64');
        %fwrite(fid,EDF.Dur,'float64');
        [n,d]=rat(EDF.Dur); fwrite(fid,[n d], 'uint32');
        fwrite(fid,EDF.NS,'uint32');
else
        H1(168+(1:16))=sprintf('%02i.%02i.%02i%02i:%02i:%02i',rem(EDF.T0([3 2 1 4 5 6]),100));
        H1(185:192)=sprintf('%-8i',EDF.HeadLen);
        H1(193:236)=EDF.reserved1;
        H1(237:244)=sprintf('%-8i',EDF.NRec);
        H1(245:252)=sprintf('%-8i',EDF.Dur);
        H1(253:256)=sprintf('%-4i',EDF.NS);
        H1(find(H1==0))=32;
        fwrite(fid,H1,'uchar');
end;        

%%%%%% generate Header 2,  NS*256 bytes 
if ~IsGDF;
        sPhysMin=32+zeros(EDF.NS,8);
        sPhysMax=32+zeros(EDF.NS,8);
        for k=1:EDF.NS,
                tmp=sprintf('%-8g',EDF.PhysMin(k));
                lt=length(tmp);
                if lt<9
                        sPhysMin(k,1:lt)=tmp;
                else
                        if any(upper(tmp)=='E') | find(tmp=='.')>8,
                                fprintf(EDF.FILE.stderr,'Error SDFOPEN-W: PhysMin(%i) does not fit into header\n', k);
                        else
                                sPhysMin(k,:)=tmp(1:8);
                        end;
                end;
                tmp=sprintf('%-8g',EDF.PhysMax(k));
                lt=length(tmp);
                if lt<9
                        sPhysMax(k,1:lt)=tmp;
                else
                        if any(upper(tmp)=='E') | find(tmp=='.')>8,
                                fprintf(EDF.FILE.stderr,'Error SDFOPEN-W: PhysMin(%i) does not fit into header\n', k);
                        else
                                sPhysMax(k,:)=tmp(1:8);
                        end;
                end;
        end;
        sPhysMin=reshape(sprintf('%-8.1f',EDF.PhysMin)',8,EDF.NS)';
        sPhysMax=reshape(sprintf('%-8.1f',EDF.PhysMax)',8,EDF.NS)';
        
        idx1=cumsum([0 H2idx]);
        idx2=EDF.NS*idx1;
        h2=setstr(32*ones(EDF.NS,256));
        size(h2);
        h2(:,idx1(1)+1:idx1(2))=EDF.Label;
        h2(:,idx1(2)+1:idx1(3))=EDF.Transducer;
        h2(:,idx1(3)+1:idx1(4))=EDF.PhysDim;
        %h2(:,idx1(4)+1:idx1(5))=sPhysMin;
        %h2(:,idx1(5)+1:idx1(6))=sPhysMax;
        h2(:,idx1(4)+1:idx1(5))=sPhysMin;
        h2(:,idx1(5)+1:idx1(6))=sPhysMax;
        h2(:,idx1(6)+1:idx1(7))=reshape(sprintf('%-8i',EDF.DigMin)',8,EDF.NS)';
        h2(:,idx1(7)+1:idx1(8))=reshape(sprintf('%-8i',EDF.DigMax)',8,EDF.NS)';
        h2(:,idx1(8)+1:idx1(9))=EDF.PreFilt;
        h2(:,idx1(9)+1:idx1(10))=reshape(sprintf('%-8i',EDF.SPR)',8,EDF.NS)';
        h2(h2==0)=32;
        for k=1:length(H2idx);
                fwrite(fid,h2(:,idx1(k)+1:idx1(k+1))','uchar');
        end;
else
        fwrite(fid, EDF.Label','16*uchar');
        fwrite(fid, EDF.Transducer','80*uchar');
        fwrite(fid, EDF.PhysDim','8*uchar');
        fwrite(fid, EDF.PhysMin,'float64');
        fwrite(fid, EDF.PhysMax,'float64');
        fwrite(fid, EDF.DigMin,'int64');
        fwrite(fid, EDF.DigMax,'int64');
        fwrite(fid, EDF.PreFilt','80*uchar');
        fwrite(fid, EDF.SPR,'uint32');
        fwrite(fid, EDF.GDFTYP,'uint32');
        fprintf(fid,'%c',32*ones(32,EDF.NS));
end;

tmp = ftell(EDF.FILE.FID);
if tmp ~= (256 * (EDF.NS+1)) 
        fprintf(1,'Warning %s-WRITE: incorrect header length %i bytes\n',upper(EDF.AS.Method),tmp);
%else   fprintf(1,'sdfopen in write mode: header info stored correctly\n');
end;        

EDF.AS.spb = sum(EDF.SPR);	% Samples per Block
11,
EDF.AS.bi = [0;cumsum(EDF.SPR)]; 
EDF.AS.BPR  = ceil(EDF.SPR.*GDFTYP_BYTE(EDF.GDFTYP+1)'); 
EDF.AS.SAMECHANTYP = all(EDF.AS.BPR == (EDF.SPR.*GDFTYP_BYTE(EDF.GDFTYP+1)')) & all(EDF.GDFTYP(:)~=EDF.GDFTYP(1));
%EDF.AS.GDFbi= [0;cumsum(EDF.AS.BPR)];
EDF.AS.GDFbi = [0;cumsum(ceil(EDF.SPR.*GDFTYP_BYTE(EDF.GDFTYP+1)'))]; 
EDF.AS.bpb = sum(ceil(EDF.SPR.*GDFTYP_BYTE(EDF.GDFTYP+1)'));	% Bytes per Block
EDF.AS.startrec = 0;
EDF.AS.numrec = 0;
EDF.FILE.POS = 0;

else % if arg2 is not 'r' or 'w'
        fprintf(EDF.FILE.stderr,'Warning %s: Incorrect 2nd argument. Argument2 must be ''r'' or ''w''\n',upper(EDF.AS.Method));
end;        

if EDF.ErrNo>0
        fprintf(EDF.FILE.stderr,'ERROR %i SDFOPEN\n',EDF.ErrNo);
end;
