% sdfread() -  Loads selected seconds of an EDF File (European Data Format for Biosignals)
%
% Usage:
%   [S,EDF] = sdfread(EDF [,NoS [,StartPos]] )
%
% Inputs:
%   EDF      = open EDF file structure returned by sdfopen();  ?????
%   NoS      = Number of seconds to read, default = 1 (second)
%   StartPos = Starting position. If not provided, the data is read
%              continously from the EDF file. Value 'Inf' reads to file end.
%              No repositioning of the file pointer is performed.
%
% Example:
%   EDF=sdfopen(Filename,'r',CHAN,ReRefMx,TSR,OFCHK);
%   [S,EDF] = sdfread(EDF, NoS, StartPos) % read NoS seconds of data from StartPos
% 
%   [S,EDF] = sdfread(EDF, EDF.NRec*EDF.Dur); % Reads all the file data
%   [S,EDF] = sdfread(EDF, inf);              % as does this
%
% Note:
% 1) Ver > 0.75 requests NoS and StartPos in seconds. 
%    Previously (Ver <0.76) the units were Records. 
% 2) In this version, modified for eeglab(), sdfseek, sdftell, 
%    sdfeof and dt have been included in the code of the function.
%
% Author: (C) 1997-2002 by Alois Schloegl, 15 Jun 2002 #0.85, (Header reworked for 
%         eeglab format and subfunctions inserted in the code, Arnaud Delorme, 27 Dec 2002)
%
% See also: fread, SDFREAD, SDFWRITE, SDFCLOSE, SDFSEEK, SDFREWIND, SDFTELL, SDFEOF

%	Copyright (c) 1997-2002 by Alois Schloegl
%	a.schloegl@ieee.org	
%	Version 0.85
%	24. Jun. 2002

% .changelog
% 13022002  V0.84    minor bug fixed
% 24062002  V0.85    check of input arguments corrected

function [S,EDF,OFCHK] = sdfread(EDF,NoS,StartPos)
global SDF_TB_TOGGLE_CHECK

if isfield(EDF.AS,'TOGGLE');
        if EDF.AS.TOGGLE~=SDF_TB_TOGGLE_CHECK
                fprintf(2,'Warning SDFREAD: [s,EDF]=sdfread(EDF, ...) \nYou forgot to pass EDF in %i call(s) of SDFREAD\n',SDF_TB_TOGGLE_CHECK-EDF.AS.TOGGLE);
        end;
else
        EDF.AS.TOGGLE=0;
        SDF_TB_TOGGLE_CHECK=0;
end;
SDF_TB_TOGGLE_CHECK = SDF_TB_TOGGLE_CHECK+1;
EDF.AS.TOGGLE = EDF.AS.TOGGLE+1;

OptiMEM=1;             % if you get low on memory (so your CPU starts swapping), 
                       % set this to one
OptiSPEED=~OptiMEM;
MAX_BLOCK_NUMBER=16;	% max. # of blocks in case of OptiMEM

GDF=strcmp(EDF.VERSION(1:3),'GDF');

if nargin<1
        fprintf(2,'error EDFREAD: missing input arguments\n');
end; if nargin<2
        NoS=1;
end; if nargin<3
        Mode=0;
end; if nargin<4 
        NoR=EDF.NRec;
end; 

if isfield(EDF,'SIE')   % SIESTA; Layer 4
        Mode_RAW=EDF.SIE.RAW; %any(Mode==[[4:7 12:15] [4:7 12:15]+16]);
	%Mode_CAL=EDF.SIE.RR; %any(Mode==[[2 3 6 7 10 11 14 15] [2 3 6 7 10 11 14 15]+16]);
        Mode_REF=0; %any(Mode==[[8:15] [8:15]+16]);
	Mode_SEC=0; %any(Mode==[[1:2:15] [1:2:15]+16]);
        Mode_RS100=EDF.SIE.RS; %bitand(Mode,2^4)>0;
        
        EDF.AS.startrec=EDF.AS.startrec+EDF.AS.numrec;
        EDF.AS.numrec=0;
        
        % prepare input arguments
        if ~EDF.SIE.TimeUnits_Seconds;        % units of time in records      
                NoR=NoS;%chan;
                if ~(nargin<3) 
                        [EDF]=sdfseek(EDF,StartPos,'bof'); 
                end;
                if isinf(NoR), 
                        NoR=EDF.NRec-EDF.FILE.POS; 
                end;
        else         % units of time in seconds
                % Prepare input arguments, transform NoS and StartPos into units of samples 
                if (nargin>2) %nargin==3 or larger
                        %tmp = floor(StartPos/EDF.Dur);
                        if StartPos+NoS > EDF.Dur*EDF.NRec;
                                fprintf(2,'Warning EDFREAD: %s has only %i seconds\n',EDF.FileName, EDF.Dur*EDF.NRec);
                                StartPos=min(StartPos*EDF.AS.MAXSPR/EDF.Dur,EDF.AS.MAXSPR*EDF.NRec); % transformation of seconds to samples + Overflowcheck
                                NoS=EDF.AS.MAXSPR*EDF.NRec-StartPos;
                        else
                                StartPos = StartPos*EDF.AS.MAXSPR/EDF.Dur;% EDF.Block.number(4);%/EDF.AS.MAXSPR*EDF.Dur;
                                NoS=NoS*EDF.AS.MAXSPR/EDF.Dur;
                        end;
                        NoR = ceil((StartPos+NoS)/EDF.AS.MAXSPR)-floor(StartPos/EDF.AS.MAXSPR);
                        %EDF.AS.A_Block_number=0;
                elseif nargin==2 %nargin==2
                        StartPos = EDF.Block.number(4); % EDF.Block.number(4);%/EDF.AS.MAXSPR*EDF.Dur;
                        NoS = NoS * EDF.AS.MAXSPR/EDF.Dur; % 
                else  % nargin<2
                        %NoR = 1, %ceil(NoS/EDF.Dur);
                        NoS = EDF.AS.MAXSPR/EDF.Dur; %EDF.Dur; % default one second
                        StartPos = EDF.Block.number(4);%/EDF.AS.MAXSPR*EDF.Dur; % EDF.FILE.POS*EDF.Dur;
                        %StartPos = EDF.FILE.POS*EDF.Dur,
                end;
                if isinf(NoS), 
                        NoS=EDF.NRec*EDF.AS.MAXSPR-EDF.Block.number(4); 
                end;
                
                % Q? whether repositioning is required, otherwise read continously, 
                if floor(StartPos/EDF.AS.MAXSPR)~=EDF.FILE.POS; % if start not next block 
                        if floor(StartPos/EDF.AS.MAXSPR)~=floor(EDF.Block.number(1)/EDF.AS.MAXSPR); % if start not same as Block.data
                                [EDF]=sdfseek(EDF,floor(StartPos/EDF.AS.MAXSPR),'bof');
                        end;
                end;
                % Q whether last bufferblock is needed or not
                if (StartPos >= EDF.Block.number(2)) | (StartPos < EDF.Block.number(1)) | diff(EDF.Block.number(1:2))==0
                        EDF.Block.number=[0 0 0 0 ];
                        EDF.Block.data=[];
                        NoR=ceil((StartPos+NoS)/EDF.AS.MAXSPR)-floor(StartPos/EDF.AS.MAXSPR);
                else
                        NoR=ceil((StartPos+NoS-EDF.Block.number(2))/EDF.AS.MAXSPR);
                end;
                
        end;
        
        
        %clear chan;
	InChanSelect=EDF.SIE.InChanSelect;
        Mode_CHANSAME = ~all(EDF.SPR(InChanSelect)==EDF.SPR(InChanSelect(1)));
        clear Mode;
        
else % Layer 3
	%if chan==0        
                InChanSelect=1:EDF.NS;
        %else
        %        InChanSelect=chan;
	%end;

	if ischar(Mode)         %%% not OCTAVE-compatible
    		arg3=upper(Mode);
    		Mode = any(arg3=='R')*8 + any(arg3=='W')*4 + (any(arg3=='A') & ~any(arg3=='N'))*2 + any(arg3=='S'); 
	end;

        Mode_RAW=any(Mode==[[4:7 12:15] [4:7 12:15]+16]);
        EDF.SIE.RR=any(Mode==[[2 3 6 7 10 11 14 15] [2 3 6 7 10 11 14 15]+16]);
        Mode_REF=any(Mode==[[8:15] [8:15]+16]);
        Mode_SEC=any(Mode==[[1:2:15] [1:2:15]+16]);
        Mode_RS100=0; %bitand(Mode,2^4)>0;
        EDF.SIE.FILT=0;
        
	% Position file pointer and calculate number of Records
	if Mode_SEC
    		if ~all(~rem([NoR StartPos],EDF.Dur))
            		fprintf(2,'Warning EDFREAD: NoR and/or StartPos do not fit to blocklength of EDF File of %i s.\n',EDF.Dur);
	        end;
    		StartPos=StartPos/EDF.Dur;
	        NoR=NoR/EDF.Dur;
	end;
	if ~(nargin<5) 
    	        EDF=sdfseek(EDF,StartPos,'bof'); 
        else
                EDF.AS.startrec=EDF.AS.startrec+EDF.AS.numrec;
                EDF.AS.numrec=0;
	end;

	if Mode_REF
    		ReRefMx=chan;
                InChanSelect=find(any(ReRefMx,1))';
	else 
                InChanSelect=reshape(chan,1,prod(size(chan)));
                ReRefMx=eye(EDF.NS);
                ReRefMx=ReRefMx(:,chan);
	end;

        tmp = EDF.SIE.RR + Mode_REF*2;
        if tmp==0
                EDF.Calib = [zeros(1,EDF.NS);eye(EDF.NS)];
                EDF.Calib = EDF.Calib(:,chan);
                %[zeros(1,EDF.NS); eye(EDF.NS)];
        elseif tmp==1
                %EDF.Calib=EDF.Calib;
                EDF.Calib = EDF.Calib(:,chan);
                %S=[ones(size(S,1),1) S]*EDF.Calib([1 chan+1],chan);
        elseif tmp==2
                EDF.Calib=[zeros(1,length(chan)); ReRefMx];
                %S=S*ReRefMx(chan,:);
        elseif tmp==3
                EDF.Calib=(EDF.Calib*ReRefMx);
                %S=[ones(size(S,1),1) S]*(EDF.Calib([1 chan+1],chan)*ReRefMx(chan,:));
        end;
        
	if exist('OCTAVE_VERSION')~=1
        	InChanSelect=find(any(EDF.Calib(2:EDF.NS+1,:),2))';
        else
                InChanSelect=find(any(EDF.Calib(2:EDF.NS+1,:)'));
        end;
        
        Mode_CHANSAME = ~all(EDF.SPR(InChanSelect)==EDF.SPR(InChanSelect(1)));
	%if any(EDF.SPR(chan)~=EDF.SPR(chan(1))) fprintf(2,'Warning EDFREAD: channels do not have the same sampling rate\n');end;
        clear chan;
end;
        
bi=[0;cumsum(EDF.SPR)];
Records=NoR;
maxspr=EDF.AS.MAXSPR; %max(EDF.SPR(EDF.SIE.ChanSelect));
count=0;

if Mode_RAW;
        S=zeros(sum(EDF.SPR(InChanSelect)),Records);
        bi=[0;cumsum(EDF.SPR(InChanSelect))];
else
        S=zeros(maxspr*Records,length(InChanSelect));
end;        

if ~EDF.AS.SAMECHANTYP; % all(EDF.GDFTYP(:)~=EDF.GDFTYP(1))
        idx0=0;
        count=0;
        while (count<Records*EDF.AS.spb) & ~sdfeof(EDF), 
                %   for K=1:length(InChanSelect), K=InChanSelect(k);
                for k=1:EDF.NS,
                        %if GDF-format should be used 
                        %datatyp=dt(EDF.GDFTYP(k));
                        %if exist('OCTAVE_VERSION')
                        %        tmp1=(datatyp~='');
                        %else
                        %        tmp1=~isempty(datatyp);
                        %end;
                        %if tmp1
                        %        [tmp,cnt]=fread(EDF.FILE.FID,EDF.SPR(k),datatyp);
                        %else 
                        %        fprintf(2,'Error SDFREAD: Invalid SDF channel type in %s at channel %i',EDF.FileName,k);
                        %end;
                        
                        [tmp,cnt]=fread(EDF.FILE.FID,EDF.SPR(k),int16);
                        if any(InChanSelect==k), K=find(InChanSelect==k);
                                if Mode_RAW;
                                        S(bi(K)+1:bi(K+1),count+1)=tmp;
                                else 
                                        if ~Mode_CHANSAME
                                                tmp=reshape(tmp(:,ones(1,maxspr/EDF.SPR(k)))',maxspr,1);
                                        end;
                                        %disp([size(tmp) size(S) k K l maxspr])
                                        S(idx0+(1:maxspr),K)=tmp; %RS%
                                end;  
                        end;      
                        count=count+cnt;
                end; 
                idx0 = idx0 + maxspr;
        end;
        S=S(1:idx0,:);
        count=count/EDF.AS.spb;
        
        %AFIR%
        if EDF.SIE.AFIR
                %EDF.AFIR.xin=S(:,find(EDF.AFIR.channel1==EDF.SIE.InChanSelect)); 
                EDF.AFIR.xin=S(:,EDF.AFIR.channel1); 
        end;
        
        % Overflow Check for SDF
        if nargout>2, %
                OFCHK=S(:,InChanSelect);
                if EDF.SIE.TH
                        for k=1:length(InChanSelect),K=InChanSelect(k);
                                OFCHK(:,k)=(S(:,k) < EDF.SIE.THRESHOLD(K,1)) + (S(:,k)>= EDF.SIE.THRESHOLD(K,2));
                        end;
                else 
                        fprintf(2,'Error SDFREAD (#SDF): No Overflow Threshold defined\n');
                end;
                OFCHK=OFCHK*abs(EDF.Calib([InChanSelect+1],:)~=0); %RS%
        end;
        if EDF.SIE.TH
                for k=1:length(InChanSelect),K=InChanSelect(k);
                        tmp=find(((S(:,k) < EDF.SIE.THRESHOLD(K,1)) + (S(:,k)>= EDF.SIE.THRESHOLD(K,2)))>0);
        %                S(tmp,k)=NaN;
                        %S(find(((S(:,k) < EDF.SIE.THRESHOLD(K,1)) + (S(:,k)>= EDF.SIE.THRESHOLD(K,2)))>0),k)=NaN;
                end;
        end;
        
else %%%%%%%%% if all channels of same type
        
datatyp=dt(EDF.GDFTYP(1));        
if Mode_RAW;
        [S, count]=fread(EDF.FILE.FID,[EDF.AS.spb, Records],datatyp);
        
        %count = floor(count/EDF.AS.spb);
        count = (count/EDF.AS.spb);
        if count<1; fprintf(2,'Warning EDFREAD: only %3.1f blocks were read instead  of %3.1f\n',count,Records); end;

        % Overflow Check for EDF-RAW
        %if nargout>2, %
                OFCHK = S;
                if EDF.SIE.TH>1
                        for k=1:length(InChanSelect),K=InChanSelect(k);
		                [y1,EDF.Block.z1{k}]=filter([1 -1],1,tmp,EDF.Block.z1{k});
            			[y2,EDF.Block.z2{k}]=filter(ones(1,EDF.SPR(K)/EDF.Dur)/(EDF.SPR(K)/EDF.Dur),1,y1==0,EDF.Block.z2{k});
            			%y2=y2>0.079*(EDF.SPR(k)/EDF.Dur);
    				[y3,EDF.Block.z3{k}]=filter(ones(1,EDF.SPR(K)/EDF.Dur)/(EDF.SPR(K)/EDF.Dur),1,(tmp>=EDF.SIE.THRESHOLD(K,2)) | (tmp<=EDF.SIE.THRESHOLD(K,1)),EDF.Block.z3{k});
            			%y3=y3>0.094*(EDF.SPR(k)/EDF.Dur);
				OFCHK(bi(k)+1:bi(k+1),:) = y3>0.094 | y2>0.079;
			end;
                elseif EDF.SIE.TH
                        for k=1:length(InChanSelect),K=InChanSelect(k);
                                OFCHK(bi(k)+1:bi(k+1),:)=(S(bi(k)+1:bi(k+1),:) < EDF.SIE.THRESHOLD(K,1)) + (S(bi(k)+1:bi(k+1),:)>= EDF.SIE.THRESHOLD(K,2));
                        end;
                        %S(OFCHK>0)=NaN;
                elseif nargout>2, %
                        fprintf(2,'Error SDFREAD (#EDF-RAW): No Overflow Threshold defined\n');
                end;
        %end;
        
	if Mode_REF
	    fprintf(2,'Warning SDFREAD: ReReferenzing "R" is not possible in combination with RAW "W"\n');
	end;
else
        if all(EDF.SPR(InChanSelect)==EDF.AS.MAXSPR)
                if ~OptiMEM % but OptiSPEED
                        [s, count]=fread(EDF.FILE.FID,[EDF.AS.spb,Records],datatyp);
                        count = (count/EDF.AS.spb);
                        S = zeros(maxspr*count,length(InChanSelect)); %%%
                        for k=1:length(InChanSelect), K=InChanSelect(k);
                                S(:,k)=reshape(s(bi(K)+1:bi(K+1),:),maxspr*count,1); %RS%
                        end;

                else %OptiMEM % but ~OptiSPEED
                        S = zeros(maxspr*Records,length(InChanSelect));%%%
                        idx0=0; count=0;
                        for l=1:ceil(Records/MAX_BLOCK_NUMBER),
                                tmp_norr=min(MAX_BLOCK_NUMBER,Records-count); % 16 = # of blocks read 
                                [s, C]=fread(EDF.FILE.FID,[EDF.AS.spb,tmp_norr],datatyp);
                                C = C/EDF.AS.spb;
                                for k=1:length(InChanSelect), K=InChanSelect(k);
                                        S(idx0+(1:maxspr*C),k) = reshape(s(bi(K)+1:bi(K+1),:),maxspr*C,1); %RS%
                                end;
                                idx0 = idx0 + maxspr*C;
                                count=count+C;
                        end;
                        S=S(1:idx0,:);
                end;
                %AFIR%
                if EDF.SIE.AFIR
                        EDF.AFIR.xin=S(:,EDF.AFIR.channel1); 
                end;
                
                % Overflow Check for EDF-all same sampling rate
                if nargout>2, %
                        OFCHK = zeros(size(S,1),length(InChanSelect));
                        if EDF.SIE.TH>1
                    		for k=1:length(InChanSelect),K=InChanSelect(k);
		            		[y1,EDF.Block.z1{k}]=filter([1 -1],1,tmp,EDF.Block.z1{k});
            				[y2,EDF.Block.z2{k}]=filter(ones(1,EDF.SPR(K)/EDF.Dur)/(EDF.SPR(K)/EDF.Dur),1,y1==0,EDF.Block.z2{k});
            				%y2=y2>0.079*(EDF.SPR(k)/EDF.Dur);
            				[y3,EDF.Block.z3{k}]=filter(ones(1,EDF.SPR(K)/EDF.Dur)/(EDF.SPR(K)/EDF.Dur),1,(tmp>=EDF.SIE.THRESHOLD(K,2)) | (tmp<=EDF.SIE.THRESHOLD(K,1)),EDF.Block.z3{k});
            				%y3=y3>0.094*(EDF.SPR(k)/EDF.Dur);
					OFCHK(bi(k)+1:bi(k+1),:) = y3>0.094 | y2>0.079;
				end;
            		elseif EDF.SIE.TH,
                                for k=1:length(InChanSelect),K=InChanSelect(k);
                                        OFCHK(:,k)=(S(:,k) <= EDF.SIE.THRESHOLD(K,1)) + (S(:,K)>= EDF.SIE.THRESHOLD(K,2));
                                end;
                                OFCHK = OFCHK*abs(EDF.Calib([InChanSelect+1],:)~=0); %RS%
                        else 
                                fprintf(2,'Error SDFREAD (#EDF-SPR=): No Overflow Threshold defined\n');
                        end;
                end;
                
                if 0; EDF.SIE.TH
                        for k=1:length(InChanSelect),K=InChanSelect(k);
                                tmp = S(:,k);
                                tmp = find(((tmp < EDF.SIE.THRESHOLD(K,1)) + (tmp >= EDF.SIE.THRESHOLD(K,2)))>0);
                %                S(tmp,k) = NaN;
                        end;
                end;
        else
                if ~OptiMEM % but OptiSPEED
                        [s, count]=fread(EDF.FILE.FID,[EDF.AS.spb,Records],datatyp);
                        count = (count/EDF.AS.spb);
                        S = zeros(maxspr*count,length(InChanSelect));
                        for k=1:length(InChanSelect), K=InChanSelect(k);
                                tmp=reshape(s(bi(K)+1:bi(K+1),:),EDF.SPR(K)*count,1);
                                if EDF.SPR(K) == maxspr
                                        S(:,k)=tmp;%reshape(tmp(:,ones(1,maxspr/EDF.SPR(K)))',maxspr*Records,1); %RS%
                                elseif EDF.SPR(K) > maxspr
                                        if rem(EDF.SPR(K)/maxspr,1)==0
                                                S(:,k)=rs(tmp,EDF.SPR(K),maxspr);
                                        else
                                                S(:,k)=rs(tmp,EDF.SIE.T); %(:,ones(1,maxspr/EDF.SPR(K)))',maxspr*Records,1); %RS%
                                        end;
                                elseif EDF.SPR(K) < maxspr
                                        S(:,k)=reshape(tmp(:,ones(1,maxspr/EDF.SPR(K)))',maxspr*count,1); %RS%
                                end;
                        end;    
                else %OptiMEM % but ~OptiSPEED
                        S = zeros(maxspr*Records,length(InChanSelect));
                        idx0=0; count=0;
                        for l=1:ceil(Records/MAX_BLOCK_NUMBER),
                        %while count < Records
                                tmp_norr = min(MAX_BLOCK_NUMBER,Records-count); % 16 = # of blocks read 
                                [s, C]=fread(EDF.FILE.FID,[EDF.AS.spb,tmp_norr],datatyp);
                                C = C/EDF.AS.spb;
                                for k=1:length(InChanSelect), K=InChanSelect(k); 
                                        tmp0=reshape(s(bi(K)+1:bi(K+1),:),EDF.SPR(K)*C,1);
                                        if EDF.SPR(K)==maxspr
                                                S(idx0+(1:maxspr*C),k)=tmp0;%reshape(tmp(:,ones(1,maxspr/EDF.SPR(K)))',maxspr*Records,1); %RS%
                                        elseif EDF.SPR(K)>maxspr
                                                if rem(EDF.SPR(K)/maxspr,1)==0
                                                        S(idx0+(1:maxspr*C),k)=rs(tmp0,EDF.SPR(K),maxspr);
                                                else
                                                        S(idx0+(1:maxspr*C),k)=rs(tmp0,EDF.SIE.T); %(:,ones(1,maxspr/EDF.SPR(K)))',maxspr*Records,1); %RS%
                                                end;
                                        elseif EDF.SPR(K)<maxspr
                                                S(idx0+(1:maxspr*C),k)=reshape(tmp0(:,ones(1,maxspr/EDF.SPR(K)))',maxspr*C,1); %RS%
                                        end;
                                end;
                                idx0 = idx0 + maxspr*C;
                                count= count+C;
                        end;
                        S=S(1:idx0,:);
                end;
                
                %AFIR%
                if EDF.SIE.AFIR
                        %EDF.AFIR.xin=S(:,find(EDF.AFIR.channel1==InChanSelect)); 
                        EDF.AFIR.xin=S(:,EDF.AFIR.channel1); 
                end;
                
                % Overflow Check for EDF with different sampling rates
                if nargout>2, %
                        if Mode_RS100
                                fprintf(2,'Warning SDFREAD: Overflow check for EDF files when reading channels with different sampling rate and resampling, is not implemented yet\n');
                        end;
                        OFCHK = S(:,InChanSelect);
                        if EDF.SIE.TH
                                for k=1:length(InChanSelect), K=InChanSelect(k);
                                        OFCHK(:,k)=(S(:,k) < EDF.SIE.THRESHOLD(K,1)) + (S(:,K)>= EDF.SIE.THRESHOLD(K,2));
                                end;
                        else 
                                fprintf(2,'Error SDFREAD (#EDF-SPR~): No Overflow Threshold defined\n');
                        end;
                end;
                if EDF.SIE.TH
                        for k=1:length(InChanSelect),K=InChanSelect(k);
                                tmp = S(:,k);
                                tmp = find(((tmp < EDF.SIE.THRESHOLD(K,1)) + (tmp >= EDF.SIE.THRESHOLD(K,2)))>0);
                %                S(tmp,k) = NaN;
                        end;
                end;
        end;
end;
end; % SDF

EDF.AS.numrec=count;
EDF.FILE.POS = EDF.AS.startrec + EDF.AS.numrec;
if EDF.AS.numrec~=Records, 
        fprintf(2,'Warning %s: %s only %i blocks instead of %i read\n',mfilename,EDF.FILE.Name,EDF.AS.numrec,Records);
end;

%%%%% Calibration of the signal 
if EDF.SIE.RR          % Autocalib 
        if Mode_RAW 
                for k=1:EDF.NS,
                        S(bi(k)+1:bi(k+1),:)=S(bi(k)+1:bi(k+1),:)*EDF.Cal(k)+EDF.Off(k);
                end;
        else
               %S=[ones(size(S,1),1) S(:,InChanSelect)]*EDF.Calib([1 InChanSelect+1],:);
                S=[ones(size(S,1),1) S]*EDF.Calib([1 InChanSelect+1],:); % EDF.Calib must be sparse, otherwise overflow-check is incorrect. 
        end;
end;

%%%%% Removing ECG Templates
if EDF.SIE.TECG          
        if (EDF.AS.startrec+EDF.AS.numrec)~=(ftell(EDF.FILE.FID)-EDF.HeadLen)/EDF.AS.bpb;
                fprintf(2,'ERROR SDFREAD: Mode TECG requires update of EDF [S,EDF]=sdfread(EDF,...)\n');
                EDF=sdftell(EDF);
        end;
        pulse = zeros(EDF.AS.numrec*EDF.SPR(12),1);
        
        Index=[];
        while EDF.TECG.idx(EDF.TECG.idxidx) <= EDF.FILE.POS*EDF.SPR(12)
                Index=[Index EDF.TECG.idx(EDF.TECG.idxidx)-EDF.AS.startrec*EDF.SPR(12)];
                EDF.TECG.idxidx=EDF.TECG.idxidx+1;
        end;

        if ~isempty(Index)
                pulse(Index) = 1;
        end;
        
        %tmp=find(EDF.TECG.idx > EDF.AS.startrec*EDF.SPR(12) & EDF.TECG.idx <= EDF.FILE.POS*EDF.SPR(12));
        %if ~isempty(tmp)
        %        pulse(EDF.TECG.idx(tmp)-EDF.AS.startrec*EDF.AS.MAXSPR) = 1;
        %end;
        
        for i=1:size(S,2),
                [tmp,EDF.TECG.Z(:,i)] = filter(EDF.TECG.QRStemp(:,i),1,pulse,EDF.TECG.Z(:,i));
                S(:,i)=S(:,i)-tmp; % corrects the signal
        end;
end;

%%%%% Filtering
if EDF.SIE.FILT
        for k=1:size(S,2);
                [S(:,k),EDF.Filter.Z(:,k)]=filter(EDF.Filter.B,EDF.Filter.A,S(:,k),EDF.Filter.Z(:,k));
        end;
        % Resampling of Overflow Check 
        if nargout>2, %
                if EDF.SIE.TH,
                        for k=1:size(S,2);
                                [S(:,k),EDF.FilterOVG.Z(:,k)]=filter(EDF.Filter.B,EDF.Filter.A,S(:,k),EDF.FilterOVG.Z(:,k));
                        end;
                end;
        end;
end;

%AFIR%
% Implements Adaptive FIR filtering for ECG removal in EDF/SDF-tb.
% based on the Algorithm of Mikko Koivuluoma <k7320@cs.tut.fi>
if EDF.SIE.AFIR
        e=zeros(size(S,1),length(EDF.AFIR.channel2));
        
        xinre2 = [EDF.AFIR.x EDF.AFIR.xin'];
        
        ddin=[EDF.AFIR.d; S]; %(:,EDF.AFIR.channel2)];
        
        for n=1:size(EDF.AFIR.xin,1),
                EDF.AFIR.x = xinre2(n + (EDF.AFIR.nord:-1:1)); % x(1:EDF.AFIR.nord-1)];
                
                y = EDF.AFIR.w * EDF.AFIR.x';
                
                en = EDF.AFIR.x * EDF.AFIR.x' + EDF.AFIR.gamma;
                e(n,:) = ddin(n,:) - y';
                
                EDF.AFIR.w = EDF.AFIR.w + (EDF.AFIR.alfa/en) * e(n,:)' * EDF.AFIR.x;
        end;
        EDF.AFIR.d = ddin(size(ddin,1)+(1-EDF.AFIR.delay:0),:);

        S=e; %output
	%S(:,EDF.AFIR.channel2) = e; %OUTPUT
end;

%%%%% select the correct seconds
if ~EDF.SIE.RAW & EDF.SIE.TimeUnits_Seconds 
        if NoR>0
                       %EDF.Block,
		       %[StartPos,StartPos+NoS,EDF.AS.startrec*EDF.Dur]
                if (StartPos < EDF.AS.startrec*EDF.AS.MAXSPR)
                        tmp = S(size(S,1)+(1-EDF.AS.MAXSPR:0),:);
                        if (nargout>2) & EDF.SIE.TH
                                tmp0 = OFCHK(size(S,1)+(1-EDF.AS.MAXSPR:0),:);
                        end;
                        
                        EDF.Block.number(3) = StartPos;
                        EDF.Block.number(4) = (StartPos+NoS);
                        
                        S = [EDF.Block.data(floor(EDF.Block.number(3)-EDF.Block.number(1))+1:EDF.AS.MAXSPR,:); S(1:floor(EDF.Block.number(4)-EDF.AS.startrec*EDF.AS.MAXSPR),:)];
                        if (nargout>2) & EDF.SIE.TH
                                OFCHK = [EDF.Block.data(floor(EDF.Block.number(3)-EDF.Block.number(1))+1:EDF.AS.MAXSPR,:); S(1:floor(EDF.Block.number(4)-EDF.AS.startrec*EDF.AS.MAXSPR),:)];
                                EDF.Block.dataOFCHK = tmp0;
                        end;
                        
                        EDF.Block.number(1:2)=(EDF.FILE.POS+[-1 0])*EDF.AS.MAXSPR;
                        EDF.Block.data = tmp;
                else
                        EDF.Block.number(3) = StartPos;
                        EDF.Block.number(4) = (StartPos+NoS);
                        EDF.Block.number(1:2)=(EDF.FILE.POS+[-1 0])*EDF.AS.MAXSPR;
                        EDF.Block.data = S(size(S,1)+(1-EDF.AS.MAXSPR:0),:);
                        
                        %[floor(EDF.Block.number(3)-EDF.AS.startrec*EDF.AS.MAXSPR)+ 1,floor(EDF.Block.number(4)-EDF.AS.startrec*EDF.AS.MAXSPR)],
                        S = S(floor(EDF.Block.number(3)-EDF.AS.startrec*EDF.AS.MAXSPR)+ 1:floor(EDF.Block.number(4)-EDF.AS.startrec*EDF.AS.MAXSPR),:);
                        if (nargout>2) & EDF.SIE.TH
                                EDF.Block.dataOFCHK = OFCHK(size(S,1)+(1-EDF.AS.MAXSPR:0),:);
                                OFCHK = OFCHK(floor(EDF.Block.number(3)-EDF.AS.startrec*EDF.AS.MAXSPR)+ 1:floor(EDF.Block.number(4)-EDF.AS.startrec*EDF.AS.MAXSPR),:);
                        end;
                        %S = S(EDF.AS.MAXSPR/EDF.Dur*(rem(StartPos,EDF.Dur))+(1:NoR*EDF.AS.MAXSPR),:);
                end;
                
        else
                EDF.Block.number(3) = StartPos;
                EDF.Block.number(4) = (StartPos+NoS);
                EDF.Block.number(1:2)=(EDF.FILE.POS+[-1 0])*EDF.AS.MAXSPR;
                S = [EDF.Block.data(floor(EDF.Block.number(3)-EDF.Block.number(1))+1:floor(EDF.Block.number(4)-EDF.Block.number(1)),:)];
                if (nargout>2) & EDF.SIE.TH
                        OFCHK = [EDF.Block.dataOFCHK(floor(EDF.Block.number(3)-EDF.Block.number(1))+1:floor(EDF.Block.number(4)-EDF.Block.number(1)),:)];
                end;
                %S = EDF.Block.data(floor(StartPos*EDF.AS.MAXSPR/EDF.Dur+1:floor((StartPos+NoS)*EDF.AS.MAXSPR/EDF.Dur),:)];
        end;
end;
%end;

%%%%% Resampling
if Mode_RS100 & ~Mode_RAW
        S=rs(S,EDF.SIE.T); %RS%
        % Resampling of Overflow Check 
        if nargout>2, %
                if EDF.SIE.TH,
                        OFCHK=rs(OFCHK, EDF.SIE.T); %RS%
                end;
        end;
end;


function [EDF]=sdftell(EDF)
% EDF=sdftell(EDF_Struct)
% returns the location of the EDF_file position indicator in the specified file.  
% Position is indicated in Blocks from the beginning of the file.  If -1 is returned, 
% it indicates that the query was unsuccessful; 
% EDF_Struct is a struct obtained by sdfopen.
%
% EDF.FILE.POS contains the position of the EDF-Identifier in Blocks


%	Version 0.85
%	15 Jun 2002
%	Copyright (c) 1997-2002 by Alois Schloegl
%	a.schloegl@ieee.org	


POS = ftell(EDF.FILE.FID);
if POS<0
        [EDF.ERROR,EDF.ErrNo] = ferror(EDF.FILE.FID);
        return; 
end;
EDF.FILE.POS = (POS-EDF.HeadLen)/EDF.AS.bpb;
EDF.ERROR=[];
EDF.ErrNo=0;

if (EDF.AS.startrec+EDF.AS.numrec)~=EDF.FILE.POS
        fprintf(2,'Warning SDFTELL: File postion error in EDF/GDF/SDF-toolbox.\n')
        EDF.AS.startrec = EDF.FILE.POS;
end;        


function [status]=sdfeof(EDF)
% sdfeof(EDF)
% returns 1 if End-of-EDF-File is reached
% returns 0 otherwise

%	Copyright (c) 1997-99 by Alois Schloegl
%	a.schloegl@ieee.org	
%	Version 0.60
%	16. Aug. 1999


% status=feof(EDF.FILE.FID);  % does not work properly
%if EDF.FILE.POS~=EDF.AS.startrec+EDF.AS.numrec;
        
status=(EDF.FILE.POS>=EDF.NRec);



function [EDF]=sdfseek(EDF,offset,origin)
% [EDF]=sdfseek(EDF,offset,origin)
% Currently, offset and origin are the number of (EDF) records. 
% EDF.status contains the value of the older version of [status]=gdfseek(...) 
%
% See also: FSEEK, SDFREAD, SDFWRITE, SDFCLOSE, SDFREWIND, GSFTELL, SDFEOF

%	Copyright (c) 1997-2002 by Alois Schloegl
%	a.schloegl@ieee.org	
%	Version 0.85
%	15. Juni 2002

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


if strcmp(origin,'bof')
	origin=-1;        
elseif strcmp(origin,'cof')
	origin=0;        
elseif strcmp(origin,'eof')
	origin=1;        
end;

if origin==-1 
	EDF.FILE.POS = offset;
        OFFSET = EDF.AS.bpb*offset;
        status = fseek(EDF.FILE.FID,EDF.HeadLen+OFFSET,-1);
elseif origin==0 
	EDF.FILE.POS = EDF.FILE.POS + offset;
        OFFSET = EDF.AS.bpb*offset;
        status = fseek(EDF.FILE.FID,OFFSET,0);
elseif origin==1 
	EDF.FILE.POS = EDF.NRec+offset;
        OFFSET = EDF.AS.bpb*offset;
        status = fseek(EDF.FILE.FID,OFFSET,1);
else
        fprintf(2,'error SDFSEEK: 3rd argument "%s" invalid\n',origin);
        return;
end;


EDF.AS.startrec=EDF.FILE.POS;
EDF.AS.numrec = 0;
EDF = sdftell(EDF); % not really needed, only for double check of algorithms

% Initialization of Bufferblock for random access (without EDF-blocklimits) of data 
if ~EDF.SIE.RAW & EDF.SIE.TimeUnits_Seconds
        EDF.Block.number=[0 0 0 0]; %Actual Blocknumber, start and end time of loaded block, diff(EDF.Block.number(1:2))==0 denotes no block is loaded;
        % EDF.Blcok.number(3:4) indicate start and end of the returned data, [units]=samples.
        EDF.Block.data=[];
        EDF.Block.dataOFCHK=[];
end;


if 1; %isfield(EDF,'AFIR');
        if EDF.SIE.AFIR
                EDF.AFIR.w = zeros(EDF.AFIR.nC,max(EDF.AFIR.nord));
                EDF.AFIR.x = zeros(1,EDF.AFIR.nord);
                EDF.AFIR.d = zeros(EDF.AFIR.delay,EDF.AFIR.nC);
                fprintf(2,'WARNING SDFSEEK: Repositioning deletes AFIR-filter status\n');
        end;
end;
if 1; %isfield(EDF,'Filter');
        if EDF.SIE.FILT
                [tmp,EDF.Filter.Z]=filter(EDF.Filter.B,EDF.Filter.A,zeros(length(EDF.Filter.B+1),length(EDF.SIE.ChanSelect)));
                EDF.FilterOVG.Z=EDF.Filter.Z;
                fprintf(2,'WARNING SDFSEEK: Repositioning deletes Filter status of Notch\n');
        end;
end;

if 1; %isfield(EDF,'TECG')
        if EDF.SIE.TECG
                
                fprintf(2,'WARNING SDFSEEK: Repositioning deletes TECG filter status\n');
        end;
end;
EDF.FILE.status=status;

function datatyp=dt(x)
k=1;
EDF.GDFTYP(1)=x;
if EDF.GDFTYP(k)==0
        datatyp=('uchar');
elseif EDF.GDFTYP(k)==1
        datatyp=('int8');
elseif EDF.GDFTYP(k)==2
        datatyp=('uint8');
elseif EDF.GDFTYP(k)==3
        datatyp=('int16');
elseif EDF.GDFTYP(k)==4
        datatyp=('uint16');
elseif EDF.GDFTYP(k)==5
        datatyp=('int32');
elseif EDF.GDFTYP(k)==6
        datatyp=('uint32');
elseif EDF.GDFTYP(k)==7
        datatyp=('int64');
elseif 0; EDF.GDFTYP(k)==8
        datatyp=('uint64');
elseif EDF.GDFTYP(k)==16
        datatyp=('float32');
elseif EDF.GDFTYP(k)==17
        datatyp=('float64');
elseif 0;EDF.GDFTYP(k)>255 & EDF.GDFTYP(k)< 256+64
        datatyp=(['bit' int2str(EDF.GDFTYP(k))]);
elseif 0;EDF.GDFTYP(k)>511 & EDF.GDFTYP(k)< 511+64
        datatyp=(['ubit' int2str(EDF.GDFTYP(k))]);
elseif EDF.GDFTYP(k)==256
        datatyp=('bit1');
elseif EDF.GDFTYP(k)==512
        datatyp=('ubit1');
elseif EDF.GDFTYP(k)==255+12
        datatyp=('bit12');
elseif EDF.GDFTYP(k)==511+12
        datatyp=('ubit12');
elseif EDF.GDFTYP(k)==255+22
        datatyp=('bit22');
elseif EDF.GDFTYP(k)==511+22
        datatyp=('ubit22');
elseif EDF.GDFTYP(k)==255+24
        datatyp=('bit24');
elseif EDF.GDFTYP(k)==511+24
        datatyp=('ubit24');
else 
        fprintf(2,'Error GDFREAD: Invalid GDF channel type\n');
        datatyp='';
end;
