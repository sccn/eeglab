% fmrib_fastr() -  Remove FMRI gradient artifacts from EEG data using
%   FMRI artifact Slice Template Removal (FASTR) [Niazy06] and uses optimal
%   basis sets (OBS) to reduce the residuals.  
%
%   This subtraction algorithm is based on the principles outlined by 
%   [Niazy06] with few improvements and modifications. This program 
%   constructs a unique artifact template for each slice then subtracts the 
%   artifact.  Residual artifacts are removed by constructing a matrix of 
%   the residuals, doing a PCA then fitting the first 4 PCs (Optimal basis 
%   set - OBS)to the residuals in each slice.  This procedure should not be 
%   applied for non-EEG channels (e.g. ECG) as it can remove some high frequency 
%   details from a signal (e.g. QRS complex). Adaptive noise cancellation 
%   (ANC) [Allen00] is then used.  
%
%
%   USAGE:
%   EEG=fmrib_fastr(EEG,lpf,L,window,Trigs,strig,anc_chk,tc_chk,Volumes,Slices,varargin);
%   EEG:  EEGLAB data structure
%   lpf:  low-pass filter cutoff
%   L: Interpolation folds
%   window: length of averaging window in number of artifacts
%   Trigs: An array of slice triggers locations.
%   strig: 1 for slice triggers, 0 for volume / section triggers.
%   anc_chk: 1 to do Adaptive noise cancellation
%			 0 to not.
%   tc_chk:  1 to correct for missing triggers, 0 for not
%   Volumes: FMRI volumes for use in trigger correction
%   Slices:  FMRI Slices / Volume for use in trigger correction
%   varargin{1}: relative position of slice trigger from beginning of
%       slice acquisition: 0 for exact start -> 1 for exact end
%       default=0.03;
%   varargin{2}: Channels not to perform OBS  on.
%   varargin{3}: Numer of PCs to use in OBS. use 0 to skip this step.
%                'auto' or empty for auto order selection.
%
% [Niazy06] R.K. Niazy, C.F. Beckmann, G.D. Iannetti, J.M. Brady, and 
%  S.M. Smith (2005) Removal of FMRI environment artifacts from EEG data 
%  using optimal basis sets. NeuroImage 28 (3), pages 720-737.
%
%
% [Allen00] A Method for Removing Imaging Artifact from Continuous EEG 
%   Recording during Functional MRI, P.J. Allen, O. Josephs, R. Turner.
%   NeuroImage 12, 230-239 (2000).
%
%
% 
%
%   Author:  Rami K. Niazy, FMRIB Centre, Univ. of Oxford.
%   
%   Copyright (c) 2006 University of Oxford

% Copyright (C) 2006 University of Oxford
% Author:   Rami K. Niazy, FMRIB Centre
%           rami@fmrib.ox.ac.uk
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


% 31 MAR 2006
% PCA of residual now done on HPF data at 70Hz.  This resolved overfitting
% issues, which sometimes caused removal of data.

% 10 MAR 2006
% Fixed potentially sereous bug related to command line input order of varargin

% 13 JAN 2006
% ECG Channel Artifact now not least-square fitted

% 13 JAN 2006
% Fixed Section Marker Limit for s=sections case

% 16 SEP 2005
% Fixed bug when processing non-eeg channels

% 15 SEP 2005
% Added code to lengthen the section to insure enough triggers are
% contained.

% 13 SEP 2005
% Edited use of 'find' to be backward compatible with MATLAB 6

% 02 SEP 2005
% Fixed typo bug when interpolation = 1;
% Fixed problem when using '/' with single precision numbers.

% 01 SEP 2005
% Commented out test code causing problems when not in 'auto' mode.

% 09 AUG 2005
% Automatic order selection of residual PCs added.

% 05 AUG 2005
% search window for max correlation now determined by original fs.
% Fixed Misc bugs for volume/section triggers
% ANC made optional.

% 12 MAR 2005
% Program made to work for interleaved recording / Volume triggers
% Fixed various bugs

% 23 JAN 2005
% Fixed rounding of Filter Order

% 14 JAN 2005
% Updated Help

% 23 DEC 2004
% Added use of binary
% adaptive noise cancellation
% binary corr (pcorr2)
% Updated (c) info

% 21 DEC 2004
% Fixed bug for
% sections==1
% made ANC called from file

% 17 DEC 2004
% Use CORRCOEF instead of
% CORR2 to eliminate need for
% image porc toolbox

% 14 DEC 2004
% Major Updates:
%   Added RAPCO.
%   Input for trig events replaced with an array of Trigs instead 
%   to simplify scripting.
%   Other misc bugs fixed.

% 06 OCT 2004
% Edited for possibility of L=1
% i.e. no interpolation / decimation



function EEG=fmrib_fastr(EEG,lpf,L,Window,Trigs,strig,anc_chk,tc_chk,Volumes,Slices,varargin)

EEG.data = double(EEG.data);

tic;
% Check Input & Intialize Params
% ------------------------------
if nargin < 8
    error('Incorrect number of input arguments.','fmrib_fastr() error!');
end

if tc_chk==1
    if nargin < 10
        error('Volumes and Slices must be included in input arguments - See help.','fmrib_fastr() error!');
    end
end

if nargin > 13
    error('Incorrect number of input arguments.','fmrib_fastr() error!');
end

pre_frac=0.03;
npc='auto';
exc=[];
hpf=70;
fs=EEG.srate;
SecT=60;
if isempty(lpf)
    lpf=0;
end

if anc_chk==1 & lpf==0 
    lpf=70;
end

if rem(Window,2)~=0
    Window=Window+1;
end
if strig==1
    hwinl=Window;
else
    hwinl=Window/2;
end

scount=0;
[m n]=size(EEG.data);
minorder=15;
minfac=3;
STARTFLAG=0;
LASTFLAG=0;
TH_SLOPE=2;
TH_CUMVAR=80;
TH_VAREXP=5;

if nargin>10
    pre_frac=varargin{1};
end
if nargin>11
    exc=varargin{2};
    if ~isempty(exc)
        exc=sort(exc);
        if exc(end)>m | exc(1)<1
            error('Channel to exclude from OBS out of range',...
                'fmrib_fastr() error!');
        end
    end
end
if nargin>12
    if ischar(npc)
        if ~strcmp(lower(npc),'auto')
            error('unknown PC option: %s\n',npc);
        end
    end
    npc=varargin{3};
end


% check for DSP Toolbox
%----------------------
if ~exist('firls')
   error('FASTR requires the signal processing toolbox.','fmrib_fastr() error!');
end

% Extract triggers from Event field and correct triggers if needed
% ----------------------------------------------------------------

markerch(Trigs)=-1;
markerch(n)=0;
if tc_chk==1
    if length(Trigs)~=(Slices*Volumes)
        fprintf('correcting slice triggers...\n')
        markerch=trigcorrect(markerch,Slices,Volumes,1);
        Trigs=[];
        Trigs=find(markerch==-1);
    end
end
ml=length(Trigs);

if strig==1
    if Window>floor((ml-3)/2)
        Window=floor((ml-3)/2);
        if rem(Window,2)~=0
            Window=Window-1;
        end
        hwinl=Window;
        warning('For %d slice triggers, the maximum averaging\nwindow length allowed is %d. Window length set to %d.\n',...
            ml,  Window,Window);
    end
elseif strig==0
    if Window>(ml-2)       
        Window=ml-2;
        if rem(Window,2)~=0
            Window=Window-1;
        end
        hwinl=Window/2;
        warning('For %d volume triggers, the maximum averaging\nwindow length allowed is %d. Window length set to %d.\n',...
            ml,  Window,Window);
    end
end
       
        

% Calculate artifact dimentions
%-------------------------------   
clear marker markerch peaks ;

marker=Trigs*L;
clear dmarker;

dmarker=diff(marker);


min_isi=ceil(median(dmarker));
max_isi=ceil(1.01*min_isi);
end_artlength=max_isi-min_isi;
pre_peak=round(min_isi*pre_frac+end_artlength);
post_peak=round((1-pre_frac)*min_isi);
max_postpeak=round(post_peak+end_artlength);
art_length=pre_peak+post_peak+1;
max_artlength=art_length+end_artlength;
searchw=round(3*L);
pad_sec=ceil((art_length/L)*(hwinl/2+1)/fs)*fs;


% calculating lpf filter weights
% ------------------------------
trans=0.15;
nyq=0.5*fs;

if lpf>0
    if (1+trans)*lpf/nyq > 1
        error('LPF cutoff frequency to close to Nyquist frequency.',...
            'fmrib_fastr() error!');
    end

    filtorder=round(minfac*fix(fs/lpf));

    if filtorder < minorder
        filtorder=minorder;
    end

    if rem(filtorder,2)~=0
        filtorder=filtorder+1;
    end

    f=[0 lpf/nyq lpf*(1+trans)/nyq 1];
    a=[1 1 0 0];
    lpfwts=firls(filtorder,f,a);
end
    

% calculating hpf filter weights for template formation
% -----------------------------------------------------

filtorder=round(1.2*fs*L/(hpf-10));

if rem(filtorder,2)~=0
    filtorder=filtorder+1;
end

f=[0 (hpf-10)/(nyq*L) (hpf+10)/(nyq*L) 1]; 
a=[0 0 1 1];
hpfwts=firls(filtorder,f,a);

% calculating hpf of ANC
% ----------------------
if strig==1
    Tr=1;
    while Tr<=length(Trigs)
        Trtime=Trigs(Tr+1)-Trigs(1);
        if Trtime>=fs
            break
        end
        Tr=Tr+1;
    end
    ANCf=0.75*Tr;
else
    ANCf=2;
end
filtorder=round(1.2*fs/(ANCf*(1-trans)));
if rem(filtorder,2)~=0
    filtorder=filtorder+1;
end
f=[0 ANCf*(1-trans)/nyq ANCf/nyq 1];
ANCfwts=firls(filtorder,f,a);

% Other & ANC settings
%----------------------
peaks=zeros(1,n);
peaks(Trigs)=1;
c=1;
tcount=0;
secmarker=[];
SCount=1;
sec=1;
d1=Trigs(1)-ceil(1.25*art_length/L);
if d1 <= 0
    d1=1;
end
d2=Trigs(end)+ceil(2.25*art_length/L);
if d2 > length(EEG.data)
    d2=length(EEG.data);
end
N=double(ceil(max_artlength/L));
mANC=d2-d1+1;
sections=floor(mANC/(floor(fs*SecT)))+1;
N=double(N);

% Validate section length has enough triggers and adjust
% ------------------------------------------------------
SecTFlag=0;
SecTENDFlag=0;
if sections>1
    while SecTFlag==0
        if (d1-1)+floor(SecT*fs)+pad_sec > EEG.pnts
            sections=1;
            SecTENDFlag=1;
            break;
        end
          
        secpeaks=peaks(d1:(d1-1)+floor(SecT*fs)+pad_sec);
        tmpmarker=find(secpeaks==1)*L;
        if strig==1
            if round(length(tmpmarker)/2)<1.5*Window
                SecT=SecT+0.1;
            else
                SecTFlag=1;
            end
        elseif strig==0
            if length(tmpmarker)<1.5*Window
                SecT=SecT+0.1;
            else
                SecTFlag=1;
            end
        end
    end
    if SecTENDFlag==0;
        sections=floor(mANC/(floor(fs*SecT)))+1;
    end
end

% Allocate Memory
% ---------------

fprintf('Allocating memory...\n');

if strig==1
    slice_art1=zeros(hwinl+1,art_length);
    slice_art2=zeros(hwinl+1,art_length);
else
    slice_art1=zeros(2*hwinl+1,art_length);
end
avg_art=zeros(1,art_length);
Idata=zeros(1,(floor(fs*SecT)+2*pad_sec)*L);
Iorig=zeros(1,(floor(fs*SecT)+2*pad_sec)*L);
INoise=zeros(1,(floor(fs*SecT)+2*pad_sec)*L);
fNoise=zeros(1,(floor(fs*SecT)+2*pad_sec));
Noise=zeros(1,n);
fcleanEEG=zeros(1,(floor(fs*SecT)+2*pad_sec));
cleanEEG=zeros(1,n);
tmpd=zeros(n,1);
d=zeros(mANC,1);
refs=zeros(mANC,1);
out=zeros(mANC,1);
y=zeros(mANC,1);
W=zeros(N+1,1);

% Align Slice triggers according to first channel
% -----------------------------------------------

while sec<=sections
    
    if sec==1 & sections > 1
        rempeaks=peaks((d1-1)+sec*floor(fs*SecT)+1:d2);
        remmarker=find(rempeaks==1);
        if length(remmarker) > 3*hwinl
            secpeaks=peaks(d1:(d1-1)+sec*floor(fs*SecT)+pad_sec);
        else
            secpeaks=peaks(d1:d2);
            sec=sections;
        end
    elseif sec==1 & sections==1
        secpeaks=peaks(d1:d2);
    elseif sec==sections
        secpeaks=peaks((d1-1)+(sec-1)*floor(fs*SecT)+1-pad_sec:d2);
    else
        rempeaks=peaks((d1-1)+sec*floor(fs*SecT)+1-pad_sec:d2);
        remmarker=find(rempeaks==1);
        if length(remmarker) > 3*hwinl
            secpeaks=peaks((d1-1)+(sec-1)*floor(fs*SecT)+1-pad_sec:...
                (d1-1)+sec*floor(fs*SecT)+pad_sec);
        else
            secpeaks=peaks((d1-1)+(sec-1)*floor(fs*SecT)+1-pad_sec:d2);
            sec=sections;
        end
    end
    
    tmpmarker=find(secpeaks==1)*L;
    markerl(SCount)=length(tmpmarker);
    if sec==1 | sec==sections
        secl(SCount)=length(secpeaks)-pad_sec;
    else
        secl(SCount)=length(secpeaks)-2*pad_sec;
    end
    secmarker=[secmarker tmpmarker];
    SCount=SCount+1;
    sec=sec+1;
end

% test prcorr2
try
    prcorr2(rand(1,100), rand(1,100))
catch,
    tmpfile = which('prcorr2');
    [tmppath tmpfilenoext] = fileparts(tmpfile);
    tmpfilenoext = fullfile(tmppath,tmpfilenoext);
    delete(tmpfile);
    disp( [ 'Removing file ' tmpfile ]);
    disp( [ 'Try recompiling the mex file by typing' ]);
    disp( [ 'mex ' tmpfilenoext '.c' ]);
end;

%----------------
SCount=1;
sections=length(secl);
steps=(sections+1)*m;
pcamat=zeros(floor(max(markerl)/2),pre_peak+max_postpeak+1);

for sec=1:sections
 
    if sec==1
        barth=5;
        barth_step=barth;
        Flag25=0;
        Flag50=0;
        Flag75=0;
        fprintf('\nStage 1 of 2: Slice Alignment\n')
    end
    
    if L > 1
        if sec==1 & sections > 1
            Idata=interp(EEG.data(c,d1:(d1-1)+secl(sec)+pad_sec),L,4,1);
        elseif sec==1 & sections==1
            Idata=interp(EEG.data(c,d1:d2),L,4,1);
        elseif sec==sections
            Idata=interp(EEG.data(c,(d1-1)+sum(secl(1:sec-1))+1-pad_sec:...
                d2),L,4,1);
        else
            Idata=interp(EEG.data(c,(d1-1)+sum(secl(1:sec-1))+1-pad_sec:...
                (d1-1)+sum(secl(1:sec))+pad_sec),L,4,1);
        end
    else
        if sec==1 & sections > 1
            Idata=EEG.data(c,d1:(d1-1)+secl(sec)+pad_sec);
        elseif sec==1 & sections==1
            Idata=EEG.data(c,d1:d2);
        elseif sec==sections
            Idata=EEG.data(c,(d1-1)+sum(secl(1:sec-1))+1-pad_sec:d2);
        else
            Idata=EEG.data(c,(d1-1)+sum(secl(1:sec-1))+1-pad_sec:...
                (d1-1)+sum(secl(1:sec))+pad_sec);
        end
    end
  
    if sec==1
        ml=0;
        if sections>1
            for nsec=2:sections
                if nsec==sections
                    starts=sum(markerl(1:nsec-1))+1;
                    lasts=sum(markerl(1:nsec));
                else
                    starts=sum(markerl(1:nsec-1))+2;
                    lasts=sum(markerl(1:nsec))-2;
                end
                ml=ml+(lasts-starts+1);
            end
        end
        starts=1;
        if sections>1
            lasts=markerl(SCount)-2;
        else
            lasts=markerl(SCount);
        end
        ml=ml+(lasts-starts+1);
    elseif sec==sections
        starts=sum(markerl(1:SCount-1))+1;
        lasts=sum(markerl(1:SCount));
    else
        starts=sum(markerl(1:SCount-1))+2;
        lasts=sum(markerl(1:SCount))-2;
    end
    
    
    for s=starts:lasts
        if s==starts & sec==1
            try
                slice_art(1,:)=...
                    Idata(secmarker(s)-pre_peak:secmarker(s)+post_peak);
            catch
                slice_art(1,:)=...
                    Idata(secmarker(s+1)-pre_peak:secmarker(s+1)+post_peak);
            end
        end
       
        if sec==sections & s==lasts
            try
                ppn=1;
			    for pp=secmarker(s)-searchw:secmarker(s)+searchw
                    match(ppn)=prcorr2(slice_art(1,:),...
                        Idata(pp-pre_peak:pp+post_peak));
                    ppn=ppn+1;
                end
				[CV,CP]=max(match);
				adjust=CP-(searchw+1);
				secmarker(s)=secmarker(s)+adjust;
            catch
            end
        elseif s>=starts+1; 
			ppn=1;
			for pp=secmarker(s)-searchw:secmarker(s)+searchw
                match(ppn)=prcorr2(slice_art(1,:),...
                    Idata(pp-pre_peak:pp+post_peak));
                ppn=ppn+1;
			end
			[CV,CP]=max(match);
			adjust=CP-(searchw+1);
			secmarker(s)=secmarker(s)+adjust;
        end
        
        tcount=tcount+1;
        percentdone=floor(tcount*100/ml);
        
        if floor(percentdone)>=barth
            
            if percentdone>=25 & Flag25==0
                fprintf('25%% ')
                Flag25=1;
            elseif percentdone>=50 & Flag50==0
                fprintf('50%% ')
                Flag50=1;
            elseif percentdone>=75 & Flag75==0
                fprintf('75%% ')
                Flag75=1;
            elseif percentdone==100
                fprintf('100%%\n')
            else
                fprintf('.')
            end
            
            while barth<=percentdone
                barth=barth+barth_step;
            end
            if barth>100
                barth=100;
            end
        end 
        
    end
    SCount=SCount+1;
end
secmarker2=secmarker;

% Construct Artifacts and Subtract
% ---------------------------------

for c=1:m
     
    % Progress bar Init
    if c==1
        barth=5;
        barth_step=barth;
        Flag25=0;
        Flag50=0;
        Flag75=0;
        fprintf('\nStage 2 of 2: Artifact Subtraction\n');
    end




    tmpdata=EEG.data(c,:)-mean(EEG.data(c,:));

    cleanEEG=EEG.data(c,:);

    % Process in sections of SecT seconds for memory concerns
    for sec=1:sections
        
        if L > 1
            if sec==1 & sections > 1
                Idata=interp(tmpdata(d1:(d1-1)+secl(sec)+pad_sec),L,4,1);
                Iorig=interp(EEG.data(c,d1:(d1-1)+secl(sec)+pad_sec),L,4,1);
            elseif sec==1 & sections==1
                Idata=interp(tmpdata(d1:d2),L,4,1);
                Iorig=interp(EEG.data(c,d1:d2),L,4,1);
            elseif sec==sections
                Idata=interp(tmpdata((d1-1)+...
                    sum(secl(1:sec-1))+1-pad_sec:d2),L,4,1);
                Iorig=interp(EEG.data(c,(d1-1)+...
                    sum(secl(1:sec-1))+1-pad_sec:d2),L,4,1);
            else
                Idata=interp(tmpdata((d1-1)+sum(secl(1:sec-1))+...
                    1-pad_sec:(d1-1)+sum(secl(1:sec))+pad_sec),L,4,1);
                Iorig=interp(EEG.data(c,(d1-1)+sum(secl(1:sec-1))+...
                    1-pad_sec:(d1-1)+sum(secl(1:sec))+pad_sec),L,4,1);
            end
        else
            if sec==1 & sections > 1
                Idata=tmpdata(d1:(d1-1)+secl(sec)+pad_sec);
                Iorig=EEG.data(c,d1:(d1-1)+secl(sec)+pad_sec);
            elseif sec==1 & sections==1
                Idata=tmpdata(d1:d2);
                Iorig=EEG.data(c,d1:d2);
            elseif sec==sections
                Idata=tmpdata((d1-1)+sum(secl(1:sec-1))+1-pad_sec:d2);
                Iorig=EEG.data(c,(d1-1)+sum(secl(1:sec-1))+1-pad_sec:d2);
            else
                Idata=tmpdata((d1-1)+sum(secl(1:sec-1))+1-pad_sec:...
                    (d1-1)+sum(secl(1:sec))+pad_sec);
                Iorig=EEG.data(c,(d1-1)+sum(secl(1:sec-1))+1-pad_sec:...
                    (d1-1)+sum(secl(1:sec))+pad_sec);
            end
        end
        
        INoise=zeros(1,length(Idata));
              
        %Average Artifacts & Subtract
        
        if sections==1
            starts=1;
            lasts=markerl(sec);
        elseif sec==1
            starts=1;
            lasts=markerl(sec)-2;
        elseif sec==sections
            starts=sum(markerl(1:sec-1))+2;
            lasts=sum(markerl(1:sec));
        else
            starts=sum(markerl(1:sec-1))+2;
            lasts=sum(markerl(1:sec))-2;
        end
        
        for s=starts:lasts
            
            if strig==1 % Slice Triggers
         
                if s==starts
                    art=1;
                    ssc=1;
                    for ss=starts+1:2:starts+2*hwinl+1
                        slice_art1(ssc,:)=...
                            Idata(secmarker(ss)-pre_peak:...
                            secmarker(ss)+post_peak);
                        ssc=ssc+1;
                    end
                    avg_art=mean(slice_art1,1); 
                elseif s==starts+1
                    ssc=1;
                    for ss=starts+2:2:starts+2*hwinl+2
                        slice_art2(ssc,:)=...
                            Idata(secmarker(ss)-pre_peak:...
                            secmarker(ss)+post_peak);
                        ssc=ssc+1;
                    end
                    avg_art=mean(slice_art2,1); 
                elseif ((s>(starts+hwinl+2)) & (s<=(lasts-(hwinl+2))))
                    ss=s+hwinl;
                    switch art
                        case 1
                            slice_art1=[slice_art1(2:end,:);...
                                   Idata(secmarker(ss)-pre_peak:...
                                   secmarker(ss)+post_peak)];
                            avg_art=mean(slice_art1,1);
                            art=2;
                        case 2
                            slice_art2=[slice_art2(2:end,:);...
                                    Idata(secmarker(ss)-pre_peak:...
                                    secmarker(ss)+post_peak)];
                            avg_art=mean(slice_art2,1);
                            art=1;
                    end
                end
                
            elseif strig==0  % Volume/Section Triggers
        
                if s==starts
                    art=1;
                    ssc=1;
                    for ss=starts+1:starts+2*hwinl+1
                        slice_art1(ssc,:)=...
                            Idata(secmarker(ss)-pre_peak:...
                            secmarker(ss)+post_peak);
                        ssc=ssc+1;
                    end
                    avg_art=mean(slice_art1,1); 
                elseif ((s>(starts+hwinl+2)) & (s<=(lasts-(hwinl+2))))
                    ss=s+hwinl;
                    slice_art1=[slice_art1(2:end,:);...
                           Idata(secmarker(ss)-pre_peak:...
                           secmarker(ss)+post_peak)];   
                    avg_art=mean(slice_art1,1);
                end
            end
       		
            % For first channel, find shift in artifact position to minimise 
            % sum of squared error between data and artifact template
            % - Assume same shift applies for all channels-
            % Also calculate Scale factor 'Alpha' to minimize sum of 
            % squared error
            ppn=1;
	        if s==1
                try
                    if c==1
                        for B=secmarker(s)-searchw:secmarker(s)+searchw
                                C(ppn)=prcorr2(Idata(B-pre_peak:B+post_peak),...
                                    avg_art);
                            ppn=ppn+1;
                        end
                        [CV,CP]=max(C);
						Beta=CP-(searchw+1);
						secmarker2(s)=secmarker(s)+Beta;
                    end
                   if isempty(intersect(exc,c))
                        Alpha=sum(Idata(secmarker2(s)-pre_peak:secmarker2(s)+...
                            post_peak).*avg_art)/sum(avg_art.*avg_art);
                    else
                       Alpha=1;
                    end
                        
                    INoise(secmarker2(s)-pre_peak:secmarker2(s)+...
                        post_peak)=Alpha*avg_art;
                catch
                    if c==1
                        if sec==1
                        warning...
                            ('Not enough data to remove first artifact segment');
                        STARTFLAG=1;
                        end
                    end
                end
            elseif sec==sections & s==lasts
                try
                    if c==1
                        for B=secmarker(s)-searchw:secmarker(s)+searchw
                                C(ppn)=prcorr2(Idata(B-pre_peak:B+post_peak),...
                                    avg_art);
                            ppn=ppn+1;
                        end
                        [CV,CP]=max(C);
						Beta=CP-(searchw+1);
						secmarker2(s)=secmarker(s)+Beta;
                    end
                    if isempty(intersect(exc,c))
                        Alpha=...
                            sum(Idata(secmarker2(s)-pre_peak:...
                            secmarker2(s)+post_peak).*avg_art)/...
                                sum(avg_art.*avg_art);
                    else
                       Alpha=1;
                    end
                        
                    INoise(secmarker2(s)-pre_peak:secmarker2(s)+post_peak)...
                        =Alpha*avg_art;
                catch
                    if c==1
                    warning('Not enough data to remove last artifact segment');
                    LASTFLAG=1;
                    end
                end
            else
        		if c==1
                    for B=secmarker(s)-searchw:secmarker(s)+searchw
                            C(ppn)=prcorr2(Idata(B-pre_peak:B+post_peak),...
                                avg_art);
                        ppn=ppn+1;
                    end
                    [CV,CP]=max(C);
					Beta=CP-(searchw+1);
					secmarker2(s)=secmarker(s)+Beta;
                end
                if isempty(intersect(exc,c))
                    Alpha=sum(Idata(secmarker2(s)-pre_peak:...
                        secmarker2(s)+post_peak).*avg_art)/...
                        sum(avg_art.*avg_art);
                else
                   Alpha=1;
                end
                    
                INoise(secmarker2(s)-pre_peak:secmarker2(s)+post_peak)=...
                    Alpha*avg_art;
            end
           
            
            if s==starts+1
                c;
            end
        end
        
       
         %----------PCA of residuals-------------------
         fitted_res=zeros(length(INoise),1);
                  
         if isempty(intersect(exc,c)) & npc~=0
            Ipca=filtfilt(hpfwts,1,double(Idata-INoise));
            pccount=1;
            skcount=1;
            pick=cumsum(ones(markerl(sec),1)*2+round(rand(markerl(sec),1)));
            if strig~=1
                pick=[pick(1):pick(end)];
            end
            

            for s=starts+1:lasts-1
                % construct PCAMAT
                if skcount==pick(pccount)
                    pcamat(pccount,:)=...
                        Ipca(secmarker(s)-pre_peak:...
                        secmarker(s)+max_postpeak);
                    pccount=pccount+1;
                end
                skcount=skcount+1;
            end
            
            pcamat=detrend(pcamat','constant')';
            [apc,ascore,asvar]=pca_calc(pcamat(1:(pccount-1),:)');

            oev=100*asvar/sum(asvar);
            if sec==1
                if ischar(npc)
                    d_oev=find(abs(diff(oev))<TH_SLOPE);
                    dd_oev=diff(d_oev);
                    for I=1:length(dd_oev)-3
                        if [dd_oev(I) dd_oev(I+1) dd_oev(I+2)]==[1 1 1]
                            break
                        end
                    end
                    SLOPETH_PC=d_oev(I)-1;
                    TMPTH=find(cumsum(oev)>TH_CUMVAR);
                    CUMVARTH_PC=TMPTH(1);
                    TMPTH=find(oev<TH_VAREXP);
                    VAREXPTH_PC=TMPTH(1)-1;
                    pcs=floor(mean([SLOPETH_PC CUMVARTH_PC VAREXPTH_PC]));
                    fprintf('\n%d residual PCs will be removed from channel %d\n . If you get an error "line 746 of fmrib_fastr: index exceeds matrix dimensions" it means there is an inconsistency in your TR triggers, either the TR length or the number of markers',pcs,c);                    
                else
                    pcs=npc;
                end
            end

% TEST CODE
%             SPCS(c)=SLOPETH_PC;
%             CPCS(c)=CUMVARTH_PC;
%             VPCS(c)=VAREXPTH_PC;
%             PCS(c)=pcs;
            
           
            if strig==0
                papc=double([ascore(:,1:pcs) ones(pre_peak+max_postpeak+1,1)]);
            else
                papc=double([ascore(:,1:pcs)]);
            end
         

            minmax1=max(papc(:,1))-min(papc(:,1));
            for apc=2:pcs
                papc(:,apc)=papc(:,apc)*minmax1/...
                    (max(papc(:,apc))-min(papc(:,apc)));
            end

            for s=starts:lasts
                if s==1
                    if ~STARTFLAG
                       fitted_res(secmarker(s)-pre_peak:secmarker(s)+max_postpeak)=...
                           papc*(papc\...
                           double(Ipca(secmarker(s)-pre_peak:...
                           secmarker(s)+max_postpeak))');
                    end
                elseif s==lasts & sec==sections
                    if ~LASTFLAG
                        fitted_res(secmarker(s)-pre_peak:secmarker(s)+max_postpeak)=...
                           papc*(papc\...
                           double(Ipca(secmarker(s)-pre_peak:...
                           secmarker(s)+max_postpeak))');
                    end
                else
                    fitted_res(secmarker(s)-pre_peak:secmarker(s)+max_postpeak)=...
                           papc*(papc\...
                           double(Ipca(secmarker(s)-pre_peak:...
                           secmarker(s)+max_postpeak))');
                end
            end
            
         elseif strig==0 % not doing OBS and using volume triggers
            
            Ipca=Idata-INoise;
            papc=double(ones(pre_peak+max_postpeak+1,1));
            for s=starts:lasts
                if s==1
                    if ~STARTFLAG
                       fitted_res(secmarker(s)-pre_peak:secmarker(s)+max_postpeak)=...
                           papc*(papc\...
                           double(Ipca(secmarker(s)-pre_peak:...
                           secmarker(s)+max_postpeak))');
                    end
                elseif s==lasts & sec==sections
                    if ~LASTFLAG
                        fitted_res(secmarker(s)-pre_peak:secmarker(s)+max_postpeak)=...
                           papc*(papc\...
                           double(Ipca(secmarker(s)-pre_peak:...
                           secmarker(s)+max_postpeak))');
                    end
                else
                    fitted_res(secmarker(s)-pre_peak:secmarker(s)+max_postpeak)=...
                           papc*(papc\...
                           double(Ipca(secmarker(s)-pre_peak:...
                           secmarker(s)+max_postpeak))');
                end
            end
            
         end
        
        %-----------------end PCA Section------------------
        
        Idata=Iorig-INoise-fitted_res';
        
        if L > 1
            fcleanEEG=decimate2(Idata,L);
            fNoise=decimate2(INoise+fitted_res',L);
        else
            fcleanEEG=Idata;
            fNoise=INoise+fitted_res';
        end
           
        
        if sec==1
            if sections==1
                cleanEEG(d1:d2)=fcleanEEG;
                Noise(d1:d2)=fNoise;
            else
                cleanEEG(d1:(d1-1)+secl(sec))=fcleanEEG(1:end-pad_sec);
                Noise(d1:(d1-1)+secl(sec))=fNoise(1:end-pad_sec);
            end
        elseif sec==sections
            cleanEEG((d1-1)+sum(secl(1:sec-1))+1:d2)=...
                fcleanEEG(pad_sec+1:end);
            Noise((d1-1)+sum(secl(1:sec-1))+1:d2)=fNoise(pad_sec+1:end);
        else
            cleanEEG((d1-1)+sum(secl(1:sec-1))+1:(d1-1)+sum(secl(1:sec)))=...
                fcleanEEG(pad_sec+1:end-pad_sec);
            Noise((d1-1)+sum(secl(1:sec-1))+1:(d1-1)+sum(secl(1:sec)))=...
                fNoise(pad_sec+1:end-pad_sec);
        end  
         
        
        %Update progress bar
        scount=scount+1;
        percentdone=floor(scount*100/steps);
        if floor(percentdone)>=barth
            
            if percentdone>=25 & Flag25==0
                fprintf('25%% ')
                Flag25=1;
            elseif percentdone>=50 & Flag50==0
                fprintf('50%% ')
                Flag50=1;
            elseif percentdone>=75 & Flag75==0
                fprintf('75%% ')
                Flag75=1;
            elseif percentdone==100
                fprintf('100%%\n')
            else
                fprintf('.')
            end
            
            while barth<=percentdone
                barth=barth+barth_step;
            end
            if barth>100
                barth=100;
            end
        end
    end
    
    if lpf>0
        cleanEEG=filtfilt(lpfwts,1,double(cleanEEG));
        Noise=filtfilt(lpfwts,1,double(Noise));
    end
    
    if anc_chk==1
        % Adaptive Noise cancellation
        % ---------------------------
        refs=Noise(d1:d2)';
        tmpd=filtfilt(ANCfwts,1,double(cleanEEG))';
        d=double(tmpd(d1:d2));
        Alpha=sum(d.*refs)/sum(refs.*refs);
        refs=double(Alpha*refs);
        mu=double(0.05/(N*var(refs)));
        [out,y]=fastranc(refs,d,N,mu);
        if isinf(max(y))
            wst=sprintf('ANC Failed for channel number %d. Skipping ANC.',c);
            warning(wst);
            EEG.data(c,d1:d2)=cleanEEG(d1:d2);
        else
            EEG.data(c,d1:d2)=cleanEEG(d1:d2)-y';
        end
    else
        EEG.data(c,d1:d2)=cleanEEG(d1:d2);
    end
    
        scount=scount+1;
        percentdone=floor(scount*100/steps);
        if floor(percentdone)>=barth

            if percentdone>=25 & Flag25==0
                fprintf('25%% ')
                Flag25=1;
            elseif percentdone>=50 & Flag50==0
                fprintf('50%% ')
                Flag50=1;
            elseif percentdone>=75 & Flag75==0
                fprintf('75%% ')
                Flag75=1;
            elseif percentdone==100
                fprintf('100%%\n')
            else
                fprintf('.')
            end

            while barth<=percentdone
                barth=barth+barth_step;
            end
            if barth>100
                barth=100;
            end
        end      
end

% Calculate processing time
% --------------------------

mttoc=floor(toc/60);
sttoc=round(toc-mttoc*60);
if mttoc < 60
    fprintf('FASTR Finished in %d min %d sec.\n',mttoc,sttoc);
else
    httoc=floor(mttoc/60);
    mttoc=round(mttoc-httoc*60);
    fprintf('FASTR Finished in %d hrs %d min %d sec.\n',httoc,mttoc,sttoc);
end
return;
