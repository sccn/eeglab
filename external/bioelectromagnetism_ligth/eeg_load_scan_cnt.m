function [x,y,tex,elec,cal,oth]=loadcnt(filen,elec,start,last);

% eeg_load_scan_cnt - Load continuous NeuroScan data.
%
%	[X,Y,LAB,ELEC,CAL,OTH]=LOADCNT(FILEN,ELEC,START,LAST) elec defines 
%       which electrodes are read, size - rows data in y. If elec=[] then 
%	reads all channels data between start and end points. If elec
%	is string then only matched are read and with elec corresponding
%	channelnumbers are returned. lab is lab(elec). Text oth is ['rev  '
%	;'date ';'time ';'nchan';'rate ';'datas';'Chann';'size';'TeegT']
%	read by scanhead.
%
%	filen can be filename or file identifier. Because function scanelec,
%	scanhead are slow, variables lab, cal and oth are saved into 
%	file.mat file in directory "Headers". If *.cnt file is changed you 
%	have to remove file.mat for new parameters to be read.
%
%DIAGNOSTICS
%	Perhaps sensitive to different SynAmps parameters. Used global 
%	variable is DIRS.
%
%REFERENCES
%	Neurosoft. Structure of data headers. 1993.
%
%SEE ALSO
%	Uses scanhead, scanelec, showwait.
%
%EXAMPLES
%	                % load and show first 1000 points
%	[x,y,la,el,ca,ot]=loadcnt('x70001a1.cnt',[],1,1000);	
%	figure
%	plotdata(x,y,la,[50 50]);ot,xlabel('t (s');ylabel('100 uV');
%	print Help\loadcnt -djpeg
%
%	                % test loading large samples
%	[x,y1]=loadcnt('testec','FZ',1,30010);
%	[x,y2]=loadcnt('testec','FZ',30000,30100);[y1(30000:30010),y2]

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $

%Mention source when using or modifying these Shareware tools
%JVIR, Jussi.Virkkala@occuphealth.fi
%JVIR,  7-Feb-1999 Saving into Headers location not Data Headers.
%JVIR,  3-Feb-1999 Modified for PCWIN Matlab 5.2.
%
%   J.Virkkala 30-Jun-94
%   J.Virkkala 12-Jan-95 Added parameter oth.
%   J.Virkkala 24-Jan-95 Saving header values into file.mat.
%   J.Virkkala 22-Feb-95 SynAmps working, Strings as file name.
%   J.Virkkala  3-Mar-95 Part of ScanUtil.
%	J.Virkkala 23-May-95 Still some problems with SynAmps.

global DIRS 
% only integer positions accepted
start=round(start);last=round(last);
%*** FILENUMBER OR NAME AS PARAMETER ***
if isstr(filen),
    file=fopen(filen,'r','l');
else
    file=filen;
end;
if file==-1,errorr('hff','loadcnt');break;end
name=[filename(file) '.mat'];
if exist(name),eval(['load ' name]);else 
    % load electrodes calibration
    showwait('loading calibrations - %1.0f');
    [a,a,d]=scanelec(['lab';'bas';'sen';'cal'],'',file);
    n=size(d,1);
    nchannels=n/5;
    lab=d(2:5:n,:);
    for i=1:nchannels;
        cal(i,1)=eval(d(i*5-2,:));  % baseline
        cal(i,2)=eval(d(i*5-1,:));  % sensitivity
        cal(i,3)=eval(d(i*5,:));    % calibration
    end;
    [a,a,oth]=scanhead(['rev  ';'date ';'time ';'ncha ';'rate ';'datas';...
            'Chann';'size ';'TeegT'],file);
    bs=10*eval(oth(4,:));
    oth=putstr(mat2str(floor(eval(oth(6,:))/bs)*bs-1),oth,6,1);
    if ~exist('ele'),ele=[];end   % created by eegevent
    eval(['save ''Headers' DIRS name ''' cal lab oth ele']);   
    showwait('');
end; 
nchannels=size(lab,1);
% if all electrodes
if isempty(elec);elec=1:nchannels;end;

% if want to look match
if isstr(elec);
    e=[];
    for i=1:size(elec,1);
        p=strfind(lab,deblank(elec(i,:)));
        % if found
        if size(p)~=[0,0];
            p=p(:,1);
            if size(e,1)>size(p,1);
                p(size(e,1),1)=0;
            end;
            if size(e,1)==0;
                %JVIR	   	p=p(find(p~=e)); % with same string possible diffent text
            end;
            p=p(1,1);
            size(e);
            if p;if size(e)==[0,0];e=p;else;e=[e;p];end;end;
        end;
    end;
    % names -> corresponding numbers
    elec=e';
else
    elec=elec(:)';
end;   

M=3000;y1=[];
%*** REPEATINGS FOR LONG SAMPLES ****
%JVIR, 4-Feb-1999
y=[];
x=[];
if (last-start)>M,
    for ind=start:M:last;
        e=ind+M-1;
        if e>last,e=last;end
        [x1,y1,tex,elec,cal,oth]=loadcnt(file,elec,ind,e);
        y=[y y1];
        x=[x x1];
    end
else
    
    l=eval(oth(6,:));
    if isempty(last),last=l;end
    if size(elec)==[0,0];return;end;
    % read data
    
    %*** NEUROSCAN SYNAMPS ***
    if eval(oth(7,:))>1310000,	% ==1310720
        bs=10*nchannels;    % block size
        st=floor(start/bs)*bs;
        la=ceil(last/bs)*bs;
        
        x=1:(last-start)+1;
        fseek(file,900+nchannels*(75+st*2),-1);
        Yall=fread(file,[nchannels*10,(la-st)/10],'short');
        % only wanted 
        y(nchannels,(la-st))=0;
        for i=1:size(elec,2),
            el=elec(i);
            a=Yall((el-1)*10+1:el*10,:);
            y(el,:)=a(:)';
        end
        s=x+start-st;
        %??? some problems with index ???
        if max(s)>size(y,2);s=s-1;end
        Yall=y(:,s);y=[];		% x+start-st
        I=ones(size(elec,2),size(x,2));
        y(1:size(elec,2),1:size(x,2))=sparse(diag(cal(elec,2)))*sparse(diag(cal(elec,3)))*...
            (Yall(elec,x)-sparse(diag(cal(elec,1)))*I)./204.8;
    else
        
        %*** NEUROSCAN FORMAT ***
        % fast but needs memory
        x=1:(last-start)+1;         
        fseek(file,900+nchannels*(75+(start-1)*2),-1);
        Yall=fread(file,[nchannels,(last-start+1)],'short');
        % only wanted channels          
        I=ones(size(elec,2),size(x,2));
        y(1:size(elec,2),1:size(x,2))=sparse(diag(cal(elec,2)))*sparse(diag(cal(elec,3)))*...
            (Yall(elec,x)-sparse(diag(cal(elec,1)))*I)./204.8;
    end
    
    x=x+ones(1,size(x,2))*(start-1);
    x=x/eval(oth(5,:)); 		% xlabel in seconds
    tex=lab(elec,:);
end


%*** FILENUMBER OR NAME AS PARAMETER ***
if isstr(filen),
    fclose(file);
end;

%END OF LOADCNT