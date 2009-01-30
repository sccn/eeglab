% BATCH for data processing
% this is a TEMPLATE 

%	(C) 2004 by Alois Schloegl <a.schloegl@ieee.org>	


%--> get list of files:
fn = dir('f5*.cnt');    % this way 
fn = '\\Fdpmias009\erd-ers\4class_offline_BCI\s1_200104\*.cnt';    % or this way
fn = 'L:\BCI7B\s74\*.bkr';    % or this way

% load data
[s,HDR] = sload(fn);
%[s,HDR] = sload(fn,1:5);       % you can select certain channels, too. 

% spatial filtering e.g. Laplace or LAR
%--> optional 
%[s,M] = laplace(s,montage);   


% get number of classes
CL = [];
if HDR.EVENT.N > 0;
        HDR.Classlabel = unique(HDR.EVENT.TYP);
end;
if isfield(HDR,'Classlabel');
        CL = sort(unique(HDR.Classlabel));
end;

% get trigger information
if ~HDR.FLAG.TRIGGERED,
        if HDR.EVENT.N>0,
                TRIG = HDR.EVENT.POS;
        else
%--> select trigger channel 
                TriggerChannel = 4;             
                TRIG = gettrigger(s(:,TriggerChannel));
                s = s(:,[1:TriggerChannel-1,TriggerChannel+1:HDR.NS]);
        end
        if length(TRIG)<2,
                error('size of Artifact selection and size of Classlabels do not fit\n');
        end;                
%--> select start and end sample
        gap   =  0;
        start = -2*HDR.SampleRate;
        final =  6*HDR.SampleRate;
        T = (start:final)/HDR.SampleRate;
else
        T = (1:HDR.SPR)/HDR.SampleRate;
end;

% artifact selection 
if isfield(HDR,'ArtifactSelection')
        if length(HDR.ArtifactSelection)==length(HDR.Classlabel),
                HDR.Classlabel = HDR.Classlabel(~HDR.ArtifactSelection);
        else
                error('size of Artifact selection and size of Classlabels do not fit\n');
        end;
        if length(HDR.ArtifactSelection)==length(TRIG),
                TRIG = TRIG(~HDR.ArtifactSelection);
        else
                error('size of Artifact selection and length of Trigger do not fit\n');
        end;
else
        warning('no info on artifact selection available\n');
end;

% trigger data
clear S;
if ~HDR.FLAG.TRIGGERED,
        if isempty(CL),
                [S{1},sz]=trigg(s,TRIG,start,final,gap);
                S{1} = reshape(S{1},sz); % convert into 3-dim matrix
        else                
                for k = 1:length(CL),
                        [S{k},sz]=trigg(s,TRIG(CL(k)==HDR.Classlabel),start,final,gap);
                        S{k} = reshape(S{k},sz); % convert into 3-dim matrix
                end;
        end;
else
        tmp  = reshape(s',[size(s,1),HDR.SPR,HDR.NRec]);
        if isempty(CL),
                S{1}=tmp;
        else
                for k = 1:length(CL),
                        S{k}=tmp(:,:,CL(k)==HDR.Classlabel);
                end;
        end;
end;

% data processing 
for k=1:length(CL),
        R{k} = statistic(S{k},3);
        R{k}.T = T;
        R{k}.datatype = 'MEAN+SEM';
end;


% display 
for k=1:6,hf(k)=subplot(2,3,k);end;        % generate subplots
for k=1:length(CL),
        plota(R{k},hf((-2:0)+k*3));
end;

