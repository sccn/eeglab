%
%
%
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %
%      <                                                      > %  
%      <                    DISCLAIMER:                       > %
%      <                                                      > %
%      <  THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. > %
%      <  THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   > %
%      <                    OFFICIAL USE.                     > %
%      <                                                      > %
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<> %
%
%
%

function [setup,sensorIndex,sensorInfo] = readresfile(folder,SETUP);

%INPUTS--------------------------------------------------------------------
%
%   folder:   the directory and filename of the .ds data set that is
%             to be read. If only the folder is inputted, SETUP = 1.
%
%   SETUP:  If SETUP = 0, when the program is run, 'setup' structure information will
%           not be displayed automatically. 
%           If SETUP = 1, when the program is run, 'setup' structure information will
%           automatically be displayed. 
%
%OUTPUTS-------------------------------------------------------------------
%   setup:  a structure of head information consisting of date, time, run
%           name, run title, subject, run description, operator, number of
%           channels, number of samples, sample rate, number of trials,
%           duration, pretrigger points, sensor filename, head zeroing,
%           and number of filters.
%   
%   sensorIndex: a structure with sensor information consisting of EEG sensors,
%           MEG sensors, reference sensors, and other sensors.
%   
%   sensorInfo:   a structure with gain and offset information consisting of
%           proper gain, Q gain, io gain, io offset, and index.

if nargin == 1
    SETUP=1;
end
cd(folder);
[path,rootname] = fileparts(folder);
fid4file = [rootname,'.res4'];
[fid,message] = fopen(fid4file,'rb','s');


%DATE/TIME-----------------------------------------------------------------
fseek(fid,778,-1);
tim = char(fread(fid,255,'char'));
fseek(fid,1033,-1);
dat = char(fread(fid,255,'char'));
time = [tim]';
date = [dat]';

if SETUP
    fprintf('\n Date Collected: %s',dat);
    fprintf('\n Time Collected: %s',tim);
    drawnow
end

%RUN NAME-----------------------------------------------------------------
fseek(fid,1360,-1);
r_name = char(fread(fid,32,'char'));
run_name = [r_name]';
if SETUP 
    fprintf('\n Run Name: %s',r_name);
    drawnow
end

%RUN TITLE----------------------------------------------------------------
fseek(fid,1392,-1);
r_title = char(fread(fid,256,'char'));
run_title = [r_title]';
if SETUP 
    fprintf('\n Run Title: %s',r_title);
    drawnow
end

%SUBJECT------------------------------------------------------------------
fseek(fid,1712,-1);
subj = char(fread(fid,32,'char'));
subject = [subj]';
if SETUP 
    fprintf('\n Subject: %s',subj);
    drawnow
end
    

%RUN DESCRIPTION-----------------------------------------------------------
fseek(fid,1836,-1);
setup = fread(fid,1,'int32');
fseek(fid,1844,-1);
run_description = char(fread(fid,setup,'char'));
if SETUP 
    fprintf('\n Run Description: %s',run_description);
    drawnow
end

%OPERATOR-----------------------------------------------------------------
fseek(fid,1744,-1);
oper = char(fread(fid,32,'char'));
operator = [oper]';
if SETUP 
    fprintf('\n Operator: %s',oper);
    drawnow
end

%NUMBER OF CHANNELS--------------------------------------------------------
fseek(fid,1292,-1);
number_channels = (fread(fid,1,'int16')');
if SETUP 
    fprintf('\n Number of Channels: %d',number_channels);
    drawnow
end

%NUMBER OF SAMPLES---------------------------------------------------------
fseek(fid,1288,-1);
number_samples = (fread(fid,1,'int32')');
if SETUP 
    fprintf('\n Number of Samples: %d',number_samples);
    drawnow
end

%SAMPLE RATE---------------------------------------------------------------
fseek(fid,1296,-1);
sample_rate = fread(fid,1,'double')';
if SETUP
    fprintf('\n Sample Rate: %g Samples per Second',sample_rate);
    drawnow
end

%NUMER OF TRIALS-----------------------------------------------------------
fseek(fid,1312,-1);
number_trials = fread(fid,1,'int16')';
if SETUP
    fprintf('\n Number of Trials: %g Trials',number_trials);
    drawnow
end

%DURATION------------------------------------------------------------------
fseek(fid,1304,-1);
dur = fread(fid,1,'double')';
duration = dur/number_trials;
if SETUP 
    fprintf('\n Duration: %g Seconds per Trial',duration);
    drawnow
end

%PRE_TRIG POINTS-----------------------------------------------------------
fseek(fid,1316,-1);
pretrig = fread(fid,1,'int32')';
if SETUP 
    fprintf('\n Pre-Trigger Points: %g Samples',pretrig);
    drawnow
end

%SENSOR FILE NAME----------------------------------------------------------
fseek(fid,1776,-1);
s_file_name = char(fread(fid,60,'char'));
sensor_file_name = [s_file_name]';
if SETUP 
    fprintf('\n Sensor File Name: %s',s_file_name);
    drawnow
end

%HEAD ZEROING--------------------------------------------------------------
fseek(fid,1348,-1);
h_zero = fread(fid,1,'int32')';
no_yes = {'no','yes'};
head_zero = no_yes{h_zero+1};
if SETUP 
    fprintf('\n Head Zeroing: %s',no_yes{h_zero+1});
    drawnow
end

%FILTERS-------------------------------------------------------------------
fseek(fid,1836,-1);
nfsetup = fread(fid,1,'int32')';

fseek(fid,1844,-1);
char(fread(fid,nfsetup,'char'));
number_filters = fread(fid,1,'int16');
if SETUP 
    fprintf('\n Number of Filters: %d\n',number_filters);
    drawnow
end
for jk = 1:number_filters
    freq=fread(fid,1,'double');
    class=fread(fid,1,'int32');
    type=fread(fid,1,'int32');
    numparam=fread(fid,1,'int16');
    params=fread(fid,numparam,'double');
    if SETUP
        fprintf(' -Filter # \t%g\n',jk)
        fprintf(' -Frequency: \t%g Hz\n',freq)
        fprintf(' -Class: \t%g\n',class)
        fprintf(' -Type: \t%g\n',type)
        if ~isempty(params)
            fprintf(' -Parameter(s): \t%g\n',params)
        end
        drawnow
    end    
end       
        
setup = struct('date',{date},'time',{time},'run_name',{run_name},'run_title',{run_title},'operator',{operator}...
    ,'number_channels',{number_channels},'number_samples',{number_samples},'sample_rate',{sample_rate},...
    'number_trials',{number_trials},'duration',{duration},'run_description',{run_description},'number_filters'...
    ,{number_filters},'subject',{subject},'pretrig',{pretrig},'sensor_file_name',{sensor_file_name},'head_zero',{head_zero});
sensorInfo = struct('proper_gain',[],'q_gain',[],'io_gain',[],'io_offset',[],'index',[],'extra',[],'channel_names',[]);

for chan = 1:number_channels 
    temp = fread(fid,32,'char');
    temp(temp>127) = 0;
    temp(temp<0) = 0;
    temp = strtok(temp,char(0));
    sensorInfo(chan).channel_names = char(temp');
end


Time = linspace(-pretrig/sample_rate,duration - (pretrig/sample_rate),number_samples);


for chan = 1:number_channels
    ftell(fid);
    sensorInfo(chan).index = fread(fid,1,'int16');
    sensorInfo(chan).extra = fread(fid,1,'int16');
    id = fread(fid,1,'int32')+1;
    sensorInfo(chan).proper_gain = fread(fid,1,'double');
    sensorInfo(chan).q_gain = fread(fid,1,'double');
    sensorInfo(chan).io_gain = fread(fid,1,'double');
    sensorInfo(chan).io_offset = fread(fid,1,'double');
    fread(fid,1,'int16');
    fread(fid,1,'int16');
    fread(fid,1,'int32');
    %fseek(fid,ftell(fid)+6,0);
  
    for pos = 1:8
        sensorInfo(chan).coil(pos).position.x = fread(fid,1,'double');
        sensorInfo(chan).coil(pos).position.y = fread(fid,1,'double');
        sensorInfo(chan).coil(pos).position.z = fread(fid,1,'double');
        fread(fid,1,'double');
        sensorInfo(chan).coil(pos).orient.x = fread(fid,1,'double');
        sensorInfo(chan).coil(pos).orient.y = fread(fid,1,'double');
        sensorInfo(chan).coil(pos).orient.z = fread(fid,1,'double');
        fread(fid,1,'double');
        fread(fid,1,'int16');
        fread(fid,1,'int32');
        fread(fid,1,'int16');
        fread(fid,1,'double');
        
        %fseek(fid,ftell(fid)+56,0);
        %fseek(fid,ftell(fid)-80,0);
    end
    
    for pos = 1:8
        sensorInfo(chan).hcoil(pos).position.x = fread(fid,1,'double');
        sensorInfo(chan).hcoil(pos).position.y = fread(fid,1,'double');
        sensorInfo(chan).hcoil(pos).position.z = fread(fid,1,'double');
        fread(fid,1,'double');
        sensorInfo(chan).hcoil(pos).orient.x = fread(fid,1,'double');
        sensorInfo(chan).hcoil(pos).orient.y = fread(fid,1,'double');
        sensorInfo(chan).hcoil(pos).orient.z = fread(fid,1,'double');
        fread(fid,1,'double');
        fread(fid,1,'int16');
        fread(fid,1,'int32');
        fread(fid,1,'int16');
        fread(fid,1,'double');
        
        %fseek(fid,ftell(fid)+56,0);
        %fseek(fid,ftell(fid)+80,0);
    end
    %fseek(fid,ftell(fid)+1288,-1);
end

megsens = find([sensorInfo.index] == 5);
eegsens = find([sensorInfo.index] == 9);
refsens = find([sensorInfo.index] == 0); 
refsens = [refsens,find([sensorInfo.index] == 1)];
a = [megsens,eegsens,refsens];
othersens = setdiff([1:setup.number_channels],a);
sensorIndex = struct('refsens',{refsens},'othersens',{othersens},'eegsens',{eegsens},'megsens',{megsens});    

%Location Coordinates of Channels in centimeters
for chan = 1:number_channels
    switch sensorInfo(chan).index
        %0=Reference Channels, 1=More Reference Channels, 5=MEG Channels  
        case {0,1,5}
            coord = [sensorInfo(chan).hcoil(:).position];
            coordx = [coord(1:2).x];
            coordy = [coord(1:2).y];
            coordz = [coord(1:2).z];
            sensorInfo(chan).location = [coordx;coordy;coordz];
            orient = [sensorInfo(chan).hcoil(:).orient];
            orientx = [orient(1).x];
            orienty = [orient(1).y];
            orientz = [orient(1).z];
            sensorInfo(chan).orientation = [orientx;orienty;orientz];
            location =  sensorInfo(chan).location';
            ps1=sum((location(1,:).*sensorInfo(chan).orientation(:,1)')')';
            sensorInfo(chan).orientation(:,1) = (sign(ps1)*ones(1,3).*sensorInfo(chan).orientation(:,1)')';
        %EEG Channels    
        case 9
            coord = [sensorInfo(chan).hcoil(:).position];
            coordx = [coord(1:2).x];
            coordy = [coord(1:2).y];
            coordz = [coord(1:2).z];
            sensorInfo(chan).location = [coordx;coordy;coordz];
    end         
end 
