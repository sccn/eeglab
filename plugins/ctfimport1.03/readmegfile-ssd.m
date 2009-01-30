%
%
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %
%      <                                                       > %  
%      <                      DISCLAIMER:                      > %
%      <                                                       > %
%      <  THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. > %
%      <  THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   > %
%      <                     OFFICIAL USE.                     > %
%      <                                                       > %
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<> %
%
%


function [read,sensorNames,sensorLocations,sensorOrientations,header] = readmegfile(folder,setup,sensorIndex,sensorInfo,CHAN,TIME,TRIALS);
%INPUTS--------------------------------------------------------------------
%
%   folder:   the directory and filename of the .ds data set that is
%             to be read. 
%
%   CHAN:   ex: [30:35] - an interval of the desired channels to be read.
%           If CHAN = 'eegsens', only eeg channels/sensorIndices will be read.
%           If CHAN = 'megsens', only meg channels/sensorIndices will be read.
%           If CHAN = 'refsens', only reference channels/sensorIndices will be read.
%           If CHAN = 'othersens' only the other channels/sensorIndices will be
%           read.
%
%   TIME:   ex. [0 5] - the desired time interval to be read.
%           If TRIALS = 'alltrials', the data for all of the trials will be
%           read.
%              
%   TRIALS: If TRIALS = n, the nth trial will be read.
%           If TRIALS = [3,5,8] (for example), trials 3,5, and 8 will be
%           read and read{1} = data for trial 3, read{2} = data for
%           trial 5, and read{3} = data for trial 8.
%           If TRIALS = [3:7] (for example), trials 3 through 7 will be
%           read.
%           If TIME = 'alltimes', the entire duration of the trial(s) will
%           be read (i.e. TIME = [1:setup.duration]).
%           
%OUTPUTS-------------------------------------------------------------------
%   read:   read contains all of the data. ex. type read{1} for it to
%           display the first set of data on the screen. 
%   
%   sensorNames: cell array of sensor names.
%   
%   sensorLocations: array of sensor locations for plotting.
%
%   header: used for writing new meg4 file.


%[setup,sensorIndex,sensorInfo] = readresfile(folder,1);    %------use this line if you only
%want to run readmegfile.m. This line will read .res4 information (i.e.
%setup, sensorIndex, and sensorInfo information) from readresfile.m.
cd(folder);
[path,rootname] = fileparts(folder);

                                         

dat4file = [rootname,'.meg4'];
[dat,message] = fopen(dat4file,'rb','s');
Time = linspace(-setup.pretrig/setup.sample_rate,setup.duration - (setup.pretrig/setup.sample_rate),setup.number_samples);

for n = 5:nargin
    chankey = eval('CHAN');
    timekey = eval('TIME');
    trialkey = eval('TRIALS');
          
    if strcmp(chankey,'megsens')
        CHAN = sensorIndex.megsens;
    end
    if strcmp(chankey,'refsens')
        CHAN = sensorIndex.refsens;
    end
    if strcmp(chankey,'eegsens')
        CHAN = sensorIndex.eegsens;
    end
    if strcmp(chankey,'othersens')
        CHAN = sensorIndex.othersens;
    end
    if strcmp(chankey,'allchans')
        CHAN = [1:setup.number_channels];
    end
    if strcmp(timekey,'alltimes')
        TIME = Time;
    end
    if strcmp(trialkey,'alltrials')
        TRIALS = [1:setup.number_trials];
    end
end

header = char(fread(dat,8,'char')');
% Apply gains and offsets
%setup.chan_names = char(setup.chan_names{:})
    channel_gain = zeros(setup.number_channels,1);
    channel_gain(sensorIndex.megsens) = ([sensorInfo(sensorIndex.megsens).proper_gain]'.*[sensorInfo(sensorIndex.megsens).q_gain]');
    channel_gain(sensorIndex.refsens) = ([sensorInfo(sensorIndex.refsens).proper_gain]'.*[sensorInfo(sensorIndex.refsens).q_gain]');
    channel_gain(sensorIndex.eegsens) = ([sensorInfo(sensorIndex.eegsens).q_gain]');
    channel_gain(sensorIndex.othersens) = ([sensorInfo(sensorIndex.othersens).q_gain]');
if nargin == 7
    n_trials = TRIALS;
        
    read = cell(length(TRIALS),1);
    
    
    %fseek(dat,8,-1);

    %smpl = sample*(TIME(end)-TIME(1));

    s_rt = 1/(Time(2) - Time(1));
    trial_size = 4*setup.number_channels*setup.number_samples;
    bet_trials = diff(n_trials)-1;
    small_trial = 4*(min(CHAN)-1)*setup.number_samples;
    large_trial = 4*(setup.number_channels-max(CHAN))*setup.number_samples; 
    
    if (max(TIME)-min(TIME))>setup.duration
        durat = setup.duration;
        disp('TIME input too large for trial...setting TIME= duration of trial...');
        pause(1.2);
        disp('...done');
        pause(1);
        fprintf('TIME= %g seconds',durat);
        drawnow
    else
        durat = (TIME(end)-TIME(1));
    end
 
    samples = round((durat)*s_rt)+1;
    intime = round((TIME(1)-Time(1))*s_rt)+1;
    channels = length([min(CHAN):max(CHAN)]);
    
    
    trl = 0;
    for trial = n_trials
        trl = trl+1;
        read{trl} = zeros(channels,samples);
        if trial == n_trials(1); %1st trial
            fseek(dat,(trial-1)*trial_size + small_trial + 4*(intime-1),0);
        else
            fseek(dat,(trial-n_trials(trl-1)-1)*trial_size + small_trial + large_trial,0);
        end
        %Read data
        %fseek(dat,ftell(dat) + 4,-1);
        read{trl} = fread(dat,[samples channels],[num2str(samples),'*int32=>int32'],4*(setup.number_samples-samples))';
        
        
        read{trl} = read{trl}(CHAN - min(CHAN)+1,:);
        
        read{trl} = diag(1./channel_gain(CHAN))*double(read{trl});
        
    end 

elseif nargin ~= 7 
    disp('Error: incorrect number of inputs');
end

for i=CHAN
    j=i-min(CHAN)+1;
    if length(sensorInfo(i).channel_names) <= 5
        sensorNames{j,1} = sensorInfo(i).channel_names;
    else 
        sensorNames{j,1} = sensorInfo(i).channel_names(1:5);
    end
    if ~isempty(sensorInfo(i).location)
        sensorLocations(j,:) = sensorInfo(i).location(:,1)';
        sensorOrientations(j,:) = sensorInfo(i).orientation(:,1)';
    end
end


