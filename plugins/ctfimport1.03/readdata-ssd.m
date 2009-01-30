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

function [read,setup,sensorNames,sensorLocations,sensorOrientations] = readdata(folder,CHAN,TIME,TRIALS);

[setup,sensorIndex,sensorInfo] = readresfile(folder);

if nargin == 1
    CHAN = 'allchans';
    TIME = 'alltimes';
    TRIALS = 'alltrials';
end

[read,sensorNames,sensorLocations,sensorOrientations] = readmegfile(folder,setup,sensorIndex,sensorInfo,CHAN,TIME,TRIALS);



% Use function [read,setup,sensorNames,sensorLocations] = readdata(folder); to read in parts 
% or all of data in a dataset. readdata.m will access readresfile.m to read 
% in header, gain/offset, and channel information. It will access 
% readmegfile.m to read in the data. 
% EXAMPLE:
% [read,sensorNames] = readdata('/data/directory/datasetname.ds');
% 
% Enter channel/time/trial parameters in on line 5 of readdata.m.
% FOR THE FUNCTION [read] = readmegfile(folder,setup,sensor,gain,'megsens','alltimes','alltrials');
% 
% INPUTS---------------------------------------------------------------------
% folder:     In readdata.m, this should remain as 'folder' everywhere in the
%             program. It is in the command line of matlab where you indicate 
%             the file and directory of the dataset you wish to view. It will
%             automatically be used in both readresfile.m and readmegfile.m.
% 
% setup:      Enter as 'setup'; setup header information acquired by readresfile.m.
%             This input is needed for the program to work.
% 
% sensorIndex:     Enter as 'sensorIndex'; sensor information acquired by readresfile.m,
%             such as megsens, othersens, refsens, and eegsens.
% 
% sensorInfo:       Emter as 'gain'; gain/offset information acquired by readresfile.m,
%             including proper_gain, io_gain, and channel locations.
%             
% CHAN:       ex: [30:35] - an interval of the desired channels to be read.
%             If CHAN = 'eegsens', only eeg channels/sensors will be read.
%             If CHAN = 'megsens', only meg channels/sensors will be read.
%             If CHAN = 'refsens', only reference channels/sensors will be read.
%             If CHAN = 'othersens' only the other channels/sensors will be
%             read.
% 
% TIME:       ex. [0 5] - seconds 0 to 5: the desired time interval to be read.
%             If TIME = 'alltimes', the entire duration of the trial(s) will
%             be read (i.e. TIME = [1:setup.duration]).
%               
% TRIALS:     If TRIALS = n, the nth trial will be read.
%             If TRIALS = [3,5,8] (for example), trials 3,5, and 8 will be
%             read and read{1} = data for trial three, read{2} = data for
%             trial 5, and read{3} = data for trial 8.
%             If TRIALS = [3:7] (for example), trials 3 through 7 will be
%             read.            
%             If TRIALS = 'alltrials', the data for all of the trials will be
%             read (i.e. TRIALS = [1:setup.duration]).
%         
% OUTPUTS--------------------------------------------------------------------
% read:   read contains all of the data. ex. type read{1} for it to display the first
%         set of data on the screen.
%
%----------------------------------------------------------------------------
%
% EXAMPLE:
% [read] = readmegfile(folder,setup,sensorIndex,sensorInfo,'megsens',[0 4],'alltrials');




