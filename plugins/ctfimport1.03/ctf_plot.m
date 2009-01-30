function ctf_plot(ctf,CHAN,TIME,TRIALS,Xhair)

% ctf_plot - plot ctf.data
%
% ctf_plot(ctf,CHAN,TIME,TRIALS,Xhair)
%
% CHAN - see ctf_channel_select for options
% TIME - see ctf_read for options (given in msec)
% TRIALS - select 1 trial to plot (the default is trial = 1)
%
% Xhair - 1 to use the crosshair GUI features (default), 0 otherwise
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%      <                                                      >
%      <                    DISCLAIMER:                       >
%      <                                                      >
%      < THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. >
%      < THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   >
%      <                    OFFICIAL USE.                     >
%      <                                                      >
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>
%

% $Revision: 1.1 $ $Date: 2009-01-30 03:49:27 $

% Copyright (C) 2004  Darren L. Weber
% 
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

% Modified: 02/2004, Darren.Weber_at_radiology.ucsf.edu
%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~exist('CHAN','var'), CHAN = 'all'; end
if ~exist('TIME','var'), TIME = 'all'; end
if ~exist('TRIALS','var'), TRIALS = 1; end
if ~exist('Xhair','var'), Xhair = 1; end


if isempty(CHAN), CHAN = 'all'; end
if isempty(TIME), TIME = 'all'; end
if isempty(TRIALS), TRIALS = 1; end
if isempty(Xhair), Xhair = 1; end


% This function calls ctf_channel_select
[CHAN,type] = ctf_channel_select(ctf,CHAN);

switch num2str(TIME),
  case 'all',
    TIME = ctf.setup.time_msec;
    TIME_index = 1:ctf.setup.number_samples;
otherwise
    fprintf('...sorry, time restrictions not implemented correctly (03/2004)\n');
    TIME = ctf.setup.time_msec;
    TIME_index = 1:ctf.setup.number_samples;
    
    
    %     % assume the input is a range of times in sec
    %     % check the range
    %     if TIME(1) > ctf.setup.time_msec(1),
    %       fprintf('...setting TIME(1) = ctf.setup.time_msec(1)\n');
    %       TIME(1) = ctf.setup.time_msec(1);
    %     end
    %     if TIME(end) > ctf.setup.time_msec(end),
    %       fprintf('...setting TIME(end) = ctf.setup.time_msec(end)\n');
    %       TIME(end) = ctf.setup.time_msec(end);
    %     end
    %     % now find the nearest indices into the samples matrix
    %     TIME_index = interp1(ctf.setup.time_msec,1:ctf.setup.number_samples,TIME,'nearest');
    %     % now ensure that the TIME array is consistent with ctf.setup.time_sec
    %     TIME = ctf.setup.time_msec(TIME_index);
end
TIME = sort(TIME);


switch num2str(TRIALS),
  case 'all',
    TRIALS = 1:ctf.setup.number_trials;
  otherwise
    % assume the input is an array of trials
end
TRIALS = unique(sort(TRIALS));

% check the input data for more than 1 trial
trials = size(ctf.data,3);
if trials > 1,
    fprintf('...plotting trial %d of %d trials in ctf.data\n',TRIALS,trials);
    data = ctf.data(:,:,TRIALS);
else
    data = ctf.data;
end


% now plot the data
switch type,
    
    case 'all',
        
        if isempty(ctf.sensor.index.meg_sens),
            fprintf('...no meg sensors\n');
        else
            figure('Name','MEG sensors');
            MEGCHAN = ctf_channel_select(ctf,'meg');
            plot(ctf.setup.time_msec(TIME_index),data(TIME_index,MEGCHAN))
            axis tight
            if exist('crosshair') & Xhair,
                crosshair;
            end
        end
        
        if isempty(ctf.sensor.index.eeg_sens),
            fprintf('...no eeg sensors\n');
        else
            figure('Name','EEG sensors');
            EEGCHAN = ctf_channel_select(ctf,'eeg');
            plot(ctf.setup.time_msec(TIME_index),data(TIME_index,EEGCHAN))
            axis tight
        end
        
        
    case 'meg',
        
        if isempty(CHAN),
            fprintf('...no meg sensors\n');
        else
            figure('Name','MEG sensors');
            plot(ctf.setup.time_msec(TIME_index),data(TIME_index,CHAN))
            axis tight
        end
        
    case 'eeg',
        
        if isempty(CHAN),
            fprintf('...no eeg sensors\n');
        else
            figure('Name','EEG sensors');
            plot(ctf.setup.time_msec(TIME_index),data(TIME_index,CHAN))
            axis tight
        end
        
    case 'ref',
        
        if isempty(CHAN),
            fprintf('...no ref sensors\n');
        else
            figure('Name','MEG REF sensors');
            plot(ctf.setup.time_msec(TIME_index),data(TIME_index,CHAN))
            axis tight
        end
        
    case 'others',
        
        if isempty(CHAN),
            fprintf('...no other sensors\n');
        else
            figure('Name','Other sensors');
            plot(ctf.setup.time_msec(TIME_index),data(TIME_index,CHAN))
            axis tight
        end
        
    otherwise
        
        figure;
        plot(ctf.setup.time_msec(TIME_index),data(TIME_index,CHAN))
        axis tight
        
end

if exist('crosshair') & Xhair,
    crosshair;
end

return
