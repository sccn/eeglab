% eegdraw() - subroutine used by eegplotold() to plot data.
%
% Author: Colin Humphries, CNL, Salk Institute, La Jolla, 7/96

% Copyright (C) Colin Humphries, CNL, Salk Institute 7/96 from eegplot()
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

% 11-07-97 fix incorrect start times -Scott Makeig
% 01-25-02 reformated help & license -ad 

function y = eegdraw(fighandle)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract variables from figure and axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

userdata = get(fighandle,'UserData');
samplerate  = userdata(1);
PLOT_TIME   = userdata(2);
spacing_var = userdata(3);
time        = userdata(4);
maxtime     = userdata(5);
axhandle    = userdata(6);
plotcolor   = userdata(7);
disp_scale  = userdata(12);

colors = ['y','w'];

data = get(axhandle,'UserData');

if samplerate <=0
   fprintf('Samplerate too small! Resetting it to 1.0.\n')
   samplerate = 1.0;
end

if PLOT_TIME < 2/samplerate
   PLOT_TIME = 2/samplerate;
   fprintf('Window width too small! Resetting it to 2 samples.\n');
end

if time >= maxtime  % fix incorrect start time
   time = max(maxtime-PLOT_TIME,0);
   fprintf('Start time too large! Resetting it to %f.\n',time)
elseif time < 0
   time = 0;
   fprintf('Start time cannot be negative! Resetting it to 0.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define internal variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[chans,frames] = size(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Label x-axis - relabel the x-axis based on 
% the new value of time.
% Labels are placed at one-second intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cla % Clear figure

Xlab = num2str(time);
for j = 1:1:PLOT_TIME
   Q = num2str(time+j);
   Xlab = str2mat(Xlab, Q);
end
set (gca, 'Ytick', 0:spacing_var:chans*spacing_var)
 set (gca, 'XTickLabels', Xlab)   % needs TickLabels with s in Matlab 4.2
set (gca, 'XTick',(0:samplerate:PLOT_TIME*samplerate))

axis([0 PLOT_TIME*samplerate 0 (chans+1)*spacing_var]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:chans                 %repeat for each channel
   if (maxtime-time>PLOT_TIME)
      F = data(chans-i+1,(time*samplerate)+1:((time+PLOT_TIME)*samplerate));
   else
      F = data(chans-i+1,(time*samplerate)+1:(maxtime*samplerate));
   end
   F = F - mean(F) + i*spacing_var;
   plot (F,'clipping','off','color',colors(plotcolor))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot scaling I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if disp_scale == 1
   ps = PLOT_TIME*samplerate;
   sv = spacing_var;
   line([1.03*ps,1.03*ps],[1.5*sv 2.5*sv],'clipping','off','color','w')
   line([1.01*ps,1.05*ps],[2.5*sv,2.5*sv],'clipping','off','color','w')
   line([1.01*ps,1.05*ps],[1.5*sv,1.5*sv],'clipping','off','color','w')
   text(1.05*ps,2*sv,num2str(round(sv)),'clipping','off')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set slider=edit - reset the value of 
% the user controls so they agree.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = findobj(fighandle,'style','slider');
D = findobj(fighandle,'style','edit');
set (D, 'string', num2str(time));
