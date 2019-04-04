% eegdrawg() - subroutine used by eegplotgold() to plot data.
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

% 4-4-97 shortened name to eegdrawq() -sm
% 4-7-97 allowed data names other than 'data' -ch 
% 01-25-02 reformated help & license -ad 

function y = eegdrawg(fighandle)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract variables from figure and axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

userdata = get(fighandle,'UserData');
samplerate = userdata(1);
PLOT_TIME = userdata(2);
spacing_var = userdata(3);
time = userdata(4);
maxtime = userdata(5);
axhandle = userdata(6);
plotcolor = userdata(7);
disp_scale = userdata(12);
colors = ['y','w'];
dataname = get(axhandle,'UserData');
eval(['global ',dataname])

cla % Clear figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define internal variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval(['[chans,frames] = size(',dataname,');']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Label x-axis
% This routine relabels the x-axis based on the new value of time.
% Labels are placed on one second intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xlab = num2str(time);
for j = 1:1:PLOT_TIME
   Q = num2str(time+j);
   Xlab = str2mat(Xlab, Q);
end
set (gca, 'Ytick', 0:spacing_var:chans*spacing_var)
set (gca, 'XTickLabels', Xlab)
set (gca, 'XTick',(0:samplerate:PLOT_TIME*samplerate))

axis([0 PLOT_TIME*samplerate 0 (chans+1)*spacing_var]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:chans                 %repeat for each channel
   if (maxtime-time>PLOT_TIME)
      eval(['F = ',dataname,'(chans-i+1,(time*samplerate)+1:((time+PLOT_TIME)*samplerate));'])
   else
      eval(['F = ',dataname,'(chans-i+1,(time*samplerate)+1:(maxtime*samplerate));'])
   end
   F = F - mean(F) + i*spacing_var;
   plot (F,'clipping','off','color',colors(plotcolor))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Scaling I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if disp_scale == 1
   line([(PLOT_TIME+.3)*samplerate,(PLOT_TIME+.3)*samplerate],[1.5*spacing_var 2.5*spacing_var],'clipping','off','color','w')
   line([(PLOT_TIME+.3)*samplerate-10,(PLOT_TIME+.3)*samplerate+10],[2.5*spacing_var,2.5*spacing_var],'clipping','off','color','w')
   line([(PLOT_TIME+.3)*samplerate-10,(PLOT_TIME+.3)*samplerate+10],[1.5*spacing_var,1.5*spacing_var],'clipping','off','color','w')
   text((PLOT_TIME+.5)*samplerate,2*spacing_var,num2str(round(spacing_var)),'clipping','off')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set slider=edit
% This routine resets the value of the user controls so that
%    they agree.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   H = findobj(fighandle,'style','slider');
   D = findobj(fighandle,'style','edit');
   set (D, 'string', num2str(time));
