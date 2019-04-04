% eegplotgold() - display EEG data in a clinical format
%
% Usage:
% >> eegplotgold('dataname', samplerate, 'chanfile', 'title', yscaling, range)
%
% Inputs:
%   'dataname' - quoted name of a desktop global variable (see Ex. below)
%   samplerate - EEG sampling rate in Hz (0 -> default 256 Hz)
%   'chanfile' - file of channel info in topoplot() style
%                                        (0 -> channel numbers)
%   'title'    - plot title string       (0 -> 'eegplotgold()')
%   yscaling   - initial y scaling factor (0 - default is 300)
%   range      - how many seconds to display in window (0 -> 10)
%
% Note: this version of eegplotgold() reguires that your data matrix 
%       be defined as a global variable before running this routine. 
%
% Example:  >> global dataname
%           >> eegplotgold('dataname')
%
% Author: Colin Humphries, CNL, Salk Institute, La Jolla, 3/97
%
% See also: eegplot(), eegplotold(), eegplotsold()

% Copyright (C) Colin Humphries, CNL, Salk Institute 3/97 from eegplotold()
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

% 4-4-97 shortened name to eegplotgold() -sm
% 5-20-97 added read of icadefs.m for MAXEEGPLOTCHANS -sm
% 8-10-97 Clarified chanfile type -sm
% 01-25-02 reformated help & license, added links -ad 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = eegplotold(dataname, samplerate, channamefile, titleval, yscaling, range)

eval (['global ',dataname])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

icadefs;

% set initial spacing

eval(['DEFAULT_SPACING = max(max(',dataname,''')-min(',dataname,'''));'])

%  spacing_var/20 = microvolts/millimeter with 21 channels
%  for n channels: 21/n * spacing_var/20 = microvolts/mm
%  for clinical data.

DEFAULT_SAMPLERATE = 256;	% default rate in Hz.
DEFAULT_PLOTTIME = 10;		% default 10 second window
DEFAULT_TITLE = 'eegplotgold()';
errorcode=0;			    % initialize error indicator

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allow for different numbers of arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
   PLOT_TIME = DEFAULT_PLOTTIME;
else
   PLOT_TIME = range;
end
if nargin < 5
   spacing_var = DEFAULT_SPACING;
else
   spacing_var = yscaling;
end
if spacing_var == 0
   spacing_var = DEFAULT_SPACING;
end

if nargin < 4
   titleval = DEFAULT_TITLE;
end
if titleval == 0
   titleval = DEFAULT_TITLE;
end
if nargin < 3 
	channamefile = 0;
end
if nargin < 2
   samplerate = DEFAULT_SAMPLERATE;
end
if samplerate == 0,
	samplerate = DEFAULT_SAMPLERATE;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define internal variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SAMPLE_RATE = samplerate;
time = 0;
eval(['[chans,frames] = size(',dataname,');'])		%size of data matrix

maxtime = frames / samplerate;       %size of matrix in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the channel names 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if channamefile ~=0,		% read file of channel names
	chid = fopen(channamefile,'r');
	if chid <3,
		fprintf('plotdata: cannot open file %s.\n',channamefile);
		errorcode=2;
		channamefile = 0;
	else
		fprintf('Chan info file %s opened\n',channamefile);
	end

    icadefs;   % read MAXEEGPLOTCHANS from icadefs.m
	if errorcode==0,
		channames = fscanf(chid,'%s',[6 MAXEEGPLOTCHANS]);
		channames = channames';
    	[r c] = size(channames);
		for i=1:r
			for j=1:c
				if channames(i,j)=='.',
					channames(i,j)=' ';
				end
			end
		end
		% fprintf('%d channel names read from file.\n',r);
		if (r>chans)
			fprintf('Using first %d names.\n',chans);
			channames = channames(1:chans,:);
		end
		if (r<chans)
			fprintf('Only %d channel names read.\n',r);
		end
	end
  end
  if channamefile ==0, % plot channel numbers
	channames = [];
	for c=1:chans
		if c<10,
			numeric = ['   ' int2str(c)];	% four-character fields
		else
			numeric = ['  '  int2str(c)];
		end
		channames = [channames;numeric];
	end
  end; % setting channames

channames = char(channames, ' ');	% add padding element to Y labels


Xlab = num2str(time);
for j = 1:1:PLOT_TIME
   Q = num2str(time+j);
   Xlab = char(Xlab, Q);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Graph Characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   figure;				% plot a new figure	
   fighandle = gcf;
   orient landscape		% choose landscape printer mode
   hold on;
   set(gcf,'NumberTitle','off')
   if VERS >= 8.04
       set(gcf,'Name',['EEGPLOTOLD #',num2str(fighandle.Number)]);
   else
       set(gcf,'Name',['EEGPLOTOLD #',num2str(gcf)]);
   end
   set (gca, 'xgrid', 'on')				%Xaxis gridlines only
   set (gca, 'GridLineStyle','-')			%Solid grid lines
   set (gca, 'XTickLabels', Xlab)			%Use Xlab for tick labels
   set (gca, 'Box', 'on')				
   set (gca, 'XTick', time*samplerate:1.0*samplerate:PLOT_TIME*samplerate) 
   set (gca, 'Ytick', 0:spacing_var:chans*spacing_var)  % ytickspacing on channels
   set (gca, 'TickLength', [0.001 0.001])
   title(titleval)					% title is titleval
   axis([0 PLOT_TIME*samplerate 0 (chans+1)*spacing_var]);       % set axis values
   set (gca, 'YTickLabels', flipud(channames))   	% write channel names

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the selected EEG data epoch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:chans			
   if (maxtime-time>PLOT_TIME)  
      eval(['F = ',dataname,'(chans-i+1,(time*samplerate)+1:(time+PLOT_TIME*samplerate));'])
   else
      eval(['F = ',dataname,'(chans-i+1,(time*samplerate)+1:(maxtime*samplerate));'])
   end
   F = F - mean(F) + i*spacing_var;
   plot (F,'clipping','off')
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Scaling I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

line([(PLOT_TIME+.3)*samplerate,(PLOT_TIME+.3)*samplerate],[1.5*spacing_var 2.5*spacing_var],'clipping','off','color','w')
line([(PLOT_TIME+.3)*samplerate-10,(PLOT_TIME+.3)*samplerate+10],[2.5*spacing_var,2.5*spacing_var],'clipping','off','color','w')
line([(PLOT_TIME+.3)*samplerate-10,(PLOT_TIME+.3)*samplerate+10],[1.5*spacing_var,1.5*spacing_var],'clipping','off','color','w')
text((PLOT_TIME+.5)*samplerate,2*spacing_var,num2str(round(spacing_var)),'clipping','off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Control Routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   slider_position = [.125 .030 .3 .024];   % position of user-controlled slider
   edit_position = [.65 .025 .1 .05];       % position of edit box
   slider_position2 = [.8 .03 .1 .024];
   b1_position = [.175 .022 .09 .047];
   b2_position = [.29 .022 .09 .047];
   b3_position = [.125 .022 .045 .047];
   b4_position = [.385 .022 .045 .047];

   Max_Space = 1;
   Min_Space = DEFAULT_SPACING*2;
%   User_Data_Mat = [data;zeros(1,length(data))];

axhandle = gca;
   User_Data_Mat(1) = samplerate;
   User_Data_Mat(2) = PLOT_TIME;
   User_Data_Mat(3) = spacing_var;
   User_Data_Mat(4) = time;
   User_Data_Mat(5) = maxtime;
   User_Data_Mat(6) = axhandle;
   User_Data_Mat(7) = 1;  % color
   User_Data_Mat(8) = frames;
   User_Data_Mat(9) = chans;
   User_Data_Mat(12) = 1;

   tstring1 = 'data1973 = get(gcf,''UserData'');';
   tstring2 = 'set(gcf,''UserData'',data1973);';
   tstring3 = 'eegdrawgv(gcf);';

   TIMESTRING = [tstring1,'if (data1973(4)-data1973(2))<0;','data1973(4) = 0;','else;','data1973(4) = data1973(4) - data1973(2);','end;',tstring2,tstring3,'clear data1973'];

   hb = uicontrol('Style','PushButton','Units','Normalized','position',b1_position,'String','PREV','Callback',TIMESTRING);

   TIMESTRING = [tstring1,'if (data1973(4)+data1973(2))>=data1973(5);','data1973(4)=data1973(4);','else;','data1973(4) = data1973(4) + data1973(2);','end;',tstring2,tstring3,'clear data1973'];

   hf = uicontrol('Style','PushButton','Units','Normalized','position',b2_position,'String','NEXT','Callback',TIMESTRING);

   TIMESTRING = [tstring1,'if (data1973(4)-1)<0;','data1973(4) = 0;','else;','data1973(4) = data1973(4) - 1;','end;',tstring2,tstring3,'clear data1973'];

   hbos = uicontrol('Style','PushButton','Units','Normalized','position',b3_position,'String','<<','Callback',TIMESTRING);

   TIMESTRING = [tstring1,'if (data1973(4)+1)>=data1973(5);','data1973(4)=data1973(4);','else;','data1973(4) = data1973(4) + 1;','end;',tstring2,tstring3,'clear data1973'];

   hfos = uicontrol('Style','PushButton','Units','Normalized','position',b4_position,'String','>>','Callback',TIMESTRING);

   TIMESTRING = [tstring1,'time1973 = get(gco,''string'');','time1973 = str2num(time1973);','data1973(4) = time1973;',tstring2,tstring3,'clear time1973 data1973'];

   w=uicontrol('style','edit','units','normalized','HorizontalAlignment','left','position',edit_position,'UserData',axhandle,'callback',TIMESTRING);

TIMESTRING = [tstring1,'data1973(3) = get(gco,''value'');',tstring2,tstring3,'clear time1973 data1973'];

u=uicontrol('style','slider','units','normalized','position',slider_position2,'Max',Max_Space,'Min',Min_Space,'value',spacing_var,'UserData',axhandle,'callback',TIMESTRING);

%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up ui menus
%%%%%%%%%%%%%%%%%%%%%%%%%

%Window menu:

TIMESTRING = ['fighand1973 = gcf;','delete(fighand1973);','clear fighand1973;'];

fm1 = uimenu('Label','Window');
fm2 = uimenu(fm1,'Label','Close ','UserData',fighandle,'Callback',TIMESTRING);

%Display menu:

TIMESTRING = [tstring1,'out1973 = gettext(''Input new windowlength(sec).'');','if isempty(out1973);','out1973 = 0;','else;','data1973(2) = str2num(out1973);',tstring2,tstring3,'end;','clear data1973 out1973'];

dm1 = uimenu('Label','Display');
dm2 = uimenu(dm1,'Label',' Window Length','Interruptible','on','Callback',TIMESTRING);
dm3 = uimenu(dm1,'Label',' Color');

TIMESTRING = [tstring1,'data1973(7) = 1;',tstring2,tstring3,'clear data1973;'];

dm4 = uimenu(dm3,'Label','Yellow ','UserData',axhandle,'Interruptible','on','Callback',TIMESTRING);

TIMESTRING = [tstring1,'data1973(7) = 2;',tstring2,tstring3,'clear data1973;'];

dm5 = uimenu(dm3,'Label','White ','UserData',axhandle,'Interruptible','on','Callback',TIMESTRING);

TIMESTRING = ['label1973 = gettext(''Enter new title.'');','if isempty(label1973);','label1973 = 0;','else;','title(label1973);','end;','clear label1973;'];

dm6 = uimenu(dm1,'Label','Title ','Interruptible','on','Callback',TIMESTRING);

TIMESTRING = [tstring1,'Check1973 = get(data1973(10),''checked'');','if (Check1973(1:2) == ''on'');','set(data1973(10),''Checked'',''off'');','set(data1973(6),''XGrid'',''off'');','else;','set(data1973(10),''Checked'',''on'');','set(data1973(6),''XGrid'',''on'');','end;','clear data1973 Check1973;'];

dm7 = uimenu(dm1,'Label','Grid','Checked','on','Callback',TIMESTRING);

User_Data_Mat(10) = dm7;

TIMESTRING = [tstring1,'Check1973 = get(data1973(11),''checked'');','if (Check1973(1:2) == ''on'');','set(data1973(11),''Checked'',''off'');','data1973(12)= 0;','else;','set(data1973(11),''Checked'',''on'');','data1973(12) = 1;','end;',tstring2,tstring3,'clear data1973 Check1973;'];

dm8 = uimenu(dm1,'Label','Scaling I','Checked','on','Callback',TIMESTRING);

User_Data_Mat(11) = dm8;

%Settings menu:

sm1 = uimenu('Label','Settings');

TIMESTRING = [tstring1,'Srate1973 = gettext(''Enter new samplerate'');','if isempty(Srate1973);','Srate1973 = 0;','else;','data1973(1) = str2num(Srate1973);','data1973(5) = data1973(8)/data1973(1);','data1973(4) = 0;',tstring2,tstring3,'end;','clear Srate1973 data1973'];

sm2 = uimenu(sm1,'Label','Samplerate','Interruptible','on','Callback',TIMESTRING);

%Electrodes menu:

em1 = uimenu('Label','Electrodes');

TIMESTRING = [tstring1,'ChanNamefile1973 = gettext(''Enter Electrode file to load.'');','if isempty(ChanNamefile1973);','ChanNamefile1973=0;','else;','ChanNames1973 = loadelec(ChanNamefile1973);','set(data1973(6),''YTickLabels'',flipud(ChanNames1973));','end;','clear data1973 ChanNamefile1973 ChanNames1973'];

em2 = uimenu(em1,'Label','Load Electrode File ','Interruptible','on','Callback',TIMESTRING);

TIMESTRING = [tstring1,'ChanNames1973 = makeelec(data1973(9));','if isempty(ChanNames1973);','ChanNames1973=0;','else;','set(data1973(6),''YTickLabels'',flipud(ChanNames1973));','end;','clear data1973 ChanNames1973'];

em3 = uimenu(em1,'Label','Make Electrode File ','Interruptible','on','Callback',TIMESTRING);

set(axhandle,'UserData',dataname)
set(fighandle,'UserData',User_Data_Mat)
