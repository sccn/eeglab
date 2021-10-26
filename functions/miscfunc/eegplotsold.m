% eegplotsold() - display data in a clinical format without scrolling
%
% Usage:
%  >> eegplotsold(data, srate, 'chanfile', 'title', ...
%                           yscaling, epoch, linecolor,xstart,vertmark)
%
% Inputs:
%   data       - data matrix (chans,frames) 
%   srate      - EEG sampling rate in Hz (0 -> 256 Hz)
%   'chanfile' - file of channel info, topoplot() style, (0 -> chan nos)
%   'title'    - plot title string {0 -> 'eegplotsold()'}
%   yscaling   - initial y scaling factor (0 -> 300)
%   epoch      - how many seconds to display in window (0 -> 10 sec)
%   linecolor  - color of eeg (0 -> 'y')
%   xstart     - start time of the data {0 -> 0}
%   vertmark   - vector of frames to mark with vertical lines {0 -> none}
%
% Author: Colin Humphries, CNL, Salk Institute, La Jolla, 3/97
%
% See also: eegplot(), eegplotold(), eegplotgold()

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

% 05-01-97 added xstart argument -sm
% 05-20-97 added read of icadefs.m
% 06-12-97 EPOCH -> epoch line 71 below -sm
% 8-10-97 Clarified chanfile type -sm
% 12-08-97 Added ischar(titleval) test -sm
% 02-09-98 length(data)->size(data,2) -sm
% 01-25-02 reformated help & license, added links -ad 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = eegplotsold(data, srate, channamefile, titleval, yscaling, epoch, linecolor,xstart,vertmark)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

icadefs; % read MAXEEGPLOTCHANS, DEFAULT_SRATE, DEFAULT_EPOCH from icadefs.m

if nargin < 1
	help eegplotsold	% print usage message
    PLOT_TIME = 10;
	data = 0.5*randn(8,floor(DEFAULT_SRATE*PLOT_TIME*2));	
	titleval = ['eegplotsold() example - random noise'];	% show example plot
end

[chans,frames] = size(data);		%size of data matrix

% set initial spacing
DEFAULT_SPACING = max(max(data')-min(data'));
%  spacing_var/20 = microvolts/millimeter with 21 channels
%  for n channels: 21/n * spacing_var/20 = microvolts/mm
%  for clinical data.
DEFAULT_PLOTTIME = 10;  	% default 10 seconds per window
DEFAULT_TITLE = 'eegplotsold()';
errorcode=0;				% initialize error indicator

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allow for different numbers of arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 9,
    vertmark = 0;
end
if nargin < 8,
	xstart = 0;
end

if nargin < 7
   linecolor = 'y'; % DEFAULT LINECOLOR
end
if linecolor == 0,
   linecolor = 'y'; % DEFAULT LINECOLOR
end
if nargin < 6
   PLOT_TIME = 0;
else
   PLOT_TIME = epoch;
end
if nargin < 5
   spacing_var = 0;
else
   spacing_var = yscaling;
end
if spacing_var == 0
   spacing_var = DEFAULT_SPACING;
end

if nargin < 4
   titleval = DEFAULT_TITLE;
end
if ~ischar(titleval)
  if titleval == 0
    titleval = DEFAULT_TITLE;
  else
    help eegplotsold
    return
  end
end
if nargin < 3 
	channamefile = 0;
end
if nargin < 2
   srate = DEFAULT_SRATE;
end
if srate == 0,
	srate = DEFAULT_SRATE;
end
if PLOT_TIME == 0
   PLOT_TIME = ceil(frames/srate);
   if PLOT_TIME > DEFAULT_EPOCH
     PLOT_TIME = DEFAULT_EPOCH;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define internal variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxtime = frames / srate;       %size of matrix in seconds

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
		% fprintf('Chan info file %s opened\n',channamefile);
	end
	if errorcode==0,
		channames = fscanf(chid,'%d %f %f %s',[7 MAXEEGPLOTCHANS]);
		channames = channames';
        channames = setstr(channames(:,4:7)); % convert ints to chars
    	[r c] = size(channames);
		for i=1:r     
			for j=1:c
				if channames(i,j)=='.',
					channames(i,j)=' '; % convert dots to spaces
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make matrix of x-tick labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Xlab = num2str(0);
% for j = 1:1:PLOT_TIME
%    Q = num2str(0+j);
%    Xlab = char(Xlab, Q);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Graph Characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   hold on;
   set (gca, 'xgrid', 'on')				    % Xaxis gridlines only
   set (gca, 'GridLineStyle','-')			% Solid grid lines
   % set (gca, 'XTickLabel', Xlab)			% Use Xlab for tick labels
   % set (gca, 'XTick', 0*srate:1.0*srate:PLOT_TIME*srate) 
   set (gca, 'Box', 'on')				
   set (gca, 'Ytick', 0:spacing_var:chans*spacing_var)  % ytick spacing on channels
   set (gca, 'TickLength', [0.02 0.02])
   title(titleval)					        % title is titleval
   axis([xstart xstart+PLOT_TIME 0 (chans+1)*spacing_var]);
   set (gca, 'YTickLabels', flipud(channames),'FontSize',14) % write channel names
                                            % Matlab-5 requires YTickLabels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the selected EEG data epoch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   xx = xstart:1/srate:xstart+PLOT_TIME-0.9/srate; % x-values
   xx = xx(:,1:size(data,2));
for i = 1:chans			
   F = data(chans-i+1,:);
   F = F - mean(F) + i*spacing_var;  % add offset to y-values
   plot (xx,F,'clipping','off','Color',linecolor); % channel plot with x-values
end 
if xstart<0 && xstart+PLOT_TIME > 0
   linetime = round(-xstart/srate);
   line ([linetime linetime],[1e10,-1e10]);
end

if vertmark ~= 0,
   for mark = vertmark,
       linetime = xstart + (mark-1)/srate;
       line ([linetime linetime],[-1e10,1e10],'Color','y','LineStyle',':');
   end
end
