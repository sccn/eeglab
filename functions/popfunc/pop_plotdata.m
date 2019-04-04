% pop_plotdata() - Plot average of EEG channels or independent components in
%                  a rectangular array. Else, (over)plot single trials.
% Usage:
%   >> avg = pop_plotdata(EEG, typeplot, indices, trials, title, singletrials, ydir, ylimits);
%
% Inputs:
%   EEG        - Input dataset 
%   typeplot   - Type data to plot (1=channels, 0=components) {Default:1}
%   indices    - Array of channels (or component) indices to plot 
%                     {Default: all}
%   trials     - Array of trial indices. sum specific trials in the average 
%                     {Default: all}
%   title      - Plot title. {Default: []}.
%   singletrials - [0|1], Plot average or overplot single trials 
%                      0 plot average, 1 plot single trials {Default: 0}
%   ydir       - [1|-1] y-axis polarity (pos-up = 1; neg-up = -1) 
%                {def command line-> pos-up; def GUI-> neg-up}
%   ylimits    - [ymin ymax] plotting limits {default [0 0] -> data range}
%
% Outputs:
%   avg        - [matrix] Data average
% 
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: plotdata(), eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% 01-25-02 reformated help & license -ad 
% 03-08-02 add eeglab options -ad
% 03-18-02 added title -ad & sm
% 03-30-02 added single trial capacities -ad

function [sigtmp, com] = pop_plotdata(EEG, typeplot, indices, trials, plottitle, singletrials, ydir,ylimits);

ylimits = [0 0];% default ylimits = data range

% warning signal must be in a 3D form before averaging

sigtmp = [];

com = '';
if nargin < 1
   help pop_plotdata;
   return;
end
if nargin < 2
	typeplot = 1; % 1=signal; 0=component
end
if exist('plottitle') ~= 1
	plottitle = '';
end

if nargin <3
	if typeplot % plot signal channels
		result = inputdlg2({  'Channel number(s):'  ...
                                      'Plot title:'  ...
                                      'Vertical limits ([0 0]-> data range):' ...
                                   }, ...
                                  'Channel ERPs in rect. array -- pop_plotdata()', 1, ...
                                   {['1:' int2str(EEG.nbchan)] [fastif(isempty(EEG.setname), '',[EEG.setname ' ERP'])] ['0 0'] }, ...
                                   'pop_plotdata' );
	else % plot components
		result = inputdlg2({  'Component number(s):' ...
                                      'Plot title:' ...
                                      'Vertical limits ([0 0]-> data range):' ...
                                   }, ...
                                   'Component ERPs in rect. array -- pop_plotdata()', 1, ...
                                   {   ['1:' int2str(size(EEG.icawinv,2))] [fastif(isempty(EEG.setname), '',[EEG.setname ' ERP'])] ['0 0'] }, ...
                                   'pop_plotdata' );
	end;		
	if length(result) == 0 return; end

	indices   	 = eval( [ '[' result{1} ']' ] );

	plottitle    = result{2};
 
	singletrials = 0;

	ylimits   	 = eval( [ '[' result{3} ']' ] );
        if length(ylimits) ~= 2
            ylimits = [0 0]; % use default if 2 values not given
        end
end;	
if ~(exist('trials') == 1)
	trials = 1:EEG.trials;
end;	
if exist('plottitle') ~= 1
    plottitle = '';
end;    
if exist('singletrials') ~= 1
    singletrials = 0;
end;    

if EEG.trials > 1 && singletrials == 0
    fprintf('Selecting trials and components...\n');
	if typeplot == 1
	   sigtmp = nan_mean(EEG.data(indices,:,trials),3);
 	else
	   if isempty(EEG.icasphere)
	      error('no ICA data for this set, first run ICA');
	   end;   
       tmpdata = eeg_getdatact(EEG, 'component', indices, 'trialindices', trials);
	   fprintf('Averaging...\n');
	   sigtmp = nan_mean(tmpdata,3);
	end
else
	if typeplot == 1
	   sigtmp = EEG.data(indices,:,trials);
	else
	   if isempty(EEG.icasphere)
	      error('no ICA data for this set, first run ICA');
	   end;  
       sigtmp = eeg_getdatact(EEG, 'component', indices, 'trialindices', trials);
	end
end
figure;
try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end
if exist('YDIR') ~= 1
     ydir = 1;
else
     ydir = YDIR;
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%% make the plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
sigtmp = reshape( sigtmp, size(sigtmp,1),  size(sigtmp,2)*size(sigtmp,3));
ymin = ylimits(1);
ymax = ylimits(2);

if ischar(ymin)
   ymin = str2num(ymin);
   ymax = str2num(ymax);
end

% com = sprintf('pop_plotdata(%s, %d, %s, [1:%d], ''%s'', %d, %d, [%f %f]);', ...
% inputname(1), typeplot, vararg2str(indices), EEG.trials, plottitle, singletrials,ydir,ymin,ymax);
% fprintf([com '\n']);
if ~isempty(EEG.chanlocs) && typeplot == 1
    tmpchanlocs = EEG.chanlocs;
    chanlabels = strvcat({ tmpchanlocs(indices).labels });
else
    chanlabels = num2str(indices(:));
end
plottopo( sigtmp, 'frames', EEG.pnts, 'limits', [EEG.xmin*1000 EEG.xmax*1000 ymin ymax], 'title', plottitle, 'chans', 1:size(sigtmp,1), 'ydir', ydir, 'channames', chanlabels);
%plotdata(sigtmp, EEG.pnts, [EEG.xmin*1000 EEG.xmax*1000 ymin ymax], plottitle, indices,0,0,ydir);

%
%%%%%%%%%%%%%%%%%%%%%%%%%% add figure title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if typeplot == 1
	set(gcf, 'name', 'Plot > Channel ERPs > In rect. array -- plotdata()');
else
	set(gcf, 'name', 'Plot > Component ERPs > In rect. array -- plotdata()');
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%% set y-axis direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if exist('ydir') ~= 1
  if exist('YDIR') ~= 1
    ydir = 1;  % default positive up
  else
    ydir = YDIR; % icadefs.m system-wide default
  end
end

if ydir==1
   set(gca,'ydir','normal');
else % ydir == -1
   set(gca,'ydir','reverse');
end

switch nargin
	case {0, 1, 2, 3}, com = sprintf('pop_plotdata(EEG, %d, %s, [1:%d], ''%s'', %d, %d, [%g %g]);', typeplot, vararg2str(indices), EEG.trials, plottitle, singletrials,ydir,ymin,ymax);
	case 4, com = sprintf('pop_plotdata(EEG, %d, %s, %s, ''%s'', %d, %d, [%g %g]);', typeplot, vararg2str(indices), vararg2str(trials), plottitle, singletrials,ydir,ymin,ymax);
end

fprintf([com '\n']);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = nan_mean(in, dim)
    tmpin = in;
    tmpin(find(isnan(in(:)))) = 0;
    out = sum(tmpin, dim) ./ sum(~isnan(in),dim);
