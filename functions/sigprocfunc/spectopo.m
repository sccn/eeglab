% spectopo() - Plot the mean log spectrum of a set of data epochs at
%              all channels as a bundle of traces. At specified frequencies,
%              plot the relative topographic distribution of power.
%              Uses Matlab psd() from signal processing toolbox.
% Usage:
%        >> [spectra,freqs] = spectopo(data, frames, srate, 'key1', 'val1' ...
%                                                           'key2', 'val2' ...);
%
% Inputs:
%       data   = 2D (nchans,frames*epochs); % can be single-epoch
%                 or 3D (nbchan,frames,epochs)
%       frames = frames per epoch {0 -> data length}
%       srate  = sampling rate per channel (Hz)
%
% Optional inputs:
%   'freq'     = [float vector (Hz)] vector of frequencies for topoplot() scalp maps
%                of power at all channels, or single frequency to plot component 
%                contributions at a single channel ('plotchan').
%   'chanlocs' = electrode locations file (format: >> topoplot example)
%   'limits'   = axis limits [xmin xmax ymin ymax cmin cmax]
%                To use data limits, omit final values or use nan's
%                i.e. [-100 900 nan nan -10 10], [-100 900]
%                Note that default color limits are symmetric around 0 and are
%                different for each head {defaults: all nans}
%   'title'    = [quoted string] plot title {default: none}
%   'freqfac'  = [int power of 2] approximate frequencies/Hz to compute {default: 4}
%   'percent'  = downsampling factor or approximate percentage of the data to
%                keep while computing spectra. Downsampling can be used to speed up
%                the computation. From 0 to 100 {default: 100}.
%   'maxfreq'  = [float] maximum frequency to plot. Overwrite limits x axis.
%   'reref'    = ['averef'|'off'] convert input data to average reference 
%                Default is 'off'. 
%
% Plot component contributions:
%   'weights'  = ICA unmixing matrix. 'freq' must contain a single frequency.
%                ICA maps of the N (='nicamaps') components that account for the most
%                power at the selected frequency ('freq') are plotted along with
%                the spectra of the selected channel ('plotchan') and components
%                ('icacomps').
%   'plotchan' = [integer] channel at which to compute independent conmponent
%                contributions at the selected frequency ('freq'). {[]=channel with
%                higest power at 'freq').If 0, plot RMS power at all channels. 
%   'nicamaps' = [integer] number of ICA component maps to plot (Default 4).
%   'icacomps' = [integer array] indices of ICA component spectra to plot ([]=all).
%   'icamode'  = ['normal'|'sub'] in 'sub' mode, instead of computing the spectrum of
%                individual ICA components, the function computes the spectrum of
%                the data minus their contribution { default: 'normal' }
%   'icamaps'  = [integer array] force plotting of selected ica compoment maps ([]=the
%                'nicamaps' largest).
%   'memory'   = ['low'|'high'] setting to low will use less memory for ICA component 
%                computing, but will be longer.
%
% Topoplot options:
%    opther 'key','val' options are propagated to topoplot() for map display
%    (see help topoplot())
%
% Outputs:
%        spectra = (nchans,nfreqs) power spectra (average over epochs) in dB
%        freqs   = frequencies of spectra (Hz)
%
% Notes: The old function call is still function for backward compatibility
%        >> [spectra,freqs] = spectopo(data, frames, srate, headfreqs, ...
%                               chanlocs, limits, titl, freqfac, percent);
%
% Authors: Scott Makeig, Arnaud Delorme & Marissa Westerfield, 
%          SCCN/INC/UCSD, La Jolla, 3/01 
%
% See also: timtopo(), envtopo(), tftopo(), topoplot()

% Copyright (C) 3/01 Scott Makeig & Arnaud Delorme & Marissa Westerfield, SCCN/INC/UCSD, 
% scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: not supported by cvs2svn $
% Revision 1.22  2002/07/29 22:22:27  arno
% update messages
%
% Revision 1.21  2002/07/26 02:11:02  arno
% debugging max limits
%
% Revision 1.20  2002/07/26 01:15:49  arno
% adding percentage of variance
%
% Revision 1.19  2002/07/25 20:58:33  luca
% divide sum of epoch spectra by length(epoch_subset).
%
% Revision 1.18  2002/07/25 20:44:32  luca
% *** empty log message ***
%
% Revision 1.17  2002/07/25 00:36:32  arno
% debugging
%
% Revision 1.16  2002/07/25 00:18:23  scott
% plot diag comp map lines in same colors as comp spec traces. -sm & ad
%
% Revision 1.15  2002/07/24 17:47:55  arno
% making it nicer for components
%
% Revision 1.14  2002/07/20 19:21:21  arno
% debugging
%
% Revision 1.13  2002/07/20 01:55:54  arno
% add ignore extras
%
% Revision 1.12  2002/07/20 01:51:14  arno
% back to old way of computing average of power
%
% Revision 1.11  2002/07/20 01:44:12  arno
% checking parameter
%
% Revision 1.10  2002/07/20 01:41:07  arno
% percent update
%
% Revision 1.9  2002/07/20 01:35:14  arno
% percent from 0 to 100
%
% Revision 1.8  2002/07/20 01:17:10  arno
% new version with component plotting options
%
% Revision 1.7  2002/07/18 16:00:46  arno
% adding option for not plotting channels
%
% Revision 1.6  2002/07/07 22:44:49  scott
% *** empty log message ***
%
% Revision 1.5  2002/07/07 22:42:17  scott
% help msg -sm
%
% Revision 1.4  2002/07/07 22:38:27  scott
% adding 'reref','averef' option -sm
%
% Revision 1.3  2002/04/21 00:38:25  scott
% 'Selecting randomly' -> 'Randomly selecting' -sm
%
% Revision 1.2  2002/04/18 18:19:28  arno
% adding 3D option
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 3-20-01 added limits arg -sm
% 01-25-02 reformated help & license -ad 
% 02-15-02 scaling by epoch number line 108 - ad, sm & lf
% 03-15-02 add all topoplot options -ad
% 03-18-02 downsampling factor to speed up computation -ad
% 03-27-02 downsampling factor exact calculation -ad
% 04-03-02 added axcopy -sm

% Uses: MATLAB psd(), changeunits(), topoplot(), textsc()

function [eegspecdB,freqs]=spectopo(data,frames,srate,varargin) 
	%headfreqs,chanlocs,limits,titl,freqfac, percent, varargin)

LOPLOTHZ = 1;  % low  Hz to plot
FREQFAC  = 4;  % approximate frequencies/Hz (default)

if nargin<3
   help spectopo
   return
end
if nargin <= 3 | isstr(varargin{1})
	% 'key' 'val' sequency
	fieldlist = { 'freq'          'real'     []                        [] ;
				  'chanlocs'      ''         []                        [] ;
				  'maxfreq'       'real'     [1 srate/2]               [] ;
				  'memory'        'string'   {'low' 'high'}           'high' ;
				  'title'         'string'   []                       '';
				  'limits'        'real'     []                       [nan nan nan nan nan nan];
				  'freqfac'       'integer'  []                        FREQFAC;
				  'percent'       'real'     [0 100]                  100 ;
				  'reref'         'string'   { 'averef' 'no' }         'no' ;
				  'icamode'       'string'   { 'normal' 'sub' }        'normal' ;
				  'weights'       'real'     []                       [] ;
				  'plotchan'      'integer'  [1:size(data,1)]         [] ;
				  'nicamaps'      'integer'  []                       4 ;
				  'icacomps'      'integer'  []                       [] ;
				  'icamaps'       'integer'  []                       [] };
	
	[g varargin] = finputcheck( varargin, fieldlist, 'spectopo', 'ignore');
	if isstr(g), error(g); end;
	if ~isempty(g.maxfreq), g.limits(2) = g.maxfreq; end;
	if ~isempty(g.weights)
		if isempty(g.freq) | length(g.freq) > 2
			error('spectopo: for computing component contribution, one must specify a (single) frequency');
		end;
	end;
else
	if nargin > 3,    g.freq = varargin{1};
	else              g.freq = [];
	end;
	if nargin > 4,	  g.chanlocs = varargin{2};
	else              g.chanlocs = [];
	end;
	if nargin > 5,    g.limits = varargin{3};
	else              g.limits = [nan nan nan nan nan nan];
	end;
	if nargin > 6,    g.title = varargin{4};
	else              g.title = '';
	end;
	if nargin > 7,    g.freqfac = varargin{5};
	else              g.freqfac = FREQFAC;
	end;
	if nargin > 8,    g.percent = varargin{6};
	else              g.percent = 100;
	end;
	if nargin > 10,    g.reref = 'averef';
	else               g.reref = 'no';
	end;
	g.weights = [];
	g.icamaps = [];
end;
if g.percent > 1
	g.percent = g.percent/100; % make it from 0 to 1
end;
if ~isempty(g.freq) & isempty(g.chanlocs)
	error('spectopo: need channel location file');
end;

data = reshape(data, size(data,1), size(data,2)*size(data,3));
if frames == 0
  frames = size(data,2); % assume one epoch
end

if ~isempty(g.freq) & min(g.freq)<0
   fprintf('spectopo(): freqs must be >=0 Hz\n');
   return
end
epochs = round(size(data,2)/frames);
if frames*epochs ~= size(data,2)
   error('Spectopo: non-integer number of epochs');
end
if ~isempty(g.weights)
	g.icawinv = pinv(g.weights); % maps
	if ~isempty(g.icacomps)
		g.weights = g.weights(g.icacomps, :);
		g.icawinv = g.icawinv(:,g.icacomps);
	else 
		g.icacomps = [1:size(g.weights,1)];
	end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute channel spectra using psd()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epoch_subset = 1:epochs;
if g.percent ~= 1 & epochs > 1
    nb = round( g.percent*epochs);
    epoch_subset = zeros(1,epochs);
    while nb>0
        index = ceil(rand*epochs);
        if ~epoch_subset(index)
            epoch_subset(index) = 1;
            nb = nb-1;
        end;
    end;        
    epoch_subset = find(epoch_subset == 1);
    fprintf('Randomly selecting %d of %d data epochs for analysis...\n', length(epoch_subset),epochs);
end;
if isempty(g.weights)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% compute data spectrum
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	fprintf('Computing spectra: ')
	[eegspecdB freqs] = spectcomp( data, frames, srate, epoch_subset, g);
	fprintf('\n');
else
	% compute data spectrum
	if isempty(g.plotchan) | g.plotchan == 0
		fprintf('Computing spectra: ')
		[eegspecdB freqs] = spectcomp( data, frames, srate, epoch_subset, g);
		fprintf('\n');
	else
		fprintf('Computing spectra at specified channel: ')
		[eegspecdB freqs] = spectcomp( data(g.plotchan,:), frames, srate, epoch_subset, g);
		fprintf('\n');
	end;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% select channel and spectrum
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isempty(g.plotchan) % find channel of minimum power
		[tmp indexfreq] = min(abs(g.freq-freqs));
		[tmp g.plotchan] = min(eegspecdB(:,indexfreq));
		fprintf('Maximum power found at channel %d\n', g.plotchan);
		eegspecdBtoplot = eegspecdB(g.plotchan, :);		
	elseif g.plotchan == 0
		fprintf('Averaging power at all channels\n');
		eegspecdBtoplot = mean(eegspecdB, 1);
	else 
		eegspecdBtoplot = eegspecdB;
	end;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% compute component spectra
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	fprintf('Computing component spectra: ')
	if g.plotchan ~= 0
		newweights = diag(g.icawinv(g.plotchan,:)) * g.weights;
	else
		newweights = diag(sum(g.icawinv,1)) * g.weights;
	end;
	if strcmp(g.memory, 'high') & strcmp(g.icamode, 'normal')
		[compeegspecdB freqs] = spectcomp( newweights*data, frames, srate, epoch_subset, g);
	else % in case out of memory error, multiply conmponent sequencially
		[compeegspecdB freqs] = spectcomp( data, frames, srate, epoch_subset, g, newweights);
	end;
	fprintf('\n');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% select components to plot
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isempty(g.icamaps)
		[tmp indexfreq] = min(abs(g.freq-freqs));
		g.icafreqsval   = compeegspecdB(:, indexfreq);
		[g.icafreqsval g.icamaps] = sort(g.icafreqsval);
		if strcmp(g.icamode, 'normal')
			g.icamaps = g.icamaps(end:-1:1);
			g.icafreqsval = g.icafreqsval(end:-1:1);
		end;
		if g.nicamaps < length(g.icamaps), g.icamaps = g.icamaps(1:g.nicamaps); end;
	else 
		[tmp indexfreq] = min(abs(g.freq-freqs));
		g.icafreqsval   = compeegspecdB(g.icamaps, indexfreq);
	end;
end;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute axis and caxis g.limits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(g.limits)<1 | isnan(g.limits(1))
   g.limits(1) = LOPLOTHZ;
end

if ~isempty(g.freq)
	if length(g.limits)<2 | isnan(g.limits(2))
		maxheadfreq = max(g.freq);
		if rem(maxheadfreq,5) ~= 0
			g.limits(2) = 5*ceil(maxheadfreq/5);
		else
			g.limits(2) = maxheadfreq*1.1;
		end
	end
	
	g.freq = sort(g.freq);          % Determine topoplot frequencies
	freqidx = zeros(1,length(g.freq)); % Do not interpolate between freqs
	for f=1:length(g.freq)
		[tmp fi] = min(abs(freqs-g.freq(f)));
		freqidx(f)=fi;
	end
else 
	g.limits(2) = 50;
end;

[tmp maxfreqidx] = min(abs(g.limits(2)-freqs)); % adjust max frequency
[tmp minfreqidx] = min(abs(g.limits(1)-freqs)); % adjust min frequency

if length(g.limits)<3|isnan(g.limits(3))
	reallimits(1) = min(min(eegspecdB(:,minfreqidx:maxfreqidx)));
else 
	reallimits(1) = g.limits(3);
end
if length(g.limits)<4|isnan(g.limits(4))
	reallimits(2) = max(max(eegspecdB(:,minfreqidx:maxfreqidx)));
	dBrange = reallimits(2)-reallimits(1);   % expand range a bit beyond data g.limits
	reallimits(1) = reallimits(1)-dBrange/7;
	reallimits(2) = reallimits(2)+dBrange/7;
else 
	reallimits(2) = g.limits(4);
end

if length(g.limits)<5 % default caxis plotting g.limits
  g.limits(5) = nan;
end
if length(g.limits)<6 
  g.limits(6) = nan;
end

if isnan(g.limits(5))+isnan(g.limits(6)) == 1
   fprintf('spectopo(): limits 5 and 6 must either be given or nan\n');
   return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot spectrum of each channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(g.freq)
	specaxes = sbplot(3,4,[5 12]); 
end;

if isempty(g.weights)
	pl=plot(freqs(1:maxfreqidx),eegspecdB(:,1:maxfreqidx)');
else 
	pl=plot(freqs(1:maxfreqidx),eegspecdBtoplot(:,1:maxfreqidx)');
end;
set(pl,'LineWidth',2);
set(gca,'TickLength',[0.02 0.02]);
axis([freqs(minfreqidx) freqs(maxfreqidx) reallimits(1) reallimits(2)]);
xl=xlabel('Frequency (Hz)');
set(xl,'fontsize',16);
yl=ylabel('Rel. Power (dB)');
set(yl,'fontsize',16);
set(gca,'fontsize',16)
box off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   plot component contribution   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colrs = {'r','b','g','m','c'}; % component spectra trace colors
if ~isempty(g.weights)
	set(pl, 'linewidth', 2, 'color', 'k');
	hold on;
	for f=1:length(g.icamaps)
		colr = colrs{mod((f-1),5)+1};
		pl2=plot(freqs(1:maxfreqidx),compeegspecdB(g.icamaps(f),1:maxfreqidx)',colr);
	end
	othercomps = setdiff(1:size(compeegspecdB,1), g.icamaps);
	if ~isempty(othercomps)
		pl2=plot(freqs(1:maxfreqidx),compeegspecdB(othercomps,1:maxfreqidx)');
	end;
	if length(g.limits)<3|isnan(g.limits(3))
		newaxis = axis;
		newaxis(3) = min(newaxis(3), min(min(compeegspecdB(:,1:maxfreqidx))));
		newaxis(4) = max(newaxis(4), max(max(compeegspecdB(:,1:maxfreqidx))));
		axis(newaxis);
	end;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% indicate component contribution %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if g.plotchan ~= 0
		maxdatadb = max(eegspecdBtoplot(:,freqidx(1)));
		[tmp indexfreq] = min(abs(g.freq-freqs));
		for index = 1:length(g.icacomps)
			if strcmp(g.icamode, 'normal')
				resvar  = 100*exp(-(maxdatadb-compeegspecdB(index, indexfreq))/10*log(10));
				fprintf('Component %d percentage relative variance:%6.2f\n', g.icacomps(index), resvar);
			else
				resvar  = 100 - 100*exp(-(maxdatadb-compeegspecdB(index, indexfreq))/10*log(10));
				fprintf('Component %d percentage variance accounted for:%6.2f\n', g.icacomps(index), resvar);
			end;
		end;
	else
		fprintf(['Warning: can not compute component percentage contribution\n' ...
				 'when considering global power (i.e. when using ''plotchan'', 0)\n']);
	end;
end;

if ~isempty(g.freq)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Plot lines through channel trace bundle at each headfreq
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isempty(g.weights)
		for f=1:length(g.freq)
			hold on; 
			plot([freqs(freqidx(f)) freqs(freqidx(f))], ...
				 [min(eegspecdB(:,freqidx(f))) max(eegspecdB(:,freqidx(f)))],...
				 'k','LineWidth',2.5);
		end;
	else
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Plot vertical line at comp analysis freq
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		mincompdB = min([min(eegspecdB(:,freqidx(1))) min(compeegspecdB(:,freqidx(1)))]);
		maxcompdB = max([max(eegspecdB(:,freqidx(1))) max(compeegspecdB(:,freqidx(1)))]);
		hold on;
		plot([freqs(freqidx(1)) freqs(freqidx(1))],[mincompdB maxcompdB],'k', 'LineWidth',2.5);
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Plot connecting lines using changeunits()
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	tmpmainpos = get(gca, 'position');
	headax = zeros(1,length(g.freq));
	for f=1:length(g.freq)+length(g.icamaps)
		headax(f) = sbplot(3,length(g.freq)+length(g.icamaps),f);
		axis([-1 1 -1 1]);
		
		%axis x coords and use
		tmppos = get(headax(f), 'position');
		allaxcoords(f) = tmppos(1);
		allaxuse(f)    = 0;
	end
	large = sbplot(1,1,1);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% compute relative positions on plot
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if  ~isempty(g.weights)
		freqnormpos = tmpmainpos(1) + tmpmainpos(3)*(freqs(freqidx(1))-g.limits(1))/(g.limits(2)-g.limits(1));
		for index = 1:length(g.icamaps)+1
			[realpos(index) allaxuse] = closestplot( freqnormpos, allaxcoords, allaxuse );
		end;
	else 
		realpos = 1:length(g.freq); % indices giving order of plotting positions
	end;
	
	for f=1:length(g.freq)+length(g.icamaps)
		if f>length(g.freq) % special case of ica components
			from = changeunits([freqs(freqidx(1)),g.icafreqsval(f-1)],specaxes,large);
			%g.icafreqsval contains the sorted frequency values at the specified frequency
		else 
			from = changeunits([freqs(freqidx(f)),max(eegspecdB(:,freqidx(f)))],specaxes,large);
		end;
		pos = get(headax(realpos(f)),'position');
		to = changeunits([0,0],headax(realpos(f)),large)+[0 -min(pos(3:4))/1.7];
		hold on;
    if f<=length(g.freq)
      colr = 'k';
    else
      colr = colrs{mod((f-2),5)+1};
    end
		li(realpos(f)) = plot([from(1) to(1)],[from(2) to(2)],colr,'LineWidth',2);
		axis([0 1 0 1]);
		axis off;
	end;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Plot heads using topoplot()
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	fprintf('Plotting scalp distributions: ')
	for f=1:length(g.freq)
		axes(headax(realpos(f)));
		topodata = eegspecdB(:,freqidx(f))-mean(eegspecdB(:,freqidx(f)));
		if isnan(g.limits(5)),     maplimits = 'absmax';
		else                       maplimits = [g.limits(5) g.limits(6)];
		end;
		if length(topodata) == 1
			if ~isempty(varargin)
				topoplot(g.plotchan,g.chanlocs,'electrodes','off', ...
						 'style', 'blank', 'emarkersize1chan', 10, varargin{:});
			else
				topoplot(g.plotchan,g.chanlocs,'electrodes','off', ...
						 'style', 'blank', 'emarkersize1chan', 10);
			end
			if isstruct(g.chanlocs)
				tl=title(g.chanlocs(g.plotchan).labels);
			else
				tl=title([ 'c' int2str(g.plotchan)]);
			end;
		else	
			if ~isempty(varargin)
				topoplot(topodata,g.chanlocs,'maplimits',maplimits, varargin{:}); 
			else
				topoplot(topodata,g.chanlocs,'maplimits',maplimits); 
			end
			if f<length(g.freq)
				tl=title([num2str(freqs(freqidx(f)), '%3.1f')]);
			else
				tl=title([num2str(freqs(freqidx(f)), '%3.1f') ' Hz']);
			end
		end;
		set(tl,'fontsize',16);
		axis square;
		drawnow
		fprintf('.');
	end;
	fprintf('\n');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Plot independant components
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if ~isempty(g.weights)
		if realpos(1) == max(realpos), plotcolbar(g); end;
		% use headaxe from 2 to end (reserved earlier)
		set(li(realpos(1)), 'linewidth', 3); % make the line with the scalp topoplot thicker than others
		for index = 1:length(g.icamaps)
			axes(headax(realpos(index+1)));						
			compnum = g.icamaps(index);
			topoplot(g.icawinv(:,compnum),g.chanlocs,varargin{:}); 
			tl=title(int2str(compnum));
			set(tl,'fontsize',16);
			axis square;
			drawnow
			if realpos(index+1) == max(realpos), plotcolbar(g); end;
		end;
	else 
		plotcolbar(g);
	end;

end;

%%%%%%%%%%%%%%%%
% Draw title
%%%%%%%%%%%%%%%%
if ~isempty(g.title)
  tl = textsc(g.title,'title');
  set(tl,'fontsize',15)
end

%%%%%%%%%%%%%%%%
% Turn on axcopy
%%%%%%%%%%%%%%%%
axcopy

%%%%%%%%%%%%%%%%
% Plot color bar
%%%%%%%%%%%%%%%%
function plotcolbar(g)
	cb=cbar;
	pos = get(cb,'position');
	set(cb,'position',[pos(1) pos(2) 0.03 pos(4)]);
	set(cb,'fontsize',12);
	try
		if isnan(g.limits(5))
			ticks = get(cb,'ytick');
			[tmp zi] = find(ticks == 0);
			ticks = [ticks(1) ticks(zi) ticks(end)];
			set(cb,'ytick',ticks);
			set(cb,'yticklabel',{'-','0','+'});
		end
	catch, end; % in a single channel is plotted

%%%%%%%%%%%%%%%%%%%%%%%
% function closest plot
%%%%%%%%%%%%%%%%%%%%%%%
function [index, usedplots] = closestplot(xpos, xcentercoords, usedplots);
	notused = find(usedplots == 0);
	xcentercoords = xcentercoords(notused);
	[tmp index] = min(abs(xcentercoords-xpos));
	index = notused(index);
	usedplots(index) = 1;
	return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function computing spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eegspecdB, freqs] = spectcomp( data, frames, srate, epoch_subset, g, newweights);
	if exist('newweights') == 1
		nchans = size(newweights,1);
	else 
		nchans = size(data,1);		
	end;
	fftlength = 2^round(log(srate)/log(2))*g.freqfac;
	for c=1:nchans % scan channels or components
		if exist('newweights') == 1
			if strcmp(g.icamode, 'normal')
				tmpdata = newweights(c,:)*data; % component activity
			else % data - component contribution
				tmpdata = data(g.plotchan,:) - newweights(c,:)*data;
			end;
		else
			tmpdata = data(c,:); % channel activity
		end;
		for e=epoch_subset
			if strcmp(g.reref, 'averef')
				[tmpspec,freqs] = psd(averef(matsel(tmpdata,frames,0,1,e)),fftlength,...
									  srate,fftlength/2,fftlength/8);
			else
				[tmpspec,freqs] = psd(matsel(tmpdata,frames,0,1,e),fftlength,...
									  srate,fftlength/2,fftlength/8);
			end
			if c==1 & e==epoch_subset(1)
				eegspec = zeros(nchans,length(freqs));
			end
			eegspec(c,:) = eegspec(c,:) + tmpspec';
		end
		fprintf('.')
	end
	eegspecdB = 10*log10(eegspec/length(epoch_subset)); % convert power to dB
	return;
