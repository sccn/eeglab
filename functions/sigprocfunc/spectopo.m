% spectopo() - Plot the mean log spectrum of a set of data epochs at all channels 
%              as a bundle of traces. At specified frequencies, plot the relative 
%              topographic distribution of power. If available, uses pwelch() from 
%              the Matlab signal processing toolbox, else the EEGLAB spec() function.
%              Plots the mean spectrum for all of the supplied data, not just
%              the pre-stimulus baseline.
% Usage:
%              >> spectopo(data, frames, srate);
%              >> [spectra,freqs,speccomp,contrib,specstd] = ...
%                    spectopo(data, frames, srate, 'key1','val1', 'key2','val2' ...);
% Inputs:
%       data   = If 2-D (nchans,time_points); % may be a continuous single epoch,
%                else a set of concatenated data epochs, else a 3-D set of data 
%                epochs (nchans,frames,epochs)
%       frames = frames per epoch {default|0 -> data length}
%       srate  = sampling rate per channel (Hz)
%
% Optional 'keyword',[argument] input pairs:
%   'freq'     = [float vector (Hz)] vector of frequencies at which to plot power 
%                scalp maps, or else a single frequency at which to plot component 
%                contributions at a single channel (see also 'plotchan').
%   'chanlocs' = [electrode locations filename or EEG.chanlocs structure]. 
%                    For format, see >> topoplot example
%   'limits'   = [xmin xmax ymin ymax cmin cmax] axis limits. Sets x, y, and color 
%                axis limits. May omit final values or use NaNs.
%                   Ex: [0 60 NaN NaN -10 10], [0 60], ...
%                Default color limits are symmetric around 0 and are different 
%                for each scalp map {default|all NaN's: from the data limits}
%   'title'    = [quoted string] plot title {default: none}
%   'freqfac'  = [integer] ntimes to oversample -> frequency resolution {default: 2}
%   'nfft'     = [integer] length to zero-pad data to. Overwrites 'freqfac' above.
%   'winsize'  = [integer] window size in data points {default: from data}
%   'overlap'  = [integer] window overlap in data points {default: 0}
%   'percent'  = [float 0 to 100] percent of the data to sample for computing the 
%                spectra. Values < 100 speed up the computation. {default: 100}.
%   'freqrange' = [min max] frequency range to plot. Changes x-axis limits {default: 
%                1 Hz for the min and Nyquist (srate/2) for the max. If specified 
%                power distribution maps are plotted, the highest mapped frequency 
%                determines the max freq}.
%   'reref'    = ['averef'|'off'] convert data to average reference {default: 'off'}
%   'mapnorm'  = [float vector] If 'data' contain the activity of an independant 
%                component, this parameter should contain its scalp map. In this case
%                the spectrum amplitude will be scaled to component RMS scalp power.
%                Useful for comparing component strengths {default: none}
%   'boundaries' = data point indices of discontinuities in the signal {default: none}
%   'plot'     = ['on'|'off'] 'off' -> disable plotting {default: 'on'}
%   'rmdc'     = ['on'|'off'] 'on' -> remove DC {default: 'off'}  
%   'plotmean' = ['on'|'off'] 'on' -> plot the mean channel spectrum {default: 'off'}  
%   'plotchans' = [integer array] plot only specific channels {default: all}
%
% Optionally plot component contributions:
%   'weights'  = ICA unmixing matrix. Here, 'freq' (above) must be a single frequency.
%                ICA maps of the N ('nicamaps') components that account for the most
%                power at the selected frequency ('freq') are plotted along with
%                the spectra of the selected channel ('plotchan') and components
%                ('icacomps').
%   'plotchan' = [integer] channel at which to compute independent conmponent
%                contributions at the selected frequency ('freq'). If 0, plot RMS 
%                power at all channels. {defatul|[] -> channel with highest power 
%                at specified 'freq' (above)). Do not confuse with
%                'plotchans' which select channels for plotting.
%   'mapchans' = [int vector] channels to plot in topoplots {default: all}
%   'mapframes'= [int vector] frames to plot {default: all}
%   'nicamaps' = [integer] number of ICA component maps to plot {default: 4}.
%   'icacomps' = [integer array] indices of ICA component spectra to plot ([] -> all).
%   'icamode'  = ['normal'|'sub'] in 'sub' mode, instead of computing the spectra of
%                individual ICA components, the function computes the spectrum of
%                the data minus their contributions {default: 'normal'}
%   'icamaps'  = [integer array] force plotting of selected ICA compoment maps 
%                {default: [] = the 'nicamaps' largest contributing components}.
%   'icawinv'  = [float array] inverse component weight or mixing matrix. Normally,
%                this is computed by inverting the ICA unmixing matrix 'weights' (above).
%                However, if any components were removed from the supplied 'weights'mapchans
%                then the component maps will not be correctly drawn and the 'icawinv'
%                matrix should be supplied here {default: from component 'weights'}
%   'memory'   = ['low'|'high'] a 'low' setting will use less memory for computing 
%                component activities, will take longer {default: 'high'}
%
% Replotting options:
%   'specdata' = [freq x chan array ] spectral data
%   'freqdata' = [freq] array of frequencies
% 
% Topoplot options:
%    other 'key','val' options are propagated to topoplot() for map display
%                (See >> help topoplot)
%
% Outputs:
%        spectra  = (nchans,nfreqs) power spectra (mean power over epochs), in dB
%        freqs    = frequencies of spectra (Hz)
%        speccomp = component spectra (ncomps,nfreqs). Warning, only the function  
%                   only computes the spectrum of the components given as input using 
%                   the 'icacomps' parameter. Other component spectrum are filled  
%                   with 0.
%        contrib  = contribution of each component (if 'icamode' is 'normal', relative
%                   variance, if 'icamode' is 'sub', percent variance accounted for)
%        specstd  = spectrum standard deviation in dB
%
% Notes: The original input format is still functional for backward compatibility.
%        psd() has been replaced by pwelch() (see Matlab note 24750 on their web site)
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

% 3-20-01 added limits arg -sm
% 01-25-02 reformated help & license -ad 
% 02-15-02 scaling by epoch number line 108 - ad, sm & lf
% 03-15-02 add all topoplot options -ad
% 03-18-02 downsampling factor to speed up computation -ad
% 03-27-02 downsampling factor exact calculation -ad
% 04-03-02 added axcopy -sm

% Uses: MATLAB pwelch(), changeunits(), topoplot(), textsc()

function [eegspecdB,freqs,compeegspecdB,resvar,specstd] = spectopo(data,frames,srate,varargin) 

% formerly: ... headfreqs,chanlocs,limits,titl,freqfac, percent, varargin)

icadefs;
LOPLOTHZ = 1;  % low  Hz to plot
FREQFAC  = 2;  % approximate frequencies/Hz (default)
allcolors = { [0 0.7500 0.7500] 
              [1 0 0] 
              [0 0.5000 0] 
              [0 0 1] 
              [0.2500 0.2500 0.2500] 
              [0.7500 0.7500 0] 
              [0.7500 0 0.7500] }; % colors from real plots                };

if nargin<3
   help spectopo
   return
end
if nargin <= 3 | isstr(varargin{1})
	% 'key' 'val' sequence
	fieldlist = { 'freq'          'real'     []                        [] ;
                  'specdata'      'real'     []                        [] ;
                  'freqdata'      'real'     []                        [] ;
				  'chanlocs'      ''         []                        [] ;
				  'freqrange'     'real'     [0 srate/2]               [] ;
				  'memory'        'string'   {'low','high'}           'high' ;
				  'plot'          'string'   {'on','off'}             'on' ;
				  'plotmean'      'string'   {'on','off'}             'off' ;
				  'title'         'string'   []                       '';
				  'limits'        'real'     []                       [nan nan nan nan nan nan];
				  'freqfac'       'integer'  []                        FREQFAC;
				  'percent'       'real'     [0 100]                  100 ;
				  'reref'         'string'   { 'averef','off','no' }  'off' ;
				  'boundaries'    'integer'  []                       [] ;
				  'nfft'          'integer'  [1 Inf]                  [] ;
				  'winsize'       'integer'  [1 Inf]                  [] ;
				  'overlap'       'integer'  [1 Inf]                  0 ;
				  'icamode'       'string'   { 'normal','sub' }        'normal' ;
				  'weights'       'real'     []                       [] ;
				  'mapnorm'       'real'     []                       [] ;
				  'plotchan'      'integer'  [1:size(data,1)]         [] ;
				  'plotchans'     'integer'  [1:size(data,1)]         [] ;
				  'nicamaps'      'integer'  []                       4 ;
				  'icawinv'       'real'     []                       [] ;
				  'icacomps'      'integer'  []                       [] ;
				  'icachansind'   'integer'  []                       [1:size(data,1)] ; % deprecated
				  'icamaps'       'integer'  []                       [] ;
                  'rmdc'           'string'   {'on','off'}          'off';
				  'mapchans'      'integer'  [1:size(data,1)]         [] 
                  'mapframes'     'integer'  [1:size(data,2)]         []};
	
	[g varargin] = finputcheck( varargin, fieldlist, 'spectopo', 'ignore');
	if isstr(g), error(g); end;
	if ~isempty(g.freqrange), g.limits(1:2) = g.freqrange; end;
	if ~isempty(g.weights)
		if isempty(g.freq) | length(g.freq) > 2
            if ~isempty(get(0,'currentfigure')) & strcmp(get(gcf, 'tag'), 'spectopo'), close(gcf); end;
         error('spectopo(): for computing component contribution, one must specify a (single) frequency');
		end;
	end;
else
    if ~isnumeric(data)
       error('spectopo(): Incorrect call format (see >> help spectopo).')
    end
    if ~isnumeric(frames) | round(frames) ~= frames
       error('spectopo(): Incorrect call format (see >> help spectopo).')
    end
    if ~isnumeric(srate)  % 3rd arg must be the sampling rate in Hz
       error('spectopo(): Incorrect call format (see >> help spectopo).')
    end

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
	else               g.reref = 'off';
	end;
	g.weights = [];
	g.icamaps = [];
end;
if g.percent > 1
	g.percent = g.percent/100; % make it from 0 to 1
end;
if ~isempty(g.freq) & isempty(g.chanlocs)
	error('spectopo(): needs channel location information');
end;
if isempty(g.weights) && ~isempty(g.plotchans)
    data = data(g.plotchans,:);
    if ~isempty(g.chanlocs)
        g.chanlocs = g.chanlocs(g.plotchans);
    end;
end;

if strcmpi(g.rmdc, 'on')
    data = data - repmat(mean(data,2), [ 1 size(data,2) 1]);
end
data = reshape(data, size(data,1), size(data,2)*size(data,3));

if frames == 0
  frames = size(data,2); % assume one epoch
end

%if ~isempty(g.plotchan) & g.plotchan == 0 & strcmpi(g.icamode, 'sub')
%    if ~isempty(get(0,'currentfigure')) & strcmp(get(gcf, 'tag'), 'spectopo'), close(gcf); end;
%    error('Cannot plot data component at all channels (option not implemented)');
%end;

if ~isempty(g.freq) & min(g.freq)<0
    if ~isempty(get(0,'currentfigure')) & strcmp(get(gcf, 'tag'), 'spectopo'), close(gcf); end;
   fprintf('spectopo(): freqs must be >=0 Hz\n');
   return
end

g.chanlocs2 = g.chanlocs;
if ~isempty(g.specdata)
    eegspecdB  = g.specdata;
    freqs      = g.freqdata;
else
    epochs = round(size(data,2)/frames);
    if frames*epochs ~= size(data,2)
       error('Spectopo: non-integer number of epochs');
    end
    if ~isempty(g.weights)
        ncompsori = size(g.weights,1);
        if isempty(g.icawinv)
            g.icawinv = pinv(g.weights); % maps
        end;
        if ~isempty(g.icacomps)
            g.weights = g.weights(g.icacomps, :);
            g.icawinv = g.icawinv(:,g.icacomps);
        else 
            g.icacomps = [1:size(g.weights,1)];
        end;
    end;
    compeegspecdB = [];
    resvar = NaN;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute channel spectra using pwelch()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    epoch_subset = ones(1,epochs);
    if g.percent ~= 1 & epochs == 1
        fprintf('Selecting the first %2.1f%% of data for analysis...\n', g.percent*100);
        frames = round(size(data,2)*g.percent);
        data = data(:, 1:frames);
        g.boundaries(find(g.boundaries > frames)) = [];
        if ~isempty(g.boundaries)
            g.boundaries(end+1) = frames;
        end;            
    end;
    if g.percent ~= 1 & epochs > 1
        epoch_subset = zeros(1,epochs);
        nb = ceil( g.percent*epochs);
        while nb>0
            index = ceil(rand*epochs);
            if ~epoch_subset(index)
                epoch_subset(index) = 1;
                nb = nb-1;
            end;
        end;        
        epoch_subset = find(epoch_subset == 1);
        fprintf('Randomly selecting %d of %d data epochs for analysis...\n', length(epoch_subset),epochs);
    else
        epoch_subset = find(epoch_subset == 1);
    end;
    if isempty(g.weights)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute data spectra
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('Computing spectra')
        [eegspecdB freqs specstd] = spectcomp( data, frames, srate, epoch_subset, g);
        if ~isempty(g.mapnorm) % normalize by component map RMS power (if data contain 1 component
            disp('Scaling spectrum by component RMS of scalp map power');
            eegspecdB       = sqrt(mean(g.mapnorm.^4)) * eegspecdB;
            % the idea is to take the RMS of the component activity (compact) projected at each channel
            % spec = sqrt( power(g.mapnorm(1)*compact).^2 + power(g.mapnorm(2)*compact).^2 + ...)
            % spec = sqrt( g.mapnorm(1)^4*power(compact).^2 + g.mapnorm(1)^4*power(compact).^2 + ...)
            % spec = sqrt( g.mapnorm(1)^4 + g.mapnorm(1)^4 + ... )*power(compact)
        end;

        tmpc = find(eegspecdB(:,1)); 			     % > 0 power chans
        tmpindices = find(eegspecdB(:,1) == 0);
        if ~isempty(tmpindices)
             zchans = int2str(tmpindices); % 0-power chans
        else zchans = [];
        end;
        if length(tmpc) ~= size(eegspecdB,1)
            fprintf('\nWarning: channels [%s] have 0 values, so will be omitted from the display', ...
                       zchans);
            eegspecdB = eegspecdB(tmpc,:);
            if ~isempty(specstd),  specstd = specstd(tmpc,:); end
            if ~isempty(g.chanlocs)
                g.chanlocs2 = g.chanlocs(tmpc);
            end
        end;
        eegspecdB = 10*log10(eegspecdB);
        specstd   = 10*log10(specstd);
        fprintf('\n');
    else
        % compute data spectrum
        if isempty(g.plotchan) | g.plotchan == 0
            fprintf('Computing spectra')
            [eegspecdB freqs specstd] = spectcomp( data, frames, srate, epoch_subset, g);	
            fprintf('\n'); % log below
        else
            fprintf('Computing spectra at specified channel')
            g.reref = 'off';
            [eegspecdB freqs specstd] = spectcomp( data(g.plotchan,:), frames, srate, epoch_subset, g);
            fprintf('\n'); % log below
        end;
        g.reref = 'off';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % select channels and spectra
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(g.plotchan) % find channel of minimum power
            [tmp indexfreq] = min(abs(g.freq-freqs));
            [tmp g.plotchan] = min(eegspecdB(:,indexfreq));
            fprintf('Channel %d has maximum power at %g\n', g.plotchan, g.freq);
            eegspecdBtoplot = eegspecdB(g.plotchan, :);		
        elseif g.plotchan == 0
            fprintf('Computing RMS power at all channels\n');
            eegspecdBtoplot = sqrt(mean(eegspecdB.^2,1)); % RMS before log as for components
        else 
            eegspecdBtoplot = eegspecdB;
        end;
        tmpc = find(eegspecdB(:,1)); 			% > 0 power chans
        zchans = int2str(find(eegspecdB(:,1) == 0)); 	% 0-power chans
        if length(tmpc) ~= size(eegspecdB,1)
            fprintf('\nWarning: channels [%s] have 0 values, so will be omitted from the display', ...
                       zchans);
            eegspecdB = eegspecdB(tmpc,:);
            if ~isempty(specstd),  specstd = specstd(tmpc,:); end
            if ~isempty(g.chanlocs)
                g.chanlocs2 = g.chanlocs(tmpc);
            end
        end;
        specstd   = 10*log10(specstd);
        eegspecdB = 10*log10(eegspecdB);
        eegspecdBtoplot = 10*log10(eegspecdBtoplot);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute component spectra
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        newweights = g.weights;
        if strcmp(g.memory, 'high') & strcmp(g.icamode, 'normal')
            fprintf('Computing component spectra: ')
            [compeegspecdB freqs] = spectcomp( newweights*data(:,:), frames, srate, epoch_subset, g);
        else % in case out of memory error, multiply conmponent sequencially
            if strcmp(g.icamode, 'sub') & ~isempty(g.plotchan) & g.plotchan == 0
                % scan all electrodes
                fprintf('Computing component spectra at each channel: ')
                for index = 1:size(data,1)
                    g.plotchan = index;
                    [compeegspecdB(:,:,index) freqs] = spectcomp( data, frames, srate, epoch_subset, g, newweights);
                end;
                g.plotchan = 0;
            else
                fprintf('Computing component spectra: ')
                [compeegspecdB freqs] = spectcomp( data, frames, srate, epoch_subset, g, newweights);
            end;
        end;
        fprintf('\n');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rescale spectra with respect to projection (rms = root mean square)
        % all channels: component_i power = rms(inverseweigths(component_i)^2)*power(activation_component_i);
        % one channel:  component_i power = inverseweigths(channel_j,component_i)^2)*power(activation_component_i);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmpi(g.icamode, 'normal')
            for index = 1:size(compeegspecdB,1) 
                if g.plotchan == 0 % normalize by component scalp map power
                    compeegspecdB(index,:) = 10*log10( sqrt(mean(g.icawinv(:,index).^4)) * compeegspecdB(index,:) );
                else 
                    compeegspecdB(index,:) = 10*log10( g.icawinv(g.plotchan,index)^2 * compeegspecdB(index,:) );
                end;
            end;
        else % already spectrum of data-components
            compeegspecdB = 10*log10( compeegspecdB );
        end;

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
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute axis and caxis g.limits
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

else % no freq specified

    if isnan(g.limits(2))
        g.limits(2) = srate/2;
    end;
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
% plot spectrum of each channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(g.plot, 'on')
    mainfig = gca; axis off;
    if ~isempty(g.freq)
        specaxes = sbplot(3,4,[5 12], 'ax', mainfig); 
    end;
    
    if isempty(g.weights)
        %pl=plot(freqs(1:maxfreqidx),eegspecdB(:,1:maxfreqidx)'); % old command
        if strcmpi(g.plotmean, 'on'), specdata = mean(eegspecdB,1); % average channels
        else                          specdata = eegspecdB;
        end;
        for index = 1:size(specdata,1) % scan channels
            tmpcol  = allcolors{mod(index, length(allcolors))+1};
            command = [ 'disp(''Channel ' int2str(index) ''')' ];
            pl(index)=plot(freqs(1:maxfreqidx),specdata(index,1:maxfreqidx)', ...
                           'color', tmpcol, 'ButtonDownFcn', command); hold on;
        end;
    else 
        for index = 1:size(eegspecdBtoplot,1)
            tmpcol  = allcolors{mod(index, length(allcolors))+1};
            command = [ 'disp(''Channel ' int2str(g.plotchan(index)) ''')' ];
            pl(index)=plot(freqs(1:maxfreqidx),eegspecdBtoplot(index,1:maxfreqidx)', ...
                           'color', tmpcol, 'ButtonDownFcn', command); hold on;
        end;
    end;
    set(pl,'LineWidth',2);
    set(gca,'TickLength',[0.02 0.02]);
    try, 
        axis([freqs(minfreqidx) freqs(maxfreqidx) reallimits(1) reallimits(2)]);
    catch, disp('Could not adjust axis'); end;
    xl=xlabel('Frequency (Hz)');
    set(xl,'fontsize',AXES_FONTSIZE_L);
    % yl=ylabel('Rel. Power (dB)');
    yl=ylabel('Power 10*log_{10}(\muV^{2}/Hz)');
    set(yl,'fontsize',AXES_FONTSIZE_L);
    set(gca,'fontsize',AXES_FONTSIZE_L)
    box off;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   plot component contribution   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colrs = {'r','b','g','m','c'}; % component spectra trace colors
if ~isempty(g.weights)
    if strcmpi(g.plot, 'on')
        if strcmpi(g.icamode, 'sub')
            set(pl,'LineWidth',5, 'color', 'k');
        else 
            set(pl, 'linewidth', 2, 'color', 'k');
        end;
        
        hold on;
        for f=1:length(g.icamaps)
            colr = colrs{mod((f-1),5)+1};
            command = [ 'disp(''Component ' int2str(g.icamaps(f)) ''')' ];
            pl2(index)=plot(freqs(1:maxfreqidx),compeegspecdB(g.icamaps(f),1:maxfreqidx)', ...
                            'color', colr, 'ButtonDownFcn', command); hold on;
        end
        othercomps = setdiff_bc(1:size(compeegspecdB,1), g.icamaps);
        if ~isempty(othercomps)
            for index = 1:length(othercomps)
                tmpcol  = allcolors{mod(index, length(allcolors))+1};
                command = [ 'disp(''Component ' int2str(othercomps(index)) ''')' ];
                pl(index)=plot(freqs(1:maxfreqidx),compeegspecdB(othercomps(index),1:maxfreqidx)', ...
                               'color', tmpcol, 'ButtonDownFcn', command); hold on;
            end;
        end;
        if length(g.limits)<3|isnan(g.limits(3))
            newaxis = axis;
            newaxis(3) = min(newaxis(3), min(min(compeegspecdB(:,1:maxfreqidx))));
            newaxis(4) = max(newaxis(4), max(max(compeegspecdB(:,1:maxfreqidx))));
            axis(newaxis);
        end;
	end;
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% indicate component contribution %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    maxdatadb = max(eegspecdBtoplot(:,freqidx(1)));
    [tmp indexfreq] = min(abs(g.freq-freqs));
    for index = 1:length(g.icacomps)
        if strcmp(g.icamode, 'normal')
            % note: maxdatadb = eegspecdBtoplot (RMS power of data)
            resvar(index)  = 100*exp(-(maxdatadb-compeegspecdB(index, indexfreq))/10*log(10));
            fprintf('Component %d percent relative variance:%6.2f\n', g.icacomps(index), resvar(index));
        else
            if g.plotchan == 0
                resvartmp = [];
                for chan = 1:size(eegspecdB,1) % scan channels
                    resvartmp(chan)  = 100 - 100*exp(-(eegspecdB(chan,freqidx(1))-compeegspecdB(index, indexfreq, chan))/10*log(10));
                end;
                resvar(index) = mean(resvartmp); % mean contribution for all channels
                stdvar(index) = std(resvartmp);
                fprintf('Component %d percent variance accounted for:%6.2f ± %3.2f\n', ...
                        g.icacomps(index), resvar(index), stdvar(index));
            else
                resvar(index)  = 100 - 100*exp(-(maxdatadb-compeegspecdB(index, indexfreq))/10*log(10));
                fprintf('Component %d percent variance accounted for:%6.2f\n', g.icacomps(index), resvar(index));
            end;
        end;
    end;
    
    % for icamode=sub and plotchan == 0 -> take the RMS power of all channels
    % -----------------------------------------------------------------------
    if ndims(compeegspecdB) == 3
        compeegspecdB = exp( compeegspecdB/10*log(10) );
        compeegspecdB = sqrt(mean(compeegspecdB.^2,3)); % RMS before log (dim1=comps, dim2=freqs, dim3=chans)
        compeegspecdB = 10*log10( compeegspecdB );
    end;
    
end;

if ~isempty(g.freq) &  strcmpi(g.plot, 'on')
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% plot vertical lines through channel trace bundle at each headfreq
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isempty(g.weights)
		for f=1:length(g.freq)
			hold on; 
			plot([freqs(freqidx(f)) freqs(freqidx(f))], ...
				 [min(eegspecdB(:,freqidx(f))) max(eegspecdB(:,freqidx(f)))],...
				 'k','LineWidth',2.5);
		end;
	else
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% plot vertical line at comp analysis freq
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		mincompdB = min([min(eegspecdB(:,freqidx(1))) min(compeegspecdB(:,freqidx(1)))]);
		maxcompdB = max([max(eegspecdB(:,freqidx(1))) max(compeegspecdB(:,freqidx(1)))]);
		hold on;
		plot([freqs(freqidx(1)) freqs(freqidx(1))],[mincompdB maxcompdB],'k', 'LineWidth',2.5);
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%
	% create axis for topoplot
	%%%%%%%%%%%%%%%%%%%%%%%%%%
	tmpmainpos = get(gca, 'position');
	headax = zeros(1,length(g.freq));
	for f=1:length(g.freq)+length(g.icamaps)
		headax(f) = sbplot(3,length(g.freq)+length(g.icamaps),f, 'ax', mainfig);
		axis([-1 1 -1 1]);
		
		%axis x coords and use
		tmppos = get(headax(f), 'position');
		allaxcoords(f) = tmppos(1);
		allaxuse(f)    = 0;
	end
	large = sbplot(1,1,1, 'ax', mainfig);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% compute relative positions on plot
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if  ~isempty(g.weights)
		freqnormpos = tmpmainpos(1) + tmpmainpos(3)*(freqs(freqidx(1))-g.limits(1))/(g.limits(2)-g.limits(1));
		for index = 1:length(g.icamaps)+1
			[realpos(index) allaxuse] = closestplot( freqnormpos, allaxcoords, allaxuse );
		end;
	
		% put the channel plot a liitle bit higher
		tmppos = get(headax(realpos(1)), 'position');
		tmppos(2) = tmppos(2)+0.04;
		set(headax(realpos(1)), 'position', tmppos);
	else 
		realpos = 1:length(g.freq); % indices giving order of plotting positions
	end;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% plot connecting lines using changeunits()
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for f=1:length(g.freq)+length(g.icamaps)
		if ~isempty(g.weights)
			if f>length(g.freq) % special case of ICA components
				from = changeunits([freqs(freqidx(1)),g.icafreqsval(f-1)],specaxes,large);
				%g.icafreqsval contains the sorted frequency values at the specified frequency
			else 
				from = changeunits([freqs(freqidx(f)),maxcompdB],specaxes,large);
			end;
		else
			from = changeunits([freqs(freqidx(f)),max(eegspecdB(:,freqidx(f)))],specaxes,large);
		end;
		pos = get(headax(realpos(f)),'position');
		to = changeunits([0,0],headax(realpos(f)),large)+[0 -min(pos(3:4))/2.5];
		hold on;
		if f<=length(g.freq)
			colr = 'k';
		else
			colr = colrs{mod((f-2),5)+1};
		end
		li(realpos(f)) = plot([from(1) to(1)],[from(2) to(2)],colr,'LineWidth',PLOT_LINEWIDTH_S);
		axis([0 1 0 1]);
		axis off;
	end;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% plot selected channel head using topoplot()
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	fprintf('Plotting scalp distributions: ')
	for f=1:length(g.freq)
		axes(headax(realpos(f)));
        
  		topodata = eegspecdB(:,freqidx(f))-nan_mean(eegspecdB(:,freqidx(f)));

		if isnan(g.limits(5)),     
			maplimits = 'absmax';
		else                       
			maplimits = [g.limits(5) g.limits(6)];
		end;
        
		%
		% If 1 channel in g.plotchan
		%

        if ~isempty(g.plotchan) & g.plotchan ~= 0 
			% if ~isempty(varargin) % if there are extra topoplot() flags
			%	topoplot(g.plotchan,g.chanlocs,'electrodes','off', ...
			%			 'style', 'blank', 'emarkersize1chan', 10, varargin{:});
			% else
				topoplot(g.plotchan,g.chanlocs,'electrodes','off', ...
						 'style', 'blank', 'emarkersize1chan', 10);
			% end
			if isstruct(g.chanlocs)
				tl=title(g.chanlocs(g.plotchan).labels);
			else
				tl=title([ 'c' int2str(g.plotchan)]);
			end;
            
		else % plot all channels in g.plotchans 

            if isempty(g.mapframes) | g.mapframes == 0
                g.mapframes = 1:size(eegspecdB,1); % default to plotting all chans
            end
			if ~isempty(varargin)
				topoplot(topodata(g.mapframes),g.chanlocs2,'maplimits',maplimits, varargin{:}); 
			else
				topoplot(topodata(g.mapframes),g.chanlocs2,'maplimits',maplimits); 
			end
			if f<length(g.freq)
				tl=title([num2str(freqs(freqidx(f)), '%3.1f')]);
			else
				tl=title([num2str(freqs(freqidx(f)), '%3.1f') ' Hz']);
			end
		end;
		set(tl,'fontsize',AXES_FONTSIZE_L);
		axis square;
		drawnow
		fprintf('.');
	end;
	fprintf('\n');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% plot independent components
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if ~isempty(g.weights)
		% use headaxe from 2 to end (reserved earlier)
		if realpos(1) == max(realpos), plotcolbar(g); end;

		% make the line with the scalp topoplot thicker than others
		set(li(realpos(1)), 'linewidth', 2.5); 
        if isempty(g.mapchans) | g.mapchans == 0
            g.mapchans = 1:size(g.icawinv,1); % default to plotting all chans
        end
		for index = 1:length(g.icamaps)
			axes(headax(realpos(index+1)));						
			compnum = g.icamaps(index);

			topoplot(g.icawinv(g.mapchans,compnum).^2,g.chanlocs,varargin{:}); 
			tl=title(int2str(g.icacomps(compnum)));
			set(tl,'fontsize',16);
			axis square;
			drawnow
            try,
                if strcmpi(g.icamode, 'normal')
                    set(gca, 'userdata', ['text(-0.6, -0.6, ''Rel. Var.: ' sprintf('%6.2f', resvar(g.icacomps(compnum))) ''');'] );
                else
                    set(gca, 'userdata', ['text(-0.6, -0.6, ''PVAF: ' sprintf('%6.2f', resvar(g.icacomps(compnum))) ''');'] );
                end;
            catch, end;
			if realpos(index+1) == max(realpos), plotcolbar(g); end;
		end;
	else 
		plotcolbar(g);
	end;
end;

%%%%%%%%%%%%%%%%
% Draw title
%%%%%%%%%%%%%%%%
if ~isempty(g.title) & strcmpi(g.plot, 'on')
    tl = textsc(g.title);
	set(tl,'fontsize',15)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return component spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if ~isempty(g.weights) & nargout >= 3
    tmp = compeegspecdB;
    compeegspecdB = zeros(ncompsori, size(tmp,2));
    compeegspecdB(g.icacomps,:) = tmp;
end;
    
%%%%%%%%%%%%%%%%
% Turn on axcopy (disabled to allow to click on curves)
%%%%%%%%%%%%%%%%
if strcmpi(g.plot, 'on')
    disp('Click on each trace for channel/component index');
    axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end;');
    % will not erase the commands for the curves
end;

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
			set(cb,'yticklabel',strvcat('-',' ','+'));
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
function [eegspecdB, freqs, specstd] = spectcomp( data, frames, srate, epoch_subset, g, newweights);
	if exist('newweights') == 1
		nchans = size(newweights,1);
	else 
		nchans = size(data,1);		
	end;
	%fftlength = 2^round(log(srate)/log(2))*g.freqfac;
	if isempty(g.winsize)
        winlength = max(pow2(nextpow2(frames)-3),4); %*2 since diveded by 2 later	
        winlength = min(winlength, 512);
        winlength = max(winlength, 256);
        winlength = min(winlength, frames);
    else
        winlength = g.winsize;
    end;
    if isempty(g.nfft)
        fftlength = 2^(nextpow2(winlength))*g.freqfac;
	else
        fftlength = g.nfft;
    end;
%     usepwelch = 1; 
    usepwelch = license('checkout','Signal_Toolbox'); % 5/22/2014 Ramon 
%     if ~license('checkout','Signal_Toolbox'), 
    if ~usepwelch, 
        fprintf('\nSignal processing toolbox (SPT) absent: spectrum computed using the pwelch()\n');
        fprintf('function from Octave which is suposedly 100%% compatible with the Matlab pwelch function\n');
    end;
    fprintf(' (window length %d; fft length: %d; overlap %d):\n', winlength, fftlength, g.overlap);	
        
	for c=1:nchans % scan channels or components
		if exist('newweights') == 1
			if strcmp(g.icamode, 'normal')
				tmpdata = newweights(c,:)*data; % component activity
			else % data - component contribution
                tmpdata = data(g.plotchan,:) - (g.icawinv(g.plotchan,c)*newweights(c,:))*data;
			end;
		else
			tmpdata = data(c,:); % channel activity
		end;
		if strcmp(g.reref, 'averef')
			tmpdata = averef(tmpdata);
		end;
		for e=epoch_subset
			if isempty(g.boundaries)
                if usepwelch
                    [tmpspec,freqs] = pwelch(matsel(tmpdata,frames,0,1,e),...
                                             winlength,g.overlap,fftlength,srate);
                else
                    [tmpspec,freqs] = spec(matsel(tmpdata,frames,0,1,e),fftlength,srate,...
                                           winlength,g.overlap);
                end;
				%[tmpspec,freqs] = psd(matsel(tmpdata,frames,0,1,e),fftlength,srate,...
				%					  winlength,g.overlap);
				if c==1 & e==epoch_subset(1)
					eegspec = zeros(nchans,length(freqs));
					specstd = zeros(nchans,length(freqs));
				end
				eegspec(c,:) = eegspec(c,:) + tmpspec';
				specstd(c,:) = specstd(c,:) + tmpspec'.^2;
            else
                g.boundaries = round(g.boundaries);
				for n=1:length(g.boundaries)-1
                    if g.boundaries(n+1) - g.boundaries(n) >= winlength % ignore segments of less than winlength
                        if usepwelch
                            [tmpspec,freqs] =  pwelch(tmpdata(e,g.boundaries(n)+1:g.boundaries(n+1)),...
                                                      winlength,g.overlap,fftlength,srate);
                        else
                            [tmpspec,freqs] =  spec(tmpdata(e,g.boundaries(n)+1:g.boundaries(n+1)),...
                                                    fftlength,srate,winlength,g.overlap);
                        end;
                        if exist('eegspec') ~= 1
                            eegspec = zeros(nchans,length(freqs));
                            specstd = zeros(nchans,length(freqs));
                        end
                        eegspec(c,:) = eegspec(c,:) + tmpspec'* ...
                            ((g.boundaries(n+1)-g.boundaries(n)+1)/g.boundaries(end));
                        specstd(c,:) = eegspec(c,:) + tmpspec'.^2 * ...
                            ((g.boundaries(n+1)-g.boundaries(n)+1)/g.boundaries(end));
                    end;
				end
			end;
		end
		fprintf('.')
	end
    
    n = length(epoch_subset);
	eegspecdB = eegspec/n; % normalize by the number of sections
    if n>1  % normalize standard deviation by the number of sections
        specstd   = sqrt( (specstd +  eegspec.^2/n)/(n-1) ); 
    else specstd   = [];
    end;
	return;
