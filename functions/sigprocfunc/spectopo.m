% spectopo() - Plot the mean log spectrum of a set of data epochs at
%              all channels as a bundle of traces. At specified frequencies,
%              plot the relative topographic distribution of power.
%              Uses Matlab psd() from signal processing toolbox.
% Usage:
%        >> [spectra,freqs] = spectopo(data,frames,srate,freqs,chanlocs,...
%                           limits,title,freqfaq, percent, 'key1', 'value1' ...);
% Inputs:
%        data   = 2D (nchans,frames*epochs); % can be single-epoch
%                 or 3D (nbchan,frames,epochs)
%        frames = frames per epoch {0 -> data length}
%        srate  = sampling rate per channel (Hz)
%        freqs  = vector of frequencies for topoplot() scalp maps (Hz)
%        chanlocs = electrode locations file (format: >> topoplot example)
%
% Optional inputs:
%        limits  = axis limits [xmin xmax ymin ymax cmin cmax]
%                  To use data limtis, omit final values or use nan's
%                  i.e. [-100 900 nan nan -10 10], [-100 900]
%                  Note that default color limits are symmetric around 0 and are
%                  different for each head {defaults: all nans}
%        title   = quoted plot title {default: none}
%        freqfac = int power of 2 => approximate frequencies/Hz {default: 4}
%        percent = downsampling factor or approximate percentage of the data to
%                  keep while computing spectra. Downsampling can help to speed up
%                  the computation. From 0 to 1 {default: 1}
%        'key','val' = optional topoplot() arguments (see topoplot())
%
% Outputs:
%        spectra = (nchans,nfreqs) power spectra (average over epochs) in dB
%        freqs   = frequencies of spectra (Hz)
%
% Authors: Scott Makeig & Marissa Westerfield, SCCN/INC/UCSD, La Jolla, 3/01 
%
% See also: timtopo(), envtopo(), tftopo(), topoplot()

% Copyright (C) 3/01 Scott Makeig & Marissa Westerfield, SCCN/INC/UCSD, 
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

% Uses: psd(), changeunits(), topoplot(), textsc()

function [eegspecdB,freqs]=spectopo(data,frames,srate,headfreqs,chanlocs,limits,titl,freqfac, percent, varargin)

LOPLOTHZ = 1;  % low  Hz to plot
FREQFAC  = 4;  % approximate frequencies/Hz (default)

if nargin<5
   help spectopo
   return
end

if nargin<8
  freqfac = FREQFAC;
end

if nargin<7
  titl = []; % default no title
end

if nargin<6
  limits = [nan nan nan nan nan nan]; % defaults from data 
end

if nargin<9
  percent = 1; % defaults sample 100% of the data 
else
    percent = max(percent, 0);
    percent = min(percent, 1);
end
data = reshape(data, size(data,1), size(data,2)*size(data,3));
if frames == 0
  frames = size(data,2); % assume one epoch
end

if min(headfreqs)<0
   fprintf('spectopo(): freqs must be >=0 Hz\n');
   return
end
nchans = size(data,1);
epochs = round(size(data,2)/frames);
if frames*epochs ~= size(data,2)
   error('Spectopo: non-integer number of epochs');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute channel spectra using psd()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epoch_subset = 1:epochs;
if percent ~= 1 & epochs > 1
    nb = round( percent*epochs);
    epoch_subset = zeros(1,epochs);
    while nb>0
        index = ceil(rand*epochs);
        if ~epoch_subset(index)
            epoch_subset(index) = 1;
            nb = nb-1;
        end;
    end;        
    epoch_subset = find(epoch_subset == 1);
    fprintf('Selecting randomly %d epochs\n', length(epoch_subset));
end;
fftlength = 2^round(log(srate)/log(2))*FREQFAC;
fprintf('Computing spectra: ')
for c=1:nchans
  for e=epoch_subset
    [tmpspec,freqs] = psd(matsel(data,frames,0,c,e),fftlength,...
                                    srate,fftlength/2,fftlength/8);
     if c==1 & e==epoch_subset(1)
       eegspec = zeros(nchans,length(freqs));
     end
     eegspec(c,:) = eegspec(c,:) + tmpspec';
  end
  fprintf('.')
end
eegspecdB = 10*log10(eegspec/epochs); % convert power to dB
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute axis and caxis limits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(limits)<1 | isnan(limits(1))
   limits(1) = LOPLOTHZ;
end

if length(limits)<2 | isnan(limits(2))
   maxheadfreq = max(headfreqs);
   if rem(maxheadfreq,5) ~= 0
     limits(2) = 5*ceil(maxheadfreq/5);
   else
     limits(2) = maxheadfreq*1.1;
   end
end

headfreqs = sort(headfreqs);          % Determine topoplot frequencies
freqidx = zeros(1,length(headfreqs)); % Do not interpolate between freqs
for f=1:length(headfreqs)
   [tmp fi] = min(abs(freqs-headfreqs(f)));
   freqidx(f)=fi;
end
[tmp maxfreqidx] = min(abs(limits(2)-freqs)); % adjust max frequency
[tmp minfreqidx] = min(abs(limits(1)-freqs)); % adjust min frequency

if length(limits)<3|isnan(limits(3))
  limits(3) = min(min(eegspecdB(:,minfreqidx:maxfreqidx)));
end
if length(limits)<4|isnan(limits(4))
  limits(4) = max(max(eegspecdB(:,minfreqidx:maxfreqidx)));
end
dBrange = limits(4)-limits(3);   % expand range a bit beyond data limits
limits(3) = limits(3)-dBrange/7;
limits(4) = limits(4)+dBrange/7;

if length(limits)<5 % default caxis plotting limits
  limits(5) = nan;
end
if length(limits)<6 
  limits(6) = nan;
end

if isnan(limits(5))+isnan(limits(6)) == 1
   fprintf('spectopo(): limits 5 and 6 must either be given or nan\n');
   return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot spectrum of each channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
specaxes = sbplot(3,4,[5 12]); 
% plot(freqs(1:maxfreqidx),eegspecdB(:,1:maxfreqidx)','b','LineWidth',2);
pl=plot(freqs(1:maxfreqidx),eegspecdB(:,1:maxfreqidx)');
set(pl,'LineWidth',2);
set(gca,'TickLength',[0.02 0.02]);

axis([freqs(minfreqidx) freqs(maxfreqidx) limits(3) limits(4)]);

xl=xlabel('Frequency (Hz)');
set(xl,'fontsize',16);
yl=ylabel('Rel. Power (dB)');
set(yl,'fontsize',16);
set(gca,'fontsize',16)
box off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot lines through channel trace bundle at each headfreq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for f=1:length(headfreqs)
   hold on; 
   plot([freqs(freqidx(f)) freqs(freqidx(f))], ...
        [min(eegspecdB(:,freqidx(f))) max(eegspecdB(:,freqidx(f)))],...
               'k','LineWidth',2.5);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot connecting lines using changeunits()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
headax = zeros(1,length(headfreqs));
for f=1:length(headfreqs) 
   headax(f) = sbplot(3,length(headfreqs),f);
   axis([-1 1 -1 1]);
end
large = sbplot(1,1,1);
for f=1:length(headfreqs)
   from = changeunits([freqs(freqidx(f)),max(eegspecdB(:,freqidx(f)))],...
                       specaxes,large);
   to = changeunits([0,0],headax(f),large);
   hold on;
   plot([from(1) to(1)],[from(2) to(2)],'k','LineWidth',2);
   axis([0 1 0 1]);
   axis off;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot heads using topoplot()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Plotting scalp distributions: ')
for f=1:length(headfreqs) 
   axes(headax(f));
   topodata = eegspecdB(:,freqidx(f))-mean(eegspecdB(:,freqidx(f)));
   if isnan(limits(5)),     maplimits = 'absmax';
   else                     maplimits = [limits(5) limits(6)];
   end;
   if ~isempty(varargin)
     topoplot(topodata,chanlocs,'maplimits',maplimits, varargin{:}); 
   else
     topoplot(topodata,chanlocs,'maplimits',maplimits); 
   end
   if f<length(headfreqs)
     tl=title([num2str(freqs(freqidx(f)), '%3.1f')]);
   else
     tl=title([num2str(freqs(freqidx(f)), '%3.1f') ' Hz']);
   end
   set(tl,'fontsize',16);
   axis square;
   drawnow
   fprintf('.');
end;
fprintf('\n');

%%%%%%%%%%%%%%%%
% Plot color bar
%%%%%%%%%%%%%%%%
cb=cbar;
pos = get(cb,'position');
set(cb,'position',[pos(1) pos(2) 0.03 pos(4)]);
set(cb,'fontsize',12);
if isnan(limits(5))
   ticks = get(cb,'ytick');
   [tmp zi] = find(ticks == 0);
   ticks = [ticks(1) ticks(zi) ticks(end)];
   set(cb,'ytick',ticks);
   set(cb,'yticklabel',{'-','0','+'});
end

%%%%%%%%%%%%%%%%
% Draw title
%%%%%%%%%%%%%%%%
if ~isempty(titl)
  tl = textsc(titl,'title');
  set(tl,'fontsize',15)
end

%%%%%%%%%%%%%%%%
% Turn on axcopy
%%%%%%%%%%%%%%%%
axcopy
