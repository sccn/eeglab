% tftopo()  - Generates a figure showing a selected image (e.g., an ERSP or ITC) from 
%             a supplied set of images for every scalp channel, plus topoplot() scalp 
%             maps at specified (x,y) (e.g., time,frequency) points.  Else, images the 
%             signed (selected) channel std(). Inputs may be outputs of timef(); else,
%             e.g., can be used to image a set of smoothed erpimage() images.
% Usage:
%        >> tftopo(tfdata,times,freqs, 'key1', 'val1', 'key2', val2' ...)
% Inputs:
%   tfdata    = Set of nchans time/freq ERSPs or ITCs from timef() (or any other
%               time/freq matrix), one for each channel. Size (time,freq,chans) or
%               (time,freq,chans,subjects) for grand RMS.
%   times     = Vector of times in msec from timef()
%   freqs     = Vector of frequencies in Hz from timef() 
%
% Optional inputs:
%  'timefreqs' = Array of time/frequency points to show topoplot() maps for
%                Format: size (nrows,2), each row [ms Hz]
%  'showchan'  = [integer] Channel number of tfdata to image, or 0
%                {default=0 -> image the median-signed st. dev. across channels} 
%  'chanlocs'  = ['string'|structure] Electrode locations file (for format see 
%                >> topoplot example) or structure             {default none}
%  'limits'    = Vector of plotting limits [minms maxms minhz maxhz mincaxis maxcaxis]
%                Omit, or use nan's to use tfdata limits. Ex: [nan nan -100 400];
%  'signifs'   = Significance level(s) (e.g., from timef()), for zero'ing non-significant 
%                tfdata. Size must be (1 or 2, freq, chans, subjects) or 
%                (1 or 2, freq, times, chans, subjects). If first dimension is
%                of size 1, tfdata is assumed to contain positive values   {default: none}
%  'sigthresh' = [integer], i.e. [K L] after masking time-frequency decomposition using 
%                'signifs' array, concatenate time/freq values only if more than K electrodes
%                have non-0 (significant) values. If several subject, the second value L
%                is used to concatenate subject in the same fashion. Default is [1 1].
%  'selchans'  = Channels to include in topoplot() scalp maps (and image std()) {default: all}
%  'smooth'    = [pow2] magnification and smoothing factor. Only in power of 2 (Default: 1}.
%  'mode'      = ['rms'|'ave'] average dB power ('ave') or return root-mean-square ('rms')
%                {Default: 'rms' }
%  'logfreq'   = ['on'|'off'] use logarithmic base for frequencies { Default: 'off' }
%  'vert'      = [times vector] -> plot vertical dashed lines at specified times in ms.
%                { Default: 0 }
%  'shiftimgs' = [time vector] - shift time frequency images (for several subject) by the
%                median reaction time of each subject (provided here).
%
% Note:
%  1) Additional topoplot() options can be used.
%  2) For topoplot maps, the average power (not masked by significance is used (instead
%     of the root-mean-square (RMS) for averaging electrode activity).
%  3) If several subjects (4-D tfdata input) RMS is first computed across electrodes
%     then across subjects.
%
% Authors: Scott Makeig, Arnaud Delorme & Marissa Westerfield, SCCN/INC/UCSD, La Jolla, 3/01 
%
% See also: timef(), topoplot(), spectopo(), timtopo(), envtopo(), changeunits()

% hidden parameter: 'shiftimgs' = array with one value per subject for shifting in time
%                                 time/freq images. Had to be inserted in tftopo because
%                                 the shift happen after the smoothing

% Copyright (C) Scott Makeig, Arnaud Delorme & Marissa Westerfield, SCCN/INC/UCSD, 
% La Jolla, 3/01
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
% Revision 1.63  2003/09/11 18:59:34  arno
% allowing image significance
%
% Revision 1.62  2003/08/28 17:41:21  arno
% debuging the first dimension for significance
%
% Revision 1.61  2003/08/27 21:49:48  arno
% fixing mode error
%
% Revision 1.60  2003/08/27 00:31:09  arno
% fix rms/average
%
% Revision 1.59  2003/03/08 02:26:32  arno
% document the shiftimgs option
%
% Revision 1.58  2003/03/08 02:23:54  arno
% update header
%
% Revision 1.57  2003/03/04 18:12:07  arno
% color limits debug
%
% Revision 1.56  2002/10/09 18:20:21  arno
% line thickness
%
% Revision 1.55  2002/10/09 18:12:00  arno
% new inputs logfreq, shiftdims, vert, smooth ...
%
% Revision 1.54  2002/10/08 23:55:35  arno
% debug
%
% Revision 1.53  2002/10/08 23:54:34  arno
% 'key', 'val' sequence, new args, new header
%
% Revision 1.52  2002/10/07 18:51:11  arno
% all data -> 3D, significance, RMS, subject RMS, ...
%
% Revision 1.51  2002/10/07 15:44:48  arno
% updating header and axis off
%
% Revision 1.50  2002/10/03 23:39:44  scott
% added sbplot() plotting -sm & ad
%
% Revision 1.49  2002/05/19 14:30:15  scott
% *** empty log message ***
%
% Revision 1.48  2002/05/19 14:28:55  scott
% improved out-of-bounds testing -sm
%
% Revision 1.47  2002/05/19 14:21:52  scott
% *** empty log message ***
%
% Revision 1.46  2002/05/19 14:20:34  scott
% *** empty log message ***
%
% Revision 1.45  2002/05/19 14:19:35  scott
% added cartoon of showchan if > 0  -sm
%
% Revision 1.44  2002/05/19 14:17:55  scott
% *** empty log message ***
%
% Revision 1.43  2002/05/19 14:16:11  scott
% *** empty log message ***
%
% Revision 1.42  2002/05/19 14:15:03  scott
% *** empty log message ***
%
% Revision 1.41  2002/05/19 14:12:53  scott
% *** empty log message ***
%
% Revision 1.40  2002/05/19 14:12:08  scott
% *** empty log message ***
%
% Revision 1.39  2002/05/19 14:07:12  scott
% *** empty log message ***
%
% Revision 1.38  2002/05/19 14:06:40  scott
% *** empty log message ***
%
% Revision 1.37  2002/05/19 14:05:50  scott
% testing -sm
%
% Revision 1.36  2002/05/19 14:00:45  scott
% *** empty log message ***
%
% Revision 1.35  2002/05/19 13:57:54  scott
% *** empty log message ***
%
% Revision 1.34  2002/05/19 13:52:24  scott
% showchans=0 -> image std() of selchans images -sm
%
% Revision 1.33  2002/05/19 13:35:16  scott
% *** empty log message ***
%
% Revision 1.32  2002/05/19 13:34:12  scott
% showchan==0 -> image signed st dev -sm
%
% Revision 1.31  2002/05/19 13:26:10  scott
% adjusted channel label -sm
%
% Revision 1.30  2002/05/19 03:04:46  scott
% *** empty log message ***
%
% Revision 1.29  2002/05/19 02:57:24  scott
% *** empty log message ***
%
% Revision 1.28  2002/05/19 02:55:45  scott
% *** empty log message ***
%
% Revision 1.27  2002/05/19 02:53:45  scott
% *** empty log message ***
%
% Revision 1.26  2002/05/19 02:50:40  scott
% *** empty log message ***
%
% Revision 1.25  2002/05/19 02:28:08  scott
% *** empty log message ***
%
% Revision 1.24  2002/05/19 02:24:15  scott
% *** empty log message ***
%
% Revision 1.23  2002/05/19 02:21:17  scott
% *** empty log message ***
%
% Revision 1.22  2002/05/19 02:20:26  scott
% *** empty log message ***
%
% Revision 1.21  2002/05/19 02:18:36  scott
% *** empty log message ***
%
% Revision 1.20  2002/05/19 02:15:13  scott
% adding separate scale for showchan==0 -sm
%
% Revision 1.19  2002/04/30 21:24:48  scott
% *** empty log message ***
%
% Revision 1.18  2002/04/30 21:23:58  scott
% *** empty log message ***
%
% Revision 1.17  2002/04/30 21:22:46  scott
% *** empty log message ***
%
% Revision 1.16  2002/04/30 21:21:50  scott
% *** empty log message ***
%
% Revision 1.15  2002/04/30 21:21:15  scott
% *** empty log message ***
%
% Revision 1.14  2002/04/30 21:19:05  scott
% debugging sign feature for showchans==0 -sm
%
% Revision 1.13  2002/04/30 21:17:59  scott
% fg
%
% Revision 1.12  2002/04/30 21:17:01  scott
% adding sign -sm
%
% Revision 1.11  2002/04/30 20:53:35  scott
% debugging showchans==0 option -sm
%
% Revision 1.10  2002/04/30 20:47:39  scott
% *** empty log message ***
%
% Revision 1.9  2002/04/30 20:45:56  scott
% showchan==0 -> blockave(abs(tfdata)) -sm
%
% Revision 1.8  2002/04/27 01:37:19  scott
% same -sm
%
% Revision 1.7  2002/04/27 01:26:40  scott
% updated topoplot call -sm
%
% Revision 1.6  2002/04/27 01:19:33  scott
% same -sm
%
% Revision 1.5  2002/04/27 01:13:46  scott
% same -sm
%
% Revision 1.4  2002/04/27 01:10:29  scott
% same -sm
%
% Revision 1.3  2002/04/27 01:06:00  scott
% same -sm
%
% Revision 1.2  2002/04/27 01:04:12  scott
% added handling of 3-d tftopo data -sm
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 

function tfave = tftopo(tfdata,times,freqs,varargin);
    %timefreqs,showchan,chanlocs,limits,signifs,selchans)

LINECOLOR= 'k';
LINEWIDTH = 2;
ZEROLINEWIDTH = 2.8;

if nargin<3
   help tftopo
   return
end

% reshape tfdata
% --------------
if length(size(tfdata))==2
   nchans = round(size(tfdata,2)/length(times));
   tfdata = reshape(tfdata, size(tfdata,1), length(times), nchans); 
elseif length(size(tfdata))>=3
   nchans = size(tfdata,3);
else
   help tftopo
   return
end
tfdataori = mean(tfdata,4); % for topoplot

% test inputs
% -----------
% 'key' 'val' sequence
fieldlist = { 'chanlocs'      { 'string' 'struct' }       []       '' ;
              'limits'        'real'     []                        [nan nan nan nan nan nan];
              'logfreq'       'string'   {'on' 'off' }             'off';
              'mode'          'string'   { 'ave' 'rms' }           'rms';
              'selchans'      'integer'  [1 nchans]                [1:nchans];
              'shiftimgs'     'real'     []                        [] ;
              'showchan'      'integer'  [0 nchans]                0 ;
              'signifs'       'real'     []                        [];
              'sigthresh'     'integer'  [1 Inf]                   [1 1];
              'smooth'        'real'     [0 Inf]                   1;
              'timefreqs'     'real'     []                        [];
              'vert'          'real'     [times(1) times(end)]     [0] };

[g varargin] = finputcheck( varargin, fieldlist, 'tftopo', 'ignore');
if isstr(g), error(g); end;

% setting more defaults
% ---------------------
if length(times) ~= size(tfdata,2)
   fprintf('tftopo(): tfdata columns must be a multiple of the length of times (%d)\n',...
                 length(times));
   return
end
if length(g.showchan) > 1
    error('tftopo(): showchan must be a single number');
end;
if length(g.limits)<1 | isnan(g.limits(1))
  g.limits(1) = times(1);
end
if length(g.limits)<2 | isnan(g.limits(2))
  g.limits(2) = times(end);
end
if length(g.limits)<3 | isnan(g.limits(3))
  g.limits(3) = freqs(1);
end
if length(g.limits)<4 | isnan(g.limits(4))
  g.limits(4) = freqs(end);
end
if length(g.limits)<5 | isnan(g.limits(5)) % default caxis plotting limits
  g.limits(5) = -max(abs(tfdata(:)));
  mincax = g.limits(5); 
end
if length(g.limits)<6 | isnan(g.limits(6))
    defaultlim = 1;
    if exist('mincax')
        g.limits(6) = -mincax; % avoid recalculation
    else
        g.limits(6) = max(abs(tfdata(:)));
    end
else 
    defaultlim = 0;
end
if length(g.sigthresh) == 1
    g.sigthresh(2) = 1;
end;
if g.sigthresh(1) > nchans
    error('tftopo(): ''sigthresh'' first number must be lower or equal to the number of channels');
end;
if g.sigthresh(2) > size(tfdata,4)
    error('tftopo(): ''sigthresh'' second number must be lower or equal to the number of subjects');
end;
if ~isempty(g.signifs)
    if size(g.signifs,1) > 2 | size(g.signifs,2) ~= size(tfdata,1)| ...
            (size(g.signifs,3) ~= size(tfdata,3) & size(g.signifs,4) ~= size(tfdata,3))
        fprintf('tftopo(): error in ''signifs'' array size not compatible with data size.\n');
        return
    end
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process time/freq data points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(g.timefreqs)
    if isempty(g.chanlocs)
        error('tftopo(): ''chanlocs'' must be defined to plot time/freq points');
    end;
    if min(g.timefreqs(:,2))<min(freqs) 
        fprintf('tftopo(): selected plotting frequency %g out of range.\n',min(g.timefreqs(:,2)));
        return
    end
    if max(g.timefreqs(:,2))>max(freqs) 
        fprintf('tftopo(): selected plotting frequency %g out of range.\n',max(g.timefreqs(:,2)));
        return
    end
    if min(g.timefreqs(:,1))<min(times) 
        fprintf('tftopo(): selected plotting time %g out of range.\n',min(g.timefreqs(:,1)));
        return
    end
    if max(g.timefreqs(:,1))>max(times) 
        fprintf('tftopo(): selected plotting time %g out of range.\n',max(g.timefreqs(:,1)));
        return
    end

    if 0 % USE USER-SUPPLIED SCALP MAP ORDER. A GOOD ALGORITHM FOR SELECTING
         % g.timefreqs POINT ORDER GIVING MAX UNCROSSED LINES IS DIFFICULT!
        [tmp tfi] = sort(g.timefreqs(:,1)); % sort on times
        tmp = g.timefreqs;
        for t=1:size(g.timefreqs,1)
            g.timefreqs(t,:) = tmp(tfi(t),:);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute timefreqs point indices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tfpoints = size(g.timefreqs,1);
    freqidx = zeros(1,tfpoints);
    for f=1:tfpoints
        [tmp fi] = min(abs(freqs-g.timefreqs(f,2)));
        freqidx(f)=fi;
    end
    timeidx = zeros(1,tfpoints);
    for f=1:tfpoints
        [tmp fi] = min(abs(times-g.timefreqs(f,1)));
        timeidx(f)=fi;
    end
    tfpidx = [timeidx' freqidx'];
else 
    tfpoints = 0;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zero out non-significant image features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
range = g.limits(6)-g.limits(5);
cc = jet(256);
if ~isempty(g.signifs)
    fprintf('Applying ''signifs'' mask by zeroing non-significant values\n');
    for subject = 1:size(tfdata,4)
        for elec = 1:size(tfdata,3)
            
            if size(g.signifs,1) == 2
                if ndims(g.signifs) > ndims(tfdata)
                    tmpfilt = (tfdata(:,:,elec,subject) >= squeeze(g.signifs(2,:,:,elec, subject))') | ...
                              (tfdata(:,:,elec,subject) <= squeeze(g.signifs(1,:,:,elec, subject))');
                else
                    tmpfilt = (tfdata(:,:,elec,subject) >= repmat(g.signifs(2,:,elec, subject)', [1 size(tfdata,2)])) | ...
                              (tfdata(:,:,elec,subject) <= repmat(g.signifs(1,:,elec, subject)', [1 size(tfdata,2)]));
                end;
            else
                if ndims(g.signifs) > ndims(tfdata)
                    tmpfilt = (tfdata(:,:,elec,subject) >= squeeze(g.signifs(1,:,:,elec, subject))');
                else
                    tmpfilt = (tfdata(:,:,elec,subject) >= repmat(g.signifs(1,:,elec, subject)', [1 size(tfdata,2)]));
                end;
            end;                
            tfdata(:,:,elec,subject) = tfdata(:,:,elec,subject) .* tmpfilt;
        end;
    end;
end;

%%%%%%%%%%%%%%%%
% magnify inputs
%%%%%%%%%%%%%%%%
if g.smooth ~= 1
    fprintf('Smoothing...\n');
    for index = 1:round(log2(g.smooth))
        [tfdata times freqs] = magnifytwice(tfdata, times, freqs);
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%
% Shift time/freq images
%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(g.shiftimgs)
    timestep = times(2) - times(1);
    for S = 1:size(tfdata,4)
        nbsteps = round(g.shiftimgs(S)/timestep);
        fprintf('Shifing images of subect %d by %3.3f ms or %d time steps\n', S, g.shiftimgs(S), nbsteps);
        if nbsteps < 0,  tfdata(:,-nbsteps+1:end,:,S) = tfdata(:,1:end+nbsteps,:,S);
        else             tfdata(:,1:end-nbsteps,:,S)  = tfdata(:,nbsteps+1:end,:,S);
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust plotting limits
%%%%%%%%%%%%%%%%%%%%%%%%%
[tmp minfreqidx] = min(abs(g.limits(3)-freqs)); % adjust min frequency
 g.limits(3) = freqs(minfreqidx);
[tmp maxfreqidx] = min(abs(g.limits(4)-freqs)); % adjust max frequency
 g.limits(4) = freqs(maxfreqidx);

[tmp mintimeidx] = min(abs(g.limits(1)-times)); % adjust min time
 g.limits(1) = times(mintimeidx);
[tmp maxtimeidx] = min(abs(g.limits(2)-times)); % adjust max time
 g.limits(2) = times(maxtimeidx);

mmidx = [mintimeidx maxtimeidx minfreqidx maxfreqidx];

%colormap('jet');
%c = colormap;
%cc = zeros(256,3);
%if size(c,1)==64
%    for i=1:3
%       cc(:,i) = interp(c(:,i),4);
%    end
%else
%    cc=c;
%nd
%cc(find(cc<0))=0;
%cc(find(cc>1))=1;

%if exist('g.signif')
%  minnull = round(256*(g.signif(1)-g.limits(5))/range);
%  if minnull<1
%    minnull = 1;
%  end
%  maxnull = round(256*(g.signif(2)-g.limits(5))/range);
%  if maxnull>256
%    maxnull = 256;
%  end
%  nullrange = minnull:maxnull;
%  cc(nullrange,:) = repmat(cc(128,:),length(nullrange),1);
%end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot tfdata image for specified channel or selchans std()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axis off;
colormap(cc);
curax = gca; % current plot axes to plot into
if tfpoints ~= 0
    plotdim = max(1+floor(tfpoints/2),4); % number of topoplots on top of image
    imgax = sbplot(plotdim,plotdim,[plotdim*(plotdim-1)+1,2*plotdim-1],'ax',curax);
else
    imgax = sbplot(1,1,1,'ax',curax);
end;
tftimes = mmidx(1):mmidx(2);
tffreqs = mmidx(3):mmidx(4);
if g.showchan>0 % -> image showchan data
    tfave = tfdata(tffreqs, tftimes,g.showchan);
else % g.showchan==0 -> image std() of selchans
    tfdat   = tfdata(tffreqs,tftimes,g.selchans,:);

    % average across electrodes
    fprintf('Applying RMS across channels (mask for at least %d non-zeros values at each time/freq)\n', g.sigthresh(1));
    tfdat = avedata(tfdat, 3, g.sigthresh(1), g.mode);

    % if several subject, first (RMS) averaging across subjects
    if size(tfdata,4) > 1
        fprintf('Applying RMS across subjects (mask for at least %d non-zeros values at each time/freq)\n', g.sigthresh(2));
        tfdat = avedata(tfdat, 4, g.sigthresh(2), g.mode);
    end;
    tfave = tfdat;
    
    if defaultlim
        g.limits(6) = max(max(abs(tfave)));
        g.limits(5) = -g.limits(6); % make symmetrical
    end;
end

if strcmpi(g.logfreq, 'on'), 
    logimagesc(times(tftimes),freqs(tffreqs),tfave);
    axis([g.limits(1) g.limits(2) log(g.limits(3)), log(g.limits(4))]);
else                            
    imagesc(times(tftimes),freqs(tffreqs),tfave);
    axis([g.limits(1:4)]);
end;
caxis([g.limits(5:6)]);
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title and vertical lines
%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(imgax)
xl=xlabel('Time (ms)');
set(xl,'fontsize',16);
set(gca,'yaxislocation','left')
if g.showchan>0
   % tl=title(['Channel ',int2str(g.showchan)]);
   % set(tl,'fontsize',14);
else
    if strcmpi(g.mode, 'rms')
        tl=title(['Signed channel rms']);
    else
        tl=title(['Signed channel average']);
    end;
  set(tl,'fontsize',14);
end

yl=ylabel(['Frequency (Hz)']);
set(yl,'fontsize',16);

set(gca,'fontsize',14)
set(gca,'ydir','normal');

for indtime = g.vert
    if strcmpi(g.logfreq, 'on'), 
        plot([indtime indtime],log([freqs(1) freqs(end)]),[LINECOLOR ':'],'linewidth',ZEROLINEWIDTH);
    else
        plot([indtime indtime],[freqs(1) freqs(end)],[LINECOLOR ':'],'linewidth',ZEROLINEWIDTH);
    end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot topoplot maps at specified timefreqs points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(g.timefreqs)
    wholeax  = sbplot(1,1,1,'ax',curax);
    topoaxes = zeros(1,tfpoints);
    for n=1:tfpoints
        if n<=plotdim
            topoaxes(n)=sbplot(plotdim,plotdim,n,'ax',curax);
        else
            topoaxes(n)=sbplot(plotdim,plotdim,plotdim*(n+1-plotdim),'ax',curax);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot connecting lines using changeunits()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmpi(g.logfreq, 'off')
            from = changeunits([g.timefreqs(n,:)],imgax,wholeax);
        else  
            from = changeunits([g.timefreqs(n,1) log(g.timefreqs(n,2))],imgax,wholeax);        
        end;
        to   = changeunits([0.5,0.5],topoaxes(n),wholeax);
        axes(wholeax);
        plot([from(1) to(1)],[from(2) to(2)],LINECOLOR,'linewidth',LINEWIDTH);
        hold on
        mk=plot(from(1),from(2),[LINECOLOR 'o'],'markersize',9);
        set(mk,'markerfacecolor',LINECOLOR);
        axis([0 1 0 1]);
        axis off;
    end
    
    endcaxis = 0;
    for n=1:tfpoints
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot scalp map using topoplot()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(topoaxes(n));
        scalpmap = squeeze(tfdataori(tfpidx(n,2),tfpidx(n,1),g.selchans));   

        %topoplot(scalpmap,g.chanlocs,'maplimits',[g.limits(5) g.limits(6)],...
        %            'electrodes','on');

        if ~isempty(varargin)
            topoplot(scalpmap,g.chanlocs,'electrodes','on', varargin{:}); 
        else
            topoplot(scalpmap,g.chanlocs,'electrodes','on'); 
        end;
        % 'interlimits','electrodes')
        axis square;
        hold on
        tl=title([int2str(g.timefreqs(n,1)),' ms, ',int2str(g.timefreqs(n,2)),' Hz']);
        set(tl,'fontsize',13);
        endcaxis = max(endcaxis,max(abs(caxis)));
        %caxis([g.limits(5:6)]);
    end;
    for n=1:tfpoints
        axes(topoaxes(n));
        caxis([-endcaxis endcaxis]);
        if n==tfpoints % & (mod(tfpoints,2)~=0) % image color bar by last map
            cb=cbar;
            pos = get(cb,'position');
            set(cb,'position',[pos(1:2) 0.023 pos(4)]);
        end
        drawnow
    end
end;

if g.showchan>0 & ~isempty(g.chanlocs)
     sbplot(4,4,1,'ax',imgax);
     topoplot(g.showchan,g.chanlocs,'electrodes','off', ...
                  'style', 'blank', 'emarkersize1chan', 10 );
     axis('square');
end
axcopy;


function tfdat = avedata(tfdat, dim, thresh, mode)
    tfsign  = sign(mean(tfdat,dim));
    tfmask  = sum(tfdat ~= 0,dim) >= thresh;
    if strcmpi(mode, 'rms')
        tfdat   = tfmask.*tfsign.*sqrt(mean(tfdat.*tfdat,dim)); % std of all channels
    else
        tfdat   = tfmask.*mean(tfdat,dim); % std of all channels
    end;
    
function [tfdatnew, times, freqs] = magnifytwice(tfdat, times, freqs);
    indicetimes = [floor(1:0.5:size(tfdat,1)) size(tfdat,1)];
    indicefreqs = [floor(1:0.5:size(tfdat,2)) size(tfdat,2)];
    tfdatnew = tfdat(indicetimes, indicefreqs, :, :);
    times = linspace(times(1), times(end), size(tfdat,2)*2);
    freqs = linspace(freqs(1), freqs(end), size(tfdat,1)*2);
    
    % smoothing
    gauss2 = gauss2d(3,3); 
    for S = 1:size(tfdat,4)
        for elec = 1:size(tfdat,3)
            tfdatnew(:,:,elec,S) = conv2(tfdatnew(:,:,elec,S), gauss2, 'same');
        end;
    end;

    %tfdatnew = convn(tfdatnew, gauss2, 'same'); % is equivalent to the loop for slowlier

