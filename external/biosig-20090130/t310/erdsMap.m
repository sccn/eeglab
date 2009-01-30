function erds = erdsMap(X, frames, tlimits, Fs, varwin, varargin)
% Plots ERDS maps based on either FFTs or wavelets.
%
% This function calculates time-frequency maps (ERDS maps) using either wavelets
% or FFT analysis. The result is stored in the output variable, which can be
% used by other functions to plot the map, for example.
%
% Usage:
%   erds = erdsmap(X, frames, tlimits, Fs);
%
% Input parameters:
%   X       ... Single channel data vector <1 x frames*ntrials>
%   frames  ... Frames per trial
%   tlimits ... Epoch time limits (ms) [mintime maxtime]
%   Fs      ... Sampling rate (Hz)
%   cycles  ... =0: Use FFTs (with constant window length)
%               >0: Number of cycles in each analysis wavelet
%
% Optional input parameters:
%   'detret'   ... Detrend data in time ['on'|'off']
%   'detrep'   ... Detrend data across trials ['on'|'off']
%   'winsize'  ... If cycles = 0: Data subwindow length
%                    If cycles > 0: Longest window length to use, determines
%                                   the lowest output frequency       {frames/8}
%   'timesout' ... Number of output times                                {200}
%   'padratio' ... FFT length/winframes (2^k)                              {2}
%                    Multiplies the number of output frequencies by dividing
%                    their spacing; when cycles = 0, frequency spacing is
%                    (low_freq/padratio)
%   'maxfreq'  ... Maximum frequency (Hz) to plot                         {40}
%   'baseline' ... Spectral baseline window center end-time (in ms)        {0}
%   'powbase'  ... Baseline power spectrum to normalize the data.
%   'alpha'    ... If non-zero, compute bootstrap significance level    {0.05}
%   'naccu'    ... Number of bootstrap replications to accumulate        {200}
%
% Output parameter:
%   erds ... Structure containing information about the ERDS map

% Copyright by Clemens Brunner, based on timef of EEGLAB
% $Revision: 1.1 $ $Date: 2009-01-30 06:04:44 $
% E-Mail: clemens.brunner@tugraz.at

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

if (nargin < 1)
    error('No input arguments specified.');
    return;
end;

if (min(size(X)) ~= 1 | length(X) < 2)
    error('Data must be a row or column vector.');
end;

if (nargin < 2 | isempty(frames) | isnan(frames))
    frames = NaN;  % Frames per trial
elseif (~isnumeric(frames) | length(frames) ~= 1 | frames ~= round(frames))
    error('Value of frames must be an integer.');
elseif (frames <= 0)
    error('Value of frames must be positive.');
elseif (rem(length(X), frames) ~= 0)
    error('Length of data vector must be divisible by frames.');
end;
if (isnan(frames) | isempty(frames))
    frames = length(X);
end;

if (nargin < 3 | isnan(tlimits) | isempty(tlimits))
    tlimits = NaN;  % Time range of g.frames (ms)
elseif (~isnumeric(tlimits) | sum(size(tlimits)) ~= 3)
    error('Value of tlimits must be a vector containing two numbers.');
elseif (tlimits(1) >= tlimits(2))
    error('tlimits interval must be ascending.');
end;

if (nargin < 4)
    Fs = 250;
elseif (~isnumeric(Fs) | length(Fs) ~= 1)
    error('Value of srate must be a number.');
elseif (Fs <= 0)
    error('Value of srate must be positive.');
end;

if (isnan(tlimits) | isempty(tlimits))
    hlim = 1000*frames/Fs;  % Fit default tlimits to srate and frames
    tlimits = [0 hlim];
end;

framesdiff = frames - Fs*(tlimits(2)-tlimits(1))/1000;
if (abs(framesdiff) > 1)
    error('Given time limits, frames and sampling rate are incompatible.');
elseif (framesdiff ~= 0)
    tlimits(1) = tlimits(1) - 0.5*framesdiff*1000/Fs;
    tlimits(2) = tlimits(2) + 0.5*framesdiff*1000/Fs;
    fprintf('Adjusted time limits to [%.1f,%.1f] ms.\n',tlimits(1),tlimits(2));
end;

if (nargin < 5)
    varwin = 0;
elseif (~isnumeric(varwin) | length(varwin) > 2)
    error('Value of cycles must be a number.');
elseif (varwin < 0)
    error('Value of cycles must be zero or positive.');
end;

if (~isempty(varargin))
    try
        g = struct(varargin{:});
    catch
        error('Argument error in the {''param'', value} sequence');
    end;
end;

tlimits(1) = tlimits(1) - 0.5 * 1000;
tlimits(2) = tlimits(2) + 0.5 * 1000;
g.tlimits = tlimits;

g.frames = frames;
g.srate = Fs;
g.cycles = varwin(1);
if (length(varwin) > 1)
    g.cyclesfact = varwin(2);
else
    g.cyclesfact = 1;
end;

try g.title;     catch g.title = ''; end;
try g.winsize;   catch g.winsize = max(pow2(nextpow2(g.frames)-3),4); end;
try g.pad;       catch g.pad = max(pow2(nextpow2(g.winsize)),4); end;
try g.timesout;  catch g.timesout = 200; end;
try g.padratio;  catch g.padratio = 2; end;
try g.maxfreq;   catch g.maxfreq = 40; end;
try g.topovec;   catch g.topovec = []; end;
try g.alpha;     catch g.alpha = 0.05; end;
try g.powbase;   catch g.powbase = NaN; end;
try g.pboot;     catch g.pboot = NaN; end;
try g.rboot;     catch g.rboot = NaN; end;
try g.plotersp;  catch g.plotersp = 'on'; end;
try g.plotitc;   catch g.plotitc = 'off'; end;
try g.detrep;    catch g.detrep = 'off'; end;
try g.detret;    catch g.detret = 'off'; end;
try g.baseline;  catch g.baseline = 0; end;
try g.baseboot;  catch g.baseboot = 1; end;
try g.linewidth; catch g.linewidth = 2; end;
try g.naccu;     catch g.naccu = 200; end;
try g.mtaper;    catch g.mtaper = []; end;
try g.vert;      catch g.vert = []; end;
try g.type;      catch g.type = 'phasecoher'; end;
try g.phsamp;    catch g.phsamp = 'off'; end;
try g.plotphase; catch g.plotphase = 'on'; end;
try g.itcmax;    catch g.itcmax = []; end;
try g.erspmax;   catch g.erspmax = []; end;
try g.verbose;   catch g.verbose = 'on'; end;
try g.chaninfo;  catch g.chaninfo = []; end;
try g.hzdir;     catch g.hzdir = 'up'; end;
try g.cbar;      catch g.cbar = 'on'; end;
try g.margersp;  catch g.margersp = 'on'; end;
try g.meanspec;  catch g.meanspec = 'on'; end;

% testing arguments consistency
if strcmp(g.hzdir,'up')
    g.hzdir = 'normal';
elseif strcmp(g.hzdir,'down')
    g.hzdir = 'reverse';
else
    error('Unknown ''hzdir'' value - not ''up'' or ''down''.');
end;

if (~isnumeric(g.winsize) || length(g.winsize) ~= 1 || ...
    g.winsize ~= round(g.winsize))
    error('Value of winsize must be an integer number.');
elseif (g.winsize <= 0)
    error('Value of winsize must be positive.');
elseif (g.cycles == 0 && pow2(nextpow2(g.winsize)) ~= g.winsize)
    error('Value of winsize must be an integer power of 2 [1,2,4,8,16,...].');
elseif (g.winsize > g.frames)
    error('Value of winsize must be less than frames per epoch.');
end;

if (~isnumeric(g.timesout) || length(g.timesout) ~=1 || ...
    g.timesout ~= round(g.timesout))
    error('Value of timesout must be an integer number.');
elseif (g.timesout <= 0)
    error('Value of timesout must be positive.');
end;
if (g.timesout > g.frames-g.winsize)
    g.timesout = g.frames-g.winsize;
    disp(['Value of timesout must be <= frames-wins, timeout adjusted to ' ...
         int2str(g.timesout)]);
end;

if (~isnumeric(g.padratio) || length(g.padratio) ~= 1 || ...
    g.padratio ~= round(g.padratio))
    error('Value of padratio must be an integer.');
elseif (g.padratio <= 0)
    error('Value of padratio must be positive.');
elseif (pow2(nextpow2(g.padratio)) ~= g.padratio)
    error('Value of padratio must be an integer power of 2 [1,2,4,8,16,...].');
end;

if (~isnumeric(g.maxfreq) || length(g.maxfreq) ~= 1)
    error('Value of maxfreq must be a number.');
elseif (g.maxfreq <= 0)
    error('Value of maxfreq must be positive.');
elseif (g.maxfreq > Fs/2)
    myprintf(g.verbose,['Warning: value of maxfreq reduced to Nyquist rate' ...
        ' (%3.2f)\n\n'], Fs/2);
    g.maxfreq = Fs/2;
end;

if (~isnumeric(g.alpha) || length(g.alpha) ~= 1)
    error('Value of alpha must be a number.');
elseif (round(g.naccu*g.alpha) < 2)
    myprintf(g.verbose,'Value of g.alpha is out of range [%g,0.5]\n',2/g.naccu);
    g.naccu = round(2/g.alpha);
    myprintf(g.verbose,'Increasing bootstrap iterations to %d\n',g.naccu);
end;
if g.alpha>0.5 | g.alpha<=0
    error('Value of g.alpha is out of the allowed range (0.00,0.5).');
end;

if ~isnumeric(g.vert)
    error('vertical line(s) option must be a vector');
else
    if min(g.vert) < g.tlimits(1) | max(g.vert) > g.tlimits(2)
        error('vertical line(s) time out-of-bound');
    end;
end;

if ~isnan (g.rboot)
    if size(g.rboot) == [1,1]
        if g.cycles == 0
            g.rboot = g.rboot*ones(g.winsize*g.padratio/2);
        end;
    end;
end;

if ~isempty(g.mtaper) % mutitaper, inspired from Bijan Pesaran matlab function
    if length(g.mtaper) < 3
        if g.mtaper(1) * g.mtaper(2) < 1
            error('mtaper 2 first arguments'' product must be higher than 1');
        end;
        if length(g.mtaper) == 2
            g.mtaper(3) = floor( 2*g.mtaper(2)*g.mtaper(1) - 1);
        end
        if length(g.mtaper) == 3
            if g.mtaper(3) > 2 * g.mtaper(1) * g.mtaper(2) -1
                error('mtaper number too high (maximum (2*N*W-1))');
            end;
        end;
        disp(['Using ' num2str(g.mtaper(3)) ' tapers.']);
        NW = g.mtaper(1)*g.mtaper(2);   % product NW
        N  = g.mtaper(1)*g.srate;
        [e,v] = dpss(N, NW, 'calc');
        e=e(:,1:g.mtaper(3));
        g.alltapers = e;
    else
        g.alltapers = g.mtaper;
        disp('mtaper not [N W] or [N W K]; considering raw taper matrix');
    end;

    g.winsize = size(g.alltapers, 1);
    g.pad = max(pow2(nextpow2(g.winsize)),256); % pad*nextpow

    %nfk = floor([0 g.maxfreq]./g.srate.*g.pad); % not used any more
    %g.padratio = 2*nfk(2)/g.winsize;

    g.padratio = g.pad/g.winsize;

    %compute number of frequencies
    %nf = max(256, g.pad*2^nextpow2(g.winsize+1));
    %nfk = floor([0 g.maxfreq]./g.srate.*nf);
    %freqs = linspace( 0, g.maxfreq, diff(nfk));

end;

switch lower(g.plotphase)
    case { 'on', 'off' }, ;
    otherwise error('plotphase must be either on or off');
end;
switch lower(g.plotersp)
    case { 'on', 'off' }, ;
    otherwise error('plotersp must be either on or off');
end;
switch lower(g.plotitc)
    case { 'on', 'off' }, ;
    otherwise error('plotitc must be either on or off');
end;
switch lower(g.detrep)
    case { 'on', 'off' }, ;
    otherwise error('detrep must be either on or off');
end;
switch lower(g.detret)
    case { 'on', 'off' }, ;
    otherwise error('detret must be either on or off');
end;
switch lower(g.phsamp)
    case { 'on', 'off' }, ;
    otherwise error('phsamp must be either on or off');
end;
if ~isnumeric(g.linewidth)
    error('linewidth must be numeric');
end;
if ~isnumeric(g.naccu)
    error('naccu must be numeric');
end;
if ~isnumeric(g.baseline)
    error('baseline must be numeric');
end;
if isempty(g.baseline)
    g.baseline = 0;
end;
switch g.baseboot
    case {0,1}, ;
    otherwise, error('baseboot must be 0 or 1');
end;
switch g.type
    case { 'coher', 'phasecoher', 'phasecoher2' },;
    otherwise error('Type must be either ''coher'' or ''phasecoher''');
end;
switch g.cbar
    case { 'on', 'off' }, ;
    otherwise error('cbar must be either on or off');
end;
switch g.margersp
    case { 'on', 'off' }, ;
    otherwise error('margersp must be either on or off');
end;
switch g.meanspec
    case { 'on', 'off' }, ;
    otherwise error('meanspec must be either on or off');
end;

if (g.cycles == 0) %%%%%%%%%%%%%% constant window-length FFTs %%%%%%%%%%%%%%%%%%
    freqs = linspace(0, g.srate/2, g.padratio*g.winsize/2+1);
    freqs = freqs(2:end);
    win = hanning(g.winsize);

    P  = zeros(g.padratio*g.winsize/2,g.timesout); % summed power
    PP = zeros(g.padratio*g.winsize/2,g.timesout); % power
    R  = zeros(g.padratio*g.winsize/2,g.timesout); % mean coherence
    RR = zeros(g.padratio*g.winsize/2,g.timesout); % (coherence)
    Pboot = zeros(g.padratio*g.winsize/2,g.naccu); % summed bootstrap power
    Rboot = zeros(g.padratio*g.winsize/2,g.naccu); % summed bootstrap coher
    Rn = zeros(1,g.timesout);
    Rbn = 0;

    switch g.type
        case { 'coher' 'phasecoher2' },
            cumulX = zeros(g.padratio*g.winsize/2,g.timesout);
            cumulXboot = zeros(g.padratio*g.winsize/2,g.naccu);
        case 'phasecoher'
            switch g.phsamp
                case 'on'
                    cumulX = zeros(g.padratio*g.winsize/2,g.timesout);
            end;
    end;

else % %%%%%%%%%%%%%%%%%% cycles>0, Constant-Q (wavelet) DFTs %%%%%%%%%%%%%%%%%%

    freqs = g.srate*g.cycles/g.winsize*[2:2/g.padratio:g.winsize]/2;
    dispf = find(freqs <= g.maxfreq);
    freqs = freqs(dispf);

    win = dftfilt(g.winsize,g.maxfreq/g.srate,g.cycles,g.padratio,g.cyclesfact);
    P = zeros(size(win,2),g.timesout);       % summed power
    R = zeros(size(win,2),g.timesout);       % mean coherence
    PP = repmat(NaN,size(win,2),g.timesout); % initialize with NaN
    RR = repmat(NaN,size(win,2),g.timesout); % initialize with NaN
    Pboot = zeros(size(win,2),g.naccu);  % summed bootstrap power
    Rboot = zeros(size(win,2),g.naccu);  % summed bootstrap coher
    Rn = zeros(1,g.timesout);
    Rbn = 0;

    switch g.type
        case { 'coher' 'phasecoher2' },
            cumulX = zeros(size(win,2),g.timesout);
            cumulXboot = zeros(size(win,2),g.naccu);
        case 'phasecoher'
            switch g.phsamp
                case 'on'
                    cumulX = zeros(size(win,2),g.timesout);
            end;
    end;
end;

switch g.phsamp
    case 'on'
        PA = zeros(size(P,1),size(P,1),g.timesout); % NB: (freqs,freqs,times)
end;                                                %       phs   amp

wintime = 1000/g.srate*(g.winsize/2); % (1000/g.srate)*(g.winsize/2);
times = [g.tlimits(1)+wintime:(g.tlimits(2)-g.tlimits(1)-2*wintime)/...
        (g.timesout-1):g.tlimits(2)-wintime];
ERPtimes = [g.tlimits(1):(g.tlimits(2)-g.tlimits(1))/(g.frames-1):...
           g.tlimits(2)+0.000001];
ERPindices = [];
for ti=times
    [tmp indx] = min(abs(ERPtimes-ti));
    ERPindices  = [ERPindices indx];
end;
ERPtimes = ERPtimes(ERPindices); % subset of ERP frames on t/f window centers

if (g.baseline ~= 0 & ~isempty(g.baseline))
    if ~isempty(find(times < g.baseline(2) & times > g.baseline(1)))
        baseln = find(times < g.baseline(2) & times > g.baseline(1));
    end;
else
    baseln = 1:length(times); % use all times as baseline
end;
if ~isnan(g.alpha) & length(baseln)==0
    return;
end;
dispf = find(freqs <= g.maxfreq);
stp = (g.frames-g.winsize)/(g.timesout-1);

trials = length(X)/g.frames;
baselength = length(baseln);

% detrend over epochs (trials) if requested
switch g.detrep
    case 'on'
        X = reshape(X, g.frames, length(X)/g.frames);
        X = X - mean(X,2)*ones(1, length(X(:))/g.frames);
        X = X(:)';
end;

for i=1:trials

    % ERP = blockave(X,g.frames); % compute the ERP trial average

    Wn = zeros(1,g.timesout);
    for j=1:g.timesout,
        tmpX = X([1:g.winsize]+floor((j-1)*stp)+(i-1)*g.frames);
        % pull out data g.frames
        tmpX = tmpX - mean(tmpX); % remove the mean for that window
        switch g.detret, case 'on', tmpX = detrend(tmpX); end;
        if ~any(isnan(tmpX))
            if (g.cycles == 0) % FFT
                if ~isempty(g.mtaper)   % apply multitaper (no hanning window)
                    tmpXMT = fft(g.alltapers .* ...
                        (tmpX(:) * ones(1,size(g.alltapers,2))), g.pad);
                    %tmpXMT = tmpXMT(nfk(1)+1:nfk(2),:);
                    tmpXMT = tmpXMT(2:g.padratio*g.winsize/2+1,:);
                    PP(:,j) = mean(abs(tmpXMT).^2, 2);
                    % power; can also ponderate multitaper by their eigenvalues
                    tmpX = win .* tmpX(:);
                    tmpX = fft(tmpX, g.pad);
                    tmpX = tmpX(2:g.padratio*g.winsize/2+1);
                else
                    tmpX = win .* tmpX(:);
                    tmpX = fft(tmpX,g.padratio*g.winsize);
                    tmpX = tmpX(2:g.padratio*g.winsize/2+1);
                    PP(:,j) = abs(tmpX).^2; % power
                end;
            else % wavelet
                if ~isempty(g.mtaper)  % apply multitaper
                    tmpXMT = g.alltapers .* (tmpX(:) * ...
                             ones(1,size(g.alltapers,2)));
                    tmpXMT = transpose(win) * tmpXMT;
                    PP(:,j) = mean(abs(tmpXMT).^2, 2); % power
                    tmpX = transpose(win) * tmpX(:);
                else
                    tmpX = transpose(win) * tmpX(:);
                    PP(:,j) = abs(tmpX).^2; % power
                end;
            end;

            if abs(tmpX) < 1e-8
                RR(:,j) = zeros(size(RR(:,j)));
            else
                switch g.type
                    case { 'coher' },
                        RR(:,j) = tmpX;
                        cumulX(:,j) = cumulX(:,j)+abs(tmpX).^2;
                    case { 'phasecoher2' },
                        RR(:,j) = tmpX;
                        cumulX(:,j) = cumulX(:,j)+abs(tmpX);
                    case 'phasecoher',
                        RR(:,j) = tmpX ./ abs(tmpX); % normalized cross-spectral
                        switch g.phsamp
                            case 'on'
                                cumulX(:,j) = cumulX(:,j)+abs(tmpX);
                        end;
                end;
            end;
            Wn(j) = 1;
        end;

        switch g.phsamp
            case 'on' % PA (freq x freq x time)
                PA(:,:,j) = PA(:,:,j)  + (tmpX ./ abs(tmpX)) * ((PP(:,j)))';
                % cross-product: unit phase (column) times amplitude (row)
        end;
    end; % window

    if ~isnan(g.alpha)  % save surrogate data for bootstrap analysis
        j = 1;
        goodbasewins = find(Wn==1);
        if g.baseboot  % use baseline windows only
            goodbasewins = find(goodbasewins<=baselength);
        end;
        ngdbasewins = length(goodbasewins);
        if ngdbasewins>1
            while j <= g.naccu
                i = ceil(rand*ngdbasewins);
                i = goodbasewins(i);
                Pboot(:,j) = Pboot(:,j) + PP(:,i);
                Rboot(:,j) = Rboot(:,j) + RR(:,i);
                switch g.type
                    case 'coher'
                        cumulXboot(:,j) = cumulXboot(:,j)+abs(tmpX).^2;
                    case 'phasecoher2'
                        cumulXboot(:,j) = cumulXboot(:,j)+abs(tmpX);
                end;
                j = j+1;
            end;
            Rbn = Rbn + 1;
        end;
    end;  % bootstrap

    Wn = find(Wn>0);
    if length(Wn)>0
        P(:,Wn) = P(:,Wn) + PP(:,Wn); % add non-NaN windows
        R(:,Wn) = R(:,Wn) + RR(:,Wn);
        Rn(Wn) = Rn(Wn) + ones(1,length(Wn)); % count number of addends
    end;
end; % trial

% if coherence, perform the division
switch g.type
    case 'coher',
        R = R ./ ( sqrt( trials*cumulX ) );
        if ~isnan(g.alpha)
            Rboot = Rboot ./ ( sqrt( trials*cumulXboot ) );
        end;
    case 'phasecoher2',
        R = R ./ ( cumulX );
        if ~isnan(g.alpha)
            Rboot = Rboot ./ cumulXboot;
        end;
    case 'phasecoher',
        R = R ./ (ones(size(R,1),1)*Rn);
end;

switch g.phsamp
    case 'on'
        tmpcx(1,:,:) = cumulX; % allow ./ below
        for j=1:g.timesout
            PA(:,:,j) = PA(:,:,j) ./ repmat(PP(:,j)', [size(PP,1) 1]);
        end;
end;

if min(Rn) < 1
    myprintf(g.verbose,'No valid timef estimates for windows %s of %d.\n',...
        int2str(find(Rn==0)),length(Rn));
    Rn(find(Rn<1))==1;
    return;
end;

P = P ./ (ones(size(P,1),1) * Rn);

if isnan(g.powbase)
    mbase = mean(P(:,baseln),2)';
else
    mbase = g.powbase;
end;

if ~isnan( mbase )
    %P = 10 * (log10(P) - repmat(log10(mbase(1:size(P,1)))',[1 g.timesout]));
    P = P./repmat(mbase',1,size(P,2)) - 1;
else
    %P = 10 * log10(P);
end;

Rsign = sign(imag(R));
if nargout > 7
    for lp = 1:size(R,1)
        Rphase(lp,:) = rem(angle(R(lp,:)),2*pi);
    end;
    Rphase(find(Rphase>pi))  = 2*pi-Rphase(find(Rphase>pi));
    Rphase(find(Rphase<-pi)) = -2*pi-Rphase(find(Rphase<-pi));
end;

R = abs(R); % convert coherence vector to magnitude

if ~isnan(g.alpha) % if bootstrap analysis included . . .
    if Rbn>0
        i = round(g.naccu*g.alpha);
        if isnan(g.pboot)
            Pboot = Pboot / Rbn; % normalize
            if ~isnan( g.baseline )
                %Pboot = 10*(log10(Pboot) - repmat(log10(mbase)',[1 g.naccu]));
                Pboot = Pboot./repmat(mbase',1,size(P,2)) - 1;
            else
                %Pboot = 10 * log10(Pboot);
            end;
            Pboot = sort(Pboot');
            Pboot = [mean(Pboot(1:i,:)) ; mean(Pboot(g.naccu-i+1:g.naccu,:))];
        else
            Pboot = g.pboot;
        end;

        if isnan(g.rboot)
            Rboot = abs(Rboot) / Rbn;
            Rboot = sort(Rboot');
            Rboot = mean(Rboot(g.naccu-i+1:g.naccu,:));
        else
            Rboot = g.rboot;
        end;
    else
        myprintf(g.verbose,'No valid bootstrap trials...!\n');
    end;
end;

PP = P;
if ~isnan(g.alpha)
    PP(find((PP > repmat(Pboot(1,:)',[1 g.timesout])) ...
            & (PP < repmat(Pboot(2,:)',[1 g.timesout])))) = 0;
end;

erds.P = P;
erds.PP = PP;
erds.sig = Pboot;
erds.dispf = dispf;
erds.freqs = freqs;
erds.times = times;
erds.baseline = g.baseline;
erds.Pref = mbase;


% symmetric Hanning tapering function
function w = hanning(n)
if ~rem(n,2)
    w = .5*(1 - cos(2*pi*(1:n/2)'/(n+1)));
    w = [w; w(end:-1:1)];
else
    w = .5*(1 - cos(2*pi*(1:(n+1)/2)'/(n+1)));
    w = [w; w(end-1:-1:1)];
end

function myprintf(verbose, varargin)
if strcmpi(verbose, 'on')
    fprintf(varargin{:});
end;
