% pac_cont() - compute phase-amplitude coupling (power of first input
%         correlation with phase of second). There is no graphical output
%         to this function.
%
% Usage:
%   >> pac_cont(x,y,srate);
%   >> [pac timesout pvals] = pac_cont(x,y,srate,'key1', 'val1', 'key2', val2' ...);
%
% Inputs:
%    x       = [float array] 1-D data vector of size (1xtimes)
%    y       = [float array] 1-D data vector of size (1xtimes)
%    srate   = data sampling rate (Hz)
%
% Optional time information intputs:
%   'winsize'   = If cycles==0: data subwindow length (fastest, 2^n<frames);
%                 If cycles >0: *longest* window length to use. This
%                 determines the lowest output frequency. Note that this
%                 parameter is overwritten if the minimum frequency has been set
%                 manually and requires a longer time window {~frames/8}
%   'ntimesout' = Number of output times (int<frames-winframes). Enter a
%                 negative value [-S] to subsample original time by S.
%   'timesout'  = Enter an array to obtain spectral decomposition at
%                 specific time values (note: algorithm find closest time
%                 point in data and this might result in an unevenly spaced
%                 time array). Overwrite 'ntimesout'. {def: automatic}
%   'tlimits'   = [min max] time limits in ms.
%
% Optional PAC inputs:
%   'method'    = ['modulation'|'plv'|'corr'|'glm'] (see reference).
%   'freqphase' = [min max] frequency limits. Default [minfreq 50],
%                 minfreq being determined by the number of data points,
%                 cycles and sampling frequency. Use 0 for minimum frequency
%                 to compute default minfreq. You may also enter an
%                 array of frequencies for the spectral decomposition
%                 (for FFT, closest computed frequency will be returned; use
%                 'padratio' to change FFT freq. resolution).
%   'freqamp'   = [float array] array of frequencies for the second
%                 argument. 'freqs' is used for the first argument.
%                 By default it is the same as 'freqs'.
%   'filterfunc' = ['eegfilt'|'iirfilt'|'eegfftfilt'] filtering function.
%                 Default is iirfilt. Warning, filtering may dramatically
%                 affect the result. With the 'corr' method, make sure you
%                 have a large window size because each window is filtered
%                 independently.
%   'filterphase' = @f_handle. Function handle to filter the data for the 
%                 phase information. For example, @(x)iirfilt(x, 1000, 2,
%                 20). Note that 'freqphase' is ignore in this case.
%   'filteramp' = @f_handle. Function handle to filter the data for the 
%                 amplitude information. Note that 'freqamp' is ignore in 
%                 this case.
%
% Inputs for statistics:
%   'alpha'     = [float] p-value threshold. Default is none (no statistics).
%   'mcorrect'  = ['none'|'fdr'] method to correct for multiple comparison.
%                 Default is 'none'.
%   'baseline'  = [min max] baseline period for the Null distribution. Default 
%                 is the whole data range. Note that this option is ignored 
%                 for instanstaneous statistics.
%   'instantstat' = ['on'|'off'] performs statistics for each time window
%                 independently. Default is 'off'.
%   'naccu'     = [integer] number of accumulations for surrogate
%                 statistics.
%   'statlim'   = ['parametric'|'surrogate'] use a parametric methods to
%                 asseess the limit of the surrogate distribution or use
%                 the tail of the distribution ('surrogate' method)
%
% Other inputs:
%   'title'     = [string] figure title. Default is none.
%   'vert'      = [float array] array of time value for which to plot
%                 vertical lines. Default is none.
%
% Outputs:
%  pac      = Phase-amplitude coupling values.
%  timesout = vector of time indices
%  pvals    = Associated p-values
%   
% Author: Arnaud Delorme and Makoto Miyakoshi, SCCN/INC, UCSD 2012-
%
% References:
% Methods used here are introduced and compared in:
% Penny, Duzel, Miller, Ojemann. (20089). Testing for Nested Oscilations.
% J Neuro Methods. 174:50-61
%
% Modulation Index is defined in:
% Canolty, Edwards, Dalal, Soltani, Nagarajan, Kirsch, et al. (2006). Modulation index is defined in High Gamma Power Is Phase-Locked to Theta 
% Oscillations in Human Neocortex. Science. 313:1626-8.
%
% PLV (Phase locking value) is defined in:
% Lachaux, Rodriguez, Martiniere, Varela. (1999). Measuring phase synchrony
% in brain signal. Hum Brain Mapp. 8:194-208.
%
% corr (correlation) method is defined in:
% Brunce, Eckhorn. (2004). Task-related coupling from high- to
% low-frequency signals among visual cortical areas in human subdural
% recordings. Int J Psychophysiol. 51:97-116.
%
% glm (general linear model) is defined in 
% Penny, Duzel, Miller, Ojemann. (20089). Testing for Nested Oscilations.
% J Neuro Methods. 174:50-61

% Copyright (C) 2012 Arnaud Delorme, UCSD
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

function [m_raw pvals indexout] = pac_cont(X, Y, srate, varargin);

if nargin < 1
    help pac_cont;
    return;
end;

% deal with 3-D inputs
% --------------------
if ndims(X) == 3 || ndims(Y) == 3, error('Cannot process 3-D input'); end;
if size(X,1) > 1, X = X'; end;
if size(Y,1) > 1, Y = Y'; end;
if size(X,1) ~= 1 || size(Y,1) ~= 1, error('Cannot only process vector input'); end;
frame = size(X,2);
pvals = [];

g = finputcheck(varargin, ...
    { ...
    'alpha'         'real'     [0 0.2]                [];
    'baseline'      'float'    []                       [];
    'freqphase'     'real'     [0 Inf]                  [0 srate/2];
    'freqamp'       'real'     [0 Inf]                  [];
    'mcorrect'      'string'   { 'none' 'fdr' }         'none';
    'method'        'string'   { 'plv' 'modulation' 'glm' 'corr' } 'modulation';
    'naccu'         'integer'  [1 Inf]                   250;
    'instantstat'   'string'   {'on','off'}              'off';
    'newfig'        'string'   {'on','off'}              'on';
    'nofig'         'string'   {'on','off'}              'off';
    'statlim'       'string'   {'surrogate','parametric'}  'parametric';
    'timesout'      'real'     []                        []; ...
    'filterfunc'    'string'   { 'eegfilt' 'iirfilt' 'eegfiltfft' }   'eegfiltfft'; ...
    'filterphase'   ''         {}                        [];
    'filteramp'     ''         {}                        [];
    'ntimesout'     'integer'  []                        200; ...
    'tlimits'       'real'     []                        [0 frame/srate];
    'title'         'string'   []                        '';
    'vert'          'real'     []                    [];
    'winsize'       'integer'  [0 Inf]                   max(pow2(nextpow2(frame)-3),4) }, 'pac');

if isstr(g), error(g); end;

if ~isempty(g.filterphase)
     x_freqphase = feval(g.filterphase, X(:)');
else x_freqphase = feval(g.filterfunc, X(:)', srate, g.freqphase(1), g.freqphase(end));
end;
if ~isempty(g.filteramp)
     x_freqamp   = feval(g.filteramp, Y(:)');
else x_freqamp   = feval(g.filterfunc, Y(:)', srate, g.freqamp(  1), g.freqamp( end));
end;
z_phasedata = hilbert(x_freqphase);
z_ampdata   = hilbert(x_freqamp);
phase       = angle(z_phasedata);
amplitude   = abs(  z_ampdata);
z           = amplitude.*exp(i*phase); % this is the pac measure

% get time windows
% ----------------
g.verbose = 'on';
g.causal  = 'off';
[ timesout1 indexout ] = gettimes(frame, g.tlimits, g.timesout, g.winsize, g.ntimesout, g.causal, g.verbose);

% scan time windows
% -----------------
if ~isempty(g.alpha)
    m_raw = zeros(1,length(indexout));
    pvals = zeros(1,length(indexout));
end;
fprintf('Computing PAC:\n');
for iWin = 1:length(indexout)
    x_phaseEpoch = x_freqphase(indexout(iWin)+[-g.winsize/2+1:g.winsize/2]);
    x_ampEpoch   = x_freqamp(  indexout(iWin)+[-g.winsize/2+1:g.winsize/2]);
    z_phaseEpoch = z_phasedata(indexout(iWin)+[-g.winsize/2+1:g.winsize/2]);
    z_ampEpoch   = z_ampdata(  indexout(iWin)+[-g.winsize/2+1:g.winsize/2]);
    z_epoch      = z(          indexout(iWin)+[-g.winsize/2+1:g.winsize/2]);
    
    numpoints=length(x_phaseEpoch);
    if rem(iWin,10) == 0,  verboseprintf(g.verbose, ' %d',iWin); end
    if rem(iWin,120) == 0, verboseprintf(g.verbose, '\n'); end
    
    % Choose method
    % -------------
    if strcmpi(g.method, 'modulation')
        
        % Modulation index
        m_raw(iWin) = abs(sum(z_epoch))/numpoints;
        
    elseif strcmpi(g.method, 'plv')
        
        if iWin == 145
            %dsfsd; 
        end;
        
        %amplitude_filt = sgolayfilt(amplitude, 3, 101);
        if ~isempty(g.filterphase)
             amplitude_filt = feval(g.filterphase, z_ampEpoch);
        else amplitude_filt = feval(g.filterfunc , z_ampEpoch, srate, g.freqphase(1), g.freqphase(end));
        end;
        z_amplitude_filt = hilbert(amplitude_filt);
        
        phase_amp_modulation = angle(z_amplitude_filt);
        m_raw(iWin) = abs(sum(exp(i*(x_phaseEpoch - phase_amp_modulation)))/numpoints);
        
    elseif strcmpi(g.method, 'corr')
        
        if iWin == inf %145
            figure; plot(abs(z_ampdata))
            hold on; plot(x_phasedata/10000000000, 'r')
            x = X(indexout(iWin)+[-g.winsize/2+1:g.winsize/2]);
            hold on; plot(x, 'g');
            dsfsd;
        end;
        [r_ESC pval_corr] = corrcoef(x_phaseEpoch, abs(z_ampEpoch));
        m_raw(iWin)   = r_ESC(1,2);
        pvals(iWin)   = pval_corr(1,2);
        
    elseif strcmpi(g.method, 'glm')
        
        [b dev stats] = glmfit(x_phaseEpoch', abs(z_ampEpoch)', 'normal');
        GLM_beta      = stats.beta(2,1);
        pvals(iWin)   = stats.p(2,1);
        m_raw(iWin)   = b(1);
        
    end;
    
    %% compute statistics (instantaneous)
    % -----------------------------------
    if ~isempty(g.alpha) && strcmpi(g.instantstat, 'on') && ~strcmpi(g.method, 'corr') && ~strcmpi(g.method, 'glm')
        
        % compute surrogate values
        numsurrogate=g.naccu;
        minskip=srate;
        maxskip=numpoints-srate; % max variation half a second
        if maxskip < 1
            error('Window size shorter than 1 second; too short for computing surrogate data');
        end;
        skip=ceil(numpoints.*rand(numsurrogate*4,1));
        skip(skip>maxskip)=[];
        skip(skip<minskip)=[];
        skip=skip(1:numsurrogate,1);
        surrogate_m=zeros(numsurrogate,1);
        for s=1:numsurrogate
            surrogate_amplitude=[amplitude(skip(s):end) amplitude(1:skip(s)-1)]; % consider circular shifts
            surrogate_m(s)=abs(mean(surrogate_amplitude.*exp(i*phase)));
            %disp(numsurrogate-s)
        end
        
        if strcmpi(g.statlim, 'surrogate')
            pvals(iWin) = stat_surrogate_pvals(surrogate_m, m_raw(iWin), 'upper');
            %fprintf('Raw PAC is %3.2f (p-value=%1.3f)\n', m_raw(iWin), pvals(iWin));
        
        else
            % Canolty method below
            
            %% fit gaussian to surrogate data, uses normfit.m from MATLAB Statistics toolbox
            [surrogate_mean,surrogate_std]=normfit(surrogate_m);
            
            %% normalize length using surrogate data (z-score)
            m_norm_length=(abs(m_raw(iWin))-surrogate_mean)/surrogate_std;
            pvals(iWin)  = normcdf(0, m_norm_length, 1);
            m_norm_phase=angle(m_raw(iWin));
            m_norm=m_norm_length*exp(i*m_norm_phase);
            
            % compare parametric and non-parametric methods (return similar
            % results)
            if iWin == length(indexout)
                figure;
                plot(-log10(pvals));
                hold on; plot(-log10(pvals2), 'r');
            end;
        end;
    end;
    
end;
fprintf('\n');

% Computes alpha
% --------------
if ~isempty(g.alpha) && strcmpi(g.instantstat, 'off')
    if isempty(g.baseline)
        g.baseline = [ timesout1(1) timesout1(end) ];
    end;
    baselineInd = find(timesout1 >= g.baseline(1) & timesout1 <= g.baseline(end));
    
    m_raw_base = abs(m_raw(baselineInd));
    
    if strcmpi(g.statlim, 'surrogate')
        for index = 1:length(m_raw)
            pvals(index) = stat_surrogate_pvals(m_raw_base, m_raw(index), 'upper');
        end;
    else
        [surrogate_mean,surrogate_std]=normfit(m_raw_base);
        m_norm_length=(abs(m_raw)-surrogate_mean)/surrogate_std;
        pvals = normcdf(0, m_norm_length, 1);
    end;
    if strcmpi(g.mcorrect, 'fdr')
        pvals = fdr(pvals);
    end;
end;

%% plot results
% -------------
if strcmpi(g.nofig, 'on')
    return
end;
if strcmpi(g.newfig, 'on')
    figure;
end;
if ~isempty(g.alpha)
    plotcurve(timesout1, m_raw, 'maskarray', pvals < g.alpha);
else plotcurve(timesout1, m_raw);
end;
xlabel('Time (ms)');
ylabel('PAC (0 to 1)');
title(g.title);

% plot vertical lines
% -------------------
if ~isempty(g.vert)
    hold on;
    yl = ylim;
    for index = 1:length(g.vert)
        plot([g.vert(index) g.vert(index)], yl, 'g');
    end;
end;

% -------------
% gettime function identical to timefreq function
% DO NOT MODIFY
% -------------
function [ timevals, timeindices ] = gettimes(frames, tlimits, timevar, winsize, ntimevar, causal, verbose);
timevect = linspace(tlimits(1), tlimits(2), frames);
srate    = 1000*(frames-1)/(tlimits(2)-tlimits(1));

if isempty(timevar) % no pre-defined time points
    if ntimevar(1) > 0
        % generate linearly space vector
        % ------------------------------
        if (ntimevar > frames-winsize)
            ntimevar = frames-winsize;
            if ntimevar < 0
                error('Not enough data points, reduce the window size or lowest frequency');
            end;
            verboseprintf(verbose, ['Value of ''timesout'' must be <= frame-winsize, ''timesout'' adjusted to ' int2str(ntimevar) '\n']);
        end
        npoints = ntimevar(1);
        wintime = 500*winsize/srate;
        if strcmpi(causal, 'on')
             timevals = linspace(tlimits(1)+2*wintime, tlimits(2), npoints);
        else timevals = linspace(tlimits(1)+wintime, tlimits(2)-wintime, npoints);
        end;
        verboseprintf(verbose, 'Generating %d time points (%1.1f to %1.1f ms)\n', npoints, min(timevals), max(timevals));
    else
        % subsample data
        % --------------
        nsub     = -ntimevar(1);
        if strcmpi(causal, 'on')
            timeindices = [ceil(winsize+nsub):nsub:length(timevect)];
        else timeindices = [ceil(winsize/2+nsub/2):nsub:length(timevect)-ceil(winsize/2)-1];
        end;
        timevals    = timevect( timeindices ); % the conversion at line 741 leaves timeindices unchanged
        verboseprintf(verbose, 'Subsampling by %d (%1.1f to %1.1f ms)\n', nsub, min(timevals), max(timevals));
    end;
else
    timevals = timevar;
    % check boundaries
    % ----------------
    wintime = 500*winsize/srate;
    if strcmpi(causal, 'on')
         tmpind  = find( (timevals >= tlimits(1)+2*wintime-0.0001) & (timevals <= tlimits(2)) );
    else tmpind  = find( (timevals >= tlimits(1)+wintime-0.0001) & (timevals <= tlimits(2)-wintime+0.0001) );
    end;
    % 0.0001 account for numerical innacuracies on opteron computers
    if isempty(tmpind)
        error('No time points. Reduce time window or minimum frequency.');
    end;
    if  length(timevals) ~= length(tmpind)
        verboseprintf(verbose, 'Warning: %d out of %d time values were removed (now %3.2f to %3.2f ms) so the lowest\n', ...
            length(timevals)-length(tmpind), length(timevals), timevals(tmpind(1)), timevals(tmpind(end)));
        verboseprintf(verbose, '         frequency could be computed with the requested accuracy\n');
    end;
    timevals = timevals(tmpind);
end;

% find closet points in data
% --------------------------
timeindices = round(eeg_lat2point(timevals, 1, srate, tlimits, 1E-3));
if length(timeindices) < length(unique(timeindices))
    timeindices = unique_bc(timeindices)
    verboseprintf(verbose, 'Warning: duplicate times, reduce the number of output times\n');
end;
if length(unique(timeindices(2:end)-timeindices(1:end-1))) > 1
    verboseprintf(verbose, 'Finding closest points for time variable\n');
    verboseprintf(verbose, 'Time values for time/freq decomposition is not perfectly uniformly distributed\n');
else
    verboseprintf(verbose, 'Distribution of data point for time/freq decomposition is perfectly uniform\n');
end;
timevals    = timevect(timeindices);

function verboseprintf(verbose, varargin)
if strcmpi(verbose, 'on')
    fprintf(varargin{:});
end;

