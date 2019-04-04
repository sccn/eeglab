% timefreq() - compute time/frequency decomposition of data trials. This 
%              function is a compute-only function called by
%              the more complete time/frequency functions newtimef()
%              and newcrossf() which also plot timefreq() results.
%
% Usage:
%     >> [tf, freqs, times]          = timefreq(data, srate);
%     >> [tf, freqs, times, itcvals] = timefreq(data, srate, ...
%                                        'key1', 'val1', 'key2', 'val2' ...)
% Inputs:
%         data    = [float array] 2-D data array of size (times,trials)
%         srate   = sampling rate
%
% Optional inputs:
%       'cycles'  = [real] indicates the number of cycles for the 
%                   time-frequency decomposition {default: 0}
%                   if 0, use FFTs and Hanning window tapering.  
%                   or [real positive scalar] Number of cycles in each Morlet
%                   wavelet, constant across frequencies.
%                   or [cycles cycles(2)] wavelet cycles increase with 
%                   frequency starting at cycles(1) and, 
%                   if cycles(2) > 1, increasing to cycles(2) at
%                   the upper frequency,
%                   or if cycles(2) = 0, same window size at all
%                   frequencies (similar to FFT if cycles(1) = 1)
%                   or if cycles(2) = 1, not increasing (same as giving
%                   only one value for 'cycles'). This corresponds to pure
%                   wavelet with the same number of cycles at each frequencies
%                   if 0 < cycles(2) < 1, linear variation in between pure 
%                   wavelets (1) and FFT (0). The exact number of cycles
%                   at the highest frequency is indicated on the command line.
%       'wavelet' = DEPRECATED, please use 'cycles'. This function does not 
%                   support multitaper. For multitaper, use timef().
%       'wletmethod' = ['dftfilt2'|'dftfilt3'] Wavelet method/program to use.
%                   {default: 'dftfilt3'}
%                   'dftfilt'  DEPRECATED. Method used in regular timef()
%                              program. Not available any more.
%                   'dftfilt2' Morlet-variant or Hanning DFT (calls dftfilt2()
%                              to generate wavelets).
%                   'dftfilt3' Morlet wavelet or Hanning DFT (exact Tallon 
%                              Baudry). Calls dftfilt3().
%       'ffttaper' = ['none'|'hanning'|'hamming'|'blackmanharris'] FFT tapering
%                   function. Default is 'hanning'. Note that 'hamming' and 
%                   'blackmanharris' require the signal processing toolbox.
% Optional ITC type:
%        'type'   = ['coher'|'phasecoher'] Compute either linear coherence
%                   ('coher') or phase coherence ('phasecoher') also known
%                   as phase coupling factor' {default: 'phasecoher'}.
%        'subitc' = ['on'|'off'] subtract stimulus locked Inter-Trial Coherence
%                   (ITC) from x and y. This computes the  'intrinsic' coherence
%                   x and y not arising from common synchronization to
%                   experimental events. See notes. {default: 'off'}
%
% Optional detrending:
%       'detrend' = ['on'|'off'], Linearly detrend each data epoch  {'off'}
%
% Optional FFT/DFT parameters:
%      'tlimits'  = [min max] time limits in ms.
%      'winsize'  = If cycles==0 (FFT, see 'wavelet' input): data subwindow
%                   length (fastest, 2^n<frames);
%                   if cycles >0: *longest* window length to use. This
%                   determines the lowest output frequency  {~frames/8}
%      'ntimesout' = Number of output times (int<frames-winsize). Enter a
%                   negative value [-S] to subsample original time by S.
%      'timesout' = Enter an array to obtain spectral decomposition at
%                   specific time values (note: algorithm find closest time
%                   point in data and this might result in an unevenly spaced
%                   time array). Overwrite 'ntimesout'. {def: automatic}
%      'freqs'    = [min max] frequency limits. Default [minfreq srate/2],
%                   minfreq being determined by the number of data points,
%                   cycles and sampling frequency. Enter a single value
%                   to compute spectral decompisition at a single frequency
%                   (note: for FFT the closest frequency will be estimated).
%                   For wavelet, reducing the max frequency reduce
%                   the computation load.
%      'padratio' = FFTlength/winsize (2^k)                     {def: 2}
%                   Multiplies the number of output frequencies by
%                   dividing their spacing. When cycles==0, frequency
%                   spacing is (low_frequency/padratio).
%      'nfreqs'   = number of output frequencies. For FFT, closest computed
%                   frequency will be returned. Overwrite 'padratio' effects
%                   for wavelets. Default: use 'padratio'.
%     'freqscale' = ['log'|'linear'] frequency scale. Default is 'linear'.
%                   Note that for obtaining 'log' spaced freqs using FFT,
%                   closest correspondant frequencies in the 'linear' space
%                   are returned.
%     'wletmethod'= ['dftfilt2'|'dftfilt3'] Wavelet method/program to use.
%                   Default is 'dftfilt3'
%                   'dftfilt3' Morlet wavelet or Hanning DFT
%                   'dftfilt2' Morlet-variant or Hanning DFT.
%                   Note that there are differences betweeen the Hanning
%                   DFTs in the two programs.
%       'causal'  = ['on'|'off'] apply FFT or time-frequency in a causal
%                   way where only data before any given latency can 
%                   influence the spectral decomposition. (default: 'off'}
%
% Optional time warping:
%   'timestretch' = {[Refmarks], [Refframes]}
%                   Stretch amplitude and phase time-course before smoothing.
%                   Refmarks is a (trials,eventframes) matrix in which rows are
%                   event marks to time-lock to for a given trial. Each trial
%                   will be stretched so that its marked events occur at frames
%                   Refframes. If Refframes is [], the median frame, across trials,
%                   of the Refmarks latencies will be used. Both Refmarks and
%                   Refframes are given in frames in this version - will be
%                   changed to ms in future.
% Outputs:
%         tf      = complex time frequency array for all trials (freqs, 
%                   times, trials)
%         freqs   = vector of computed frequencies (Hz)
%         times   = vector of computed time points (ms)
%         itcvals = time frequency "average" for all trials (freqs, times).
%                   In the coherence case, it is the real mean of the time
%                   frequency decomposition, but in the phase coherence case
%                   (see 'type' input'), this is the mean of the normalized
%                   spectral estimate.
%
% Authors: Arnaud Delorme, Jean Hausser & Scott Makeig
%          CNL/Salk Institute 1998-2001; SCCN/INC/UCSD, La Jolla, 2002-
%          Fix FFT frequency innacuracy, bug 874 by WuQiang
%
% See also: timef(), newtimef(), crossf(), newcrossf()

% Note: it is not advised to use a FFT decomposition in a log scale. Output
%       value are accurate but plotting might not be because of the non-uniform
%       frequency output values in log-space. If you have to do it, use a
%       padratio as large as possible, or interpolate time-freq image at
%       exact log scale values before plotting.

% Copyright (C) 8/1/98  Arnaud Delorme, Sigurd Enghoff & Scott Makeig, SCCN/INC/UCSD
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

function [tmpall, freqs, timesout, itcvals] = timefreq(data, srate, varargin)

if nargin < 2
    help timefreq;
    return;
end

[chan frame trials]= size(data);
if trials == 1 && chan ~= 1
    trials = frame;
    frame  = chan;
    chan   = 1;
end
g = finputcheck(varargin, ...
    { 'ntimesout'     'integer'  []                     []; ...
    'timesout'      'real'     []                       []; ...
    'winsize'       'integer'  [0 Inf]                  []; ...
    'tlimits'       'real'     []                       []; ...
    'detrend'       'string'   {'on','off'}             'off'; ...
    'causal'        'string'   {'on','off'}             'off'; ...
    'verbose'       'string'   {'on','off'}             'on'; ...
    'freqs'         'real'     [0 Inf]                  []; ...
    'nfreqs'        'integer'  [0 Inf]                  []; ...
    'freqscale'     'string'   { 'linear','log','' }    'linear'; ...
    'ffttaper'      'string'   { 'hanning','hamming','blackmanharris','none' }  'hanning';
    'wavelet'       'real'     [0 Inf]                  0; ...
    'cycles'        {'real','integer'}    [0 Inf]       0; ...
    'padratio'      'integer'  [1 Inf]                  2; ...
    'itctype'       'string'   {'phasecoher','phasecoher2','coher'}  'phasecoher'; ...
    'subitc'        'string'   {'on','off'}             'off'; ...
    'timestretch'   'cell'     []                       {}; ...
    'wletmethod'    'string'   {'dftfilt2','dftfilt3'}    'dftfilt3'; ...
    });
if ischar(g), error(g); end
if isempty(g.freqscale), g.freqscale = 'linear'; end
if isempty(g.winsize),   g.winsize   = max(pow2(nextpow2(frame)-3),4); end
if isempty(g.ntimesout), g.ntimesout = 200; end
if isempty(g.freqs),     g.freqs     = [0 srate/2]; end
if isempty(g.tlimits),   g.tlimits   = [0 frame/srate*1000]; end

% checkin parameters
% ------------------

% Use 'wavelet' if 'cycles' undefined for backwards compatibility
if g.cycles == 0
    g.cycles = g.wavelet;
end

if (g.winsize > frame)
    error('Value of winsize must be less than frame length.');
end
if (pow2(nextpow2(g.padratio)) ~= g.padratio)
    error('Value of padratio must be an integer power of two [1,2,4,8,16,...]');
end

% finding frequency limits
% ------------------------
if g.cycles(1) ~= 0 && g.freqs(1) == 0, g.freqs(1) = srate*g.cycles(1)/g.winsize; end

% finding frequencies
% -------------------
if length(g.freqs) == 2

    % min and max
    % -----------
    if g.freqs(1) == 0 && g.cycles(1) ~= 0
        g.freqs(1) = srate*g.cycles(1)/g.winsize;
    end

    % default number of freqs using padratio
    % --------------------------------------
    if isempty(g.nfreqs)
        g.nfreqs = g.winsize/2*g.padratio+1;
        % adjust nfreqs depending on frequency range
        tmpfreqs = linspace(0, srate/2, g.nfreqs);
        tmpfreqs = tmpfreqs(2:end);  % remove DC (match the output of PSD)

        % adjust limits for FFT (only linear scale)
        if g.cycles(1) == 0 && ~strcmpi(g.freqscale, 'log')
            if ~any(tmpfreqs == g.freqs(1))
                [tmp minind] = min(abs(tmpfreqs-g.freqs(1)));
                g.freqs(1)   = tmpfreqs(minind);
                verboseprintf(g.verbose, 'Adjust min freq. to %3.2f Hz to match FFT output frequencies\n', g.freqs(1));
            end
            if ~any(tmpfreqs == g.freqs(2))
                [tmp minind] = min(abs(tmpfreqs-g.freqs(2)));
                g.freqs(2)   = tmpfreqs(minind);
                verboseprintf(g.verbose, 'Adjust max freq. to %3.2f Hz to match FFT output frequencies\n', g.freqs(2));
            end
        end

        % find number of frequencies
        % --------------------------
        g.nfreqs = length(tmpfreqs( intersect( find(tmpfreqs >= g.freqs(1)), find(tmpfreqs <= g.freqs(2)))));
        if g.freqs(1)==g.freqs(2), g.nfreqs = 1; end
    end

    % find closest freqs for FFT
    % --------------------------
    if strcmpi(g.freqscale, 'log')
        g.freqs = linspace(log(g.freqs(1)), log(g.freqs(end)), g.nfreqs);
        g.freqs = exp(g.freqs);
    else
        g.freqs = linspace(g.freqs(1), g.freqs(2), g.nfreqs); % this should be OK for FFT
        % because of the limit adjustment
    end
end
g.nfreqs = length(g.freqs);

% function for time freq initialisation
% -------------------------------------
if (g.cycles(1) == 0) %%%%%%%%%%%%%% constant window-length FFTs %%%%%%%%%%%%%%%%
    freqs = linspace(0, srate/2, g.winsize*g.padratio/2+1);
    freqs = freqs(2:end); % remove DC (match the output of PSD)
    %srate/g.winsize*[1:2/g.padratio:g.winsize]/2
    verboseprintf(g.verbose, 'Using %s FFT tapering\n', g.ffttaper);
    switch g.ffttaper
        case 'hanning',        g.win   = hanning(g.winsize);
        case 'hamming',        g.win   = hamming(g.winsize);
        case 'blackmanharris', g.win   = blackmanharris(g.winsize);
        case 'none',           g.win   = ones(g.winsize,1);
    end
else % %%%%%%%%%%%%%%%%%% Constant-Q (wavelet) DFTs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %freqs = srate*g.cycles/g.winsize*[2:2/g.padratio:g.winsize]/2;
    %g.win = dftfilt(g.winsize,g.freqs(2)/srate,g.cycles,g.padratio,g.cyclesfact);

    freqs = g.freqs;
    if length(g.cycles) == 2
        if g.cycles(2) < 1
            g.cycles = [ g.cycles(1) g.cycles(1)*g.freqs(end)/g.freqs(1)*(1-g.cycles(2))];
        end
        verboseprintf(g.verbose, 'Using %g cycles at lowest frequency to %g at highest.\n', g.cycles(1), g.cycles(2));
    elseif length(g.cycles) == 1
        verboseprintf(g.verbose, 'Using %d cycles at all frequencies.\n',g.cycles);
    else
        verboseprintf(g.verbose, 'Using user-defined cycle for each frequency\n');
    end
    if strcmp(g.wletmethod, 'dftfilt2')
        g.win    = dftfilt2(g.freqs,g.cycles,srate, g.freqscale); % uses Morlet taper by default
    elseif strcmp(g.wletmethod, 'dftfilt3')     % Default
        g.win    = dftfilt3(g.freqs,g.cycles,srate, 'cycleinc', g.freqscale); % uses Morlet taper by default
    else return
    end
    g.winsize = 0;
    for index = 1:length(g.win)
        g.winsize = max(g.winsize,length(g.win{index}));
    end
end

% compute time vector
% -------------------
[ g.timesout g.indexout ] = gettimes(frame, g.tlimits, g.timesout, g.winsize, g.ntimesout, g.causal, g.verbose);

% -------------------------------
% compute time freq decomposition
% -------------------------------
verboseprintf(g.verbose, 'The window size used is %d samples (%g ms) wide.\n',g.winsize, 1000/srate*g.winsize);
if strcmpi(g.freqscale, 'log') % fastif was having strange "function not available" messages
    scaletoprint = 'log';
else scaletoprint = 'linear';
end
verboseprintf(g.verbose, 'Estimating %d %s-spaced frequencies from %2.1f Hz to %3.1f Hz.\n', length(g.freqs), ...
    scaletoprint, g.freqs(1), g.freqs(end));
%verboseprintf(g.verbose, 'Estimating %d %s-spaced frequencies from %2.1f Hz to %3.1f Hz.\n', length(g.freqs), ...
%    fastif(strcmpi(g.freqscale, 'log'), 'log', 'linear'), g.freqs(1), g.freqs(end));

if g.cycles(1) == 0
    if 1
        % build large matrix to compute FFT
        % ---------------------------------
        indices = repmat([-g.winsize/2+1:g.winsize/2]', [1 length(g.indexout) trials]);
        indices = indices + repmat(g.indexout, [size(indices,1) 1 trials]);
        indices = indices + repmat(reshape(([1:trials]-1)*frame,1,1,trials), [size(indices,1) length(g.indexout) 1]);
        if chan > 1
            tmpall = repmat(nan,[chan length(freqs) length(g.timesout) trials]);
            tmpX = reshape(data(:,indices), [ size(data,1) size(indices)]);
            tmpX = bsxfun(@minus, tmpX, mean( tmpX, 2)); % avoids repmat - faster than tmpX = tmpX - repmat(mean(tmpX), [size(tmpX,1) 1 1]);
            tmpX = bsxfun(@times, tmpX, g.win');
            tmpX = fft(tmpX,g.padratio*g.winsize,2);
            tmpall = squeeze(tmpX(:,2:g.padratio*g.winsize/2+1,:,:));
        else
            tmpall = repmat(nan,[length(freqs) length(g.timesout) trials]);
            tmpX    = data(indices);
            tmpX = bsxfun(@minus, tmpX, mean( tmpX, 1)); % avoids repmat - faster than tmpX = tmpX - repmat(mean(tmpX), [size(tmpX,1) 1 1]);
            tmpX = bsxfun(@times, tmpX, g.win);
            %tmpX = fft(tmpX,2^ceil(log2(g.padratio*g.winsize)));
            %tmpall = tmpX(2:g.padratio*g.winsize/2+1,:,:);
            tmpX = fft(tmpX,g.padratio*g.winsize);
            tmpall = tmpX(2:g.padratio*g.winsize/2+1,:,:);
        end
    else % old iterative computation
        tmpall = repmat(nan,[length(freqs) length(g.timesout) trials]);
        verboseprintf(g.verbose, 'Processing trial (of %d):',trials);
        for trial = 1:trials
            if rem(trial,10) == 0,  verboseprintf(g.verbose, ' %d',trial); end
            if rem(trial,120) == 0, verboseprintf(g.verbose, '\n'); end
            for index = 1:length(g.indexout)
                if strcmpi(g.causal, 'off')
                    tmpX = data([-g.winsize/2+1:g.winsize/2]+g.indexout(index)+(trial-1)*frame); % 1 point imprecision
                else
                    tmpX = data([-g.winsize+1:0]+g.indexout(index)+(trial-1)*frame); % 1 point imprecision
                end

                tmpX = tmpX - mean(tmpX);
                if strcmpi(g.detrend, 'on'),
                    tmpX = detrend(tmpX);
                end

                tmpX = g.win .* tmpX(:);
                tmpX = fft(tmpX,g.padratio*g.winsize);
                tmpX = tmpX(2:g.padratio*g.winsize/2+1);
                tmpall(:,index, trial) = tmpX(:);
            end
        end
    end
else % wavelet
    if chan > 1
        % wavelets are processed in groups of the same size
        % to speed up computation. Wavelet of groups of different size
        % can be processed together but at a cost of a lot of RAM and 
        % a lot of extra computation -> not efficient
        tmpall = repmat(nan,[chan length(freqs) length(g.timesout) trials]);
        wt = [ 1 find(diff(cellfun(@length,g.win)))+1 length(g.win)+1];
        verboseprintf(g.verbose, 'Computing of %d:', length(wt));
        for ind = 1:length(wt)-1
            verboseprintf(g.verbose, '.');
            wavarray = reshape([ g.win{wt(ind):wt(ind+1)-1} ], [ length(g.win{wt(ind)}) wt(ind+1)-wt(ind) ]);
            sizewav = size(wavarray,1)-1;
            indices = repmat([-sizewav/2:sizewav/2]', [1 size(wavarray,2) length(g.indexout) trials]);
            indices = indices + repmat(reshape(g.indexout, 1,1,length(g.indexout)), [size(indices,1) size(indices,2) 1 trials]);
            indices = indices + repmat(reshape(([1:trials]-1)*frame,1,1,1,trials),  [size(indices,1) size(indices,2) size(indices,3) 1]);
            szfreqdata = [ size(data,1) size(indices) ];
            tmpX    = reshape(data(:,indices), szfreqdata);
            tmpX    = bsxfun(@minus, tmpX, mean( tmpX, 2)); % avoids repmat - faster than tmpX = tmpX - repmat(mean(tmpX), [size(tmpX,1) 1 1]);
            wavarray = reshape(wavarray, [1 size(wavarray,1) size(wavarray,2)]);
            tmpall(:,wt(ind):wt(ind+1)-1,:,:,:) = reshape(sum(bsxfun(@times, tmpX, wavarray),2), [szfreqdata(1) szfreqdata(3:end)]);
        end
        verboseprintf(g.verbose, '\n');
        %tmpall = squeeze(tmpall(1,:,:,:));
    elseif 0
        tmpall = repmat(nan,[length(freqs) length(g.timesout) trials]);
        % wavelets are processed in groups of the same size
        % to speed up computation. Wavelet of groups of different size
        % can be processed together but at a cost of a lot of RAM and 
        % a lot of extra computation -> not faster than the regular
        % iterative method
        wt = [ 1 find(diff(cellfun(@length,g.win)))+1 length(g.win)+1];
        for ind = 1:length(wt)-1
            wavarray = reshape([ g.win{wt(ind):wt(ind+1)-1} ], [ length(g.win{wt(ind)}) wt(ind+1)-wt(ind) ]);
            sizewav = size(wavarray,1)-1;
            indices = repmat([-sizewav/2:sizewav/2]', [1 size(wavarray,2) length(g.indexout) trials]);
            indices = indices + repmat(reshape(g.indexout, 1,1,length(g.indexout)), [size(indices,1) size(indices,2) 1 trials]);
            indices = indices + repmat(reshape(([1:trials]-1)*frame,1,1,1,trials), [size(indices,1) size(indices,2) size(indices,3) 1]);
            tmpX    = data(indices);
            tmpX    = bsxfun(@minus, tmpX, mean( tmpX, 1)); % avoids repmat - faster than tmpX = tmpX - repmat(mean(tmpX), [size(tmpX,1) 1 1]);
            tmpall(wt(ind):wt(ind+1)-1,:,:)  = squeeze(sum(bsxfun(@times, tmpX, wavarray),1));
        end
    elseif 0
        % wavelets are processed one by one but all windows simultaneously
        % -> not faster than the regular iterative method
        tmpall  = repmat(nan,[length(freqs) length(g.timesout) trials]);
        sizewav = length(g.win{1})-1; % max window size
        mainc   = sizewav/2;
        indices = repmat([-sizewav/2:sizewav/2]', [1 length(g.indexout) trials]);
        indices = indices + repmat(g.indexout, [size(indices,1) 1 trials]);
        indices = indices + repmat(reshape(([1:trials]-1)*frame,1,1,trials), [size(indices,1) length(g.indexout) 1]);
        
        for freqind = 1:length(g.win)
            winc = (length(g.win{freqind})-1)/2;
            wins = length(g.win{freqind})-1;
            wini = [-wins/2:wins/2]+winc+mainc-winc+1;
            tmpX = data(indices(wini,:,:));
            tmpX = bsxfun(@minus, tmpX, mean( tmpX, 1)); % avoids repmat - faster than tmpX = tmpX - repmat(mean(tmpX), [size(tmpX,1) 1 1]);
            tmpX = sum(bsxfun(@times, tmpX, g.win{freqind}'),1);
            tmpall(freqind,:,:) = tmpX;
        end
    else
        % prepare wavelet filters
        % -----------------------
        for index = 1:length(g.win)
            g.win{index} = transpose(repmat(g.win{index}, [trials 1]));
        end

        % apply filters
        % -------------
        verboseprintf(g.verbose, 'Processing time point (of %d):',length(g.timesout));
        tmpall = zeros(length(g.win), length(g.indexout), size(data,2));
        for index = 1:length(g.indexout)
            if rem(index,10) == 0,  verboseprintf(g.verbose, ' %d',index); end
            if rem(index,120) == 0, verboseprintf(g.verbose, '\n'); end
            for freqind = 1:length(g.win)
                wav = g.win{freqind}; 
                sizewav = size(wav,1)-1;
                %g.indexout(index), size(wav,1), g.freqs(freqind)
                if strcmpi(g.causal, 'off')
                    tmpX = data([-sizewav/2:sizewav/2]+g.indexout(index),:);
                else
                    tmpX = data([-sizewav:0]+g.indexout(index),:);
                end

                tmpX = tmpX - ones(size(tmpX,1),1)*mean(tmpX);
                if strcmpi(g.detrend, 'on'),
                    for trial = 1:trials
                        tmpX(:,trial) = detrend(tmpX(:,trial));
                    end
                end

                tmpX = sum(wav .* tmpX);
                tmpall( freqind, index, :) = tmpX;
            end
        end
    end
end
verboseprintf(g.verbose, '\n');

% time-warp code begins -Jean
% ---------------------------
if ~isempty(g.timestretch) && length(g.timestretch{1}) > 0

    timemarks = g.timestretch{1}';
    if isempty(g.timestretch{2}) || length(g.timestretch{2}) == 0
        timerefs = median(g.timestretch{1}',2);
    else
        timerefs = g.timestretch{2};
    end
    trials = size(tmpall,3);

    % convert timerefs to subsampled ERSP space
    % -----------------------------------------

    [dummy refsPos] = min(transpose(abs( ...
        repmat(timerefs, [1 length(g.indexout)]) - repmat(g.indexout, [length(timerefs) 1]))));
    refsPos(end+1) = 1;
    refsPos(end+1) = length(g.indexout);
    refsPos = sort(refsPos);

    for t=1:trials

        % convert timemarks to subsampled ERSP space
        % ------------------------------------------

        %[dummy pos]=min(abs(repmat(timemarks(2:7,1), [1 length(g.indexout)])-repmat(g.indexout,[6 1])));

        outOfTimeRangeTimeWarpMarkers = find(timemarks(:,t) < min(g.indexout) | timemarks(:,t) > max(g.indexout));
        
%         if ~isempty(outOfTimeRangeTimeWarpMarkers)
%             verboseprintf(g.verbose, 'Timefreq warning: time-warp latencies in epoch %d are out of time range defined for calculation of ERSP.\n', t);
%         end
        
        [dummy marksPos] = min(transpose( ...
            abs( ...
            repmat(timemarks(:,t), [1 length(g.indexout)]) ...
            - repmat(g.indexout, [size(timemarks,1) 1]) ...
            ) ...
            ));
         

        marksPos(end+1) = 1;
        marksPos(end+1) = length(g.indexout);
        marksPos = sort(marksPos);

        %now warp tmpall
        mytmpall = tmpall(:,:,t);
        r = sqrt(mytmpall.*conj(mytmpall));
        theta = angle(mytmpall);

        % So mytmpall is almost equal to r.*exp(i*theta)
        % whos marksPos refsPos

        M = timewarp(marksPos, refsPos);
        
        TSr = transpose(M*r');
        TStheta = zeros(size(theta,1), size(theta,2));

        for freqInd=1:size(TStheta,1)
            TStheta(freqInd, :) = angtimewarp(marksPos, refsPos, theta(freqInd, :));            
        end
        TStmpall = TSr.*exp(i*TStheta);

        % $$$     keyboard;

        tmpall(:,:,t) =  TStmpall;
    end
end

%time-warp ends
zerovals = tmpall == 0;
if any(reshape(zerovals, 1, prod(size(zerovals))))
    tmpall(zerovals) = Inf;
    minval = min(tmpall(:)); % remove bug
    tmpall(zerovals) = minval;
end

% compute and subtract ITC
% ------------------------
if nargout > 3 || strcmpi(g.subitc, 'on')
    itcvals = newtimefitc(tmpall, g.itctype);
end
if strcmpi(g.subitc, 'on')
    %a = gcf; figure; imagesc(abs(itcvals)); cbar; figure(a);
    if ndims(tmpall) <= 3
         tmpall = (tmpall - abs(tmpall) .* repmat(itcvals, [1 1 trials])) ./ abs(tmpall);
    else tmpall = (tmpall - abs(tmpall) .* repmat(itcvals, [1 1 1 trials])) ./ abs(tmpall);
    end
end

% find closest output frequencies
% -------------------------------
if length(g.freqs) ~= length(freqs) || any(g.freqs ~= freqs)
    allindices = zeros(1,length(g.freqs));
    for index = 1:length(g.freqs)
        [dum ind] = min(abs(freqs-g.freqs(index)));
        allindices(index) = ind;
    end
    verboseprintf(g.verbose, 'finding closest frequencies: %d freqs removed\n', length(freqs)-length(allindices));
    freqs = freqs(allindices);
    if ndims(tmpall) <= 3
         tmpall = tmpall(allindices,:,:);
    else tmpall = tmpall(:,allindices,:,:);
    end
    if nargout > 3 || strcmpi(g.subitc, 'on')
        if ndims(tmpall) <= 3
             itcvals = itcvals(allindices,:,:);
        else itcvals = itcvals(:,allindices,:,:);
        end
    end
end

timesout = g.timesout;

%figure; imagesc(abs(sum(itcvals,3))); cbar;
return;

function w = hanning(n)
if ~rem(n,2)
    w = .5*(1 - cos(2*pi*(1:n/2)'/(n+1)));
    w = [w; w(end:-1:1)];
else
    w = .5*(1 - cos(2*pi*(1:(n+1)/2)'/(n+1)));
    w = [w; w(end-1:-1:1)];
end

% get time points
% ---------------
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
            end
            verboseprintf(verbose, ['Value of ''timesout'' must be <= frame-winsize, ''timesout'' adjusted to ' int2str(ntimevar) '\n']);
        end
        npoints = ntimevar(1);
        wintime = 500*winsize/srate;
        if strcmpi(causal, 'on')
             timevals = linspace(tlimits(1)+2*wintime, tlimits(2), npoints);
        else timevals = linspace(tlimits(1)+wintime, tlimits(2)-wintime, npoints);
        end
        verboseprintf(verbose, 'Generating %d time points (%1.1f to %1.1f ms)\n', npoints, min(timevals), max(timevals));
    else
        % subsample data
        % --------------
        nsub     = -ntimevar(1);
        if strcmpi(causal, 'on')
             timeindices = [ceil(winsize+nsub):nsub:length(timevect)];
        else timeindices = [ceil(winsize/2+nsub/2):nsub:length(timevect)-ceil(winsize/2)-1];
        end
        timevals    = timevect( timeindices ); % the conversion at line 741 leaves timeindices unchanged
        verboseprintf(verbose, 'Subsampling by %d (%1.1f to %1.1f ms)\n', nsub, min(timevals), max(timevals));
    end
else
    timevals = timevar;
    % check boundaries
    % ----------------
    wintime = 500*winsize/srate;
    if strcmpi(causal, 'on')
         tmpind  = find( (timevals >= tlimits(1)+2*wintime-0.0001) & (timevals <= tlimits(2)) ); 
    else tmpind  = find( (timevals >= tlimits(1)+wintime-0.0001) & (timevals <= tlimits(2)-wintime+0.0001) ); 
    end
    % 0.0001 account for numerical innacuracies on opteron computers
    if isempty(tmpind)
        error('No time points. Reduce time window or minimum frequency.');
    end
    if  length(timevals) ~= length(tmpind)
        verboseprintf(verbose, 'Warning: %d out of %d time values were removed (now %3.2f to %3.2f ms) so the lowest\n', ...
            length(timevals)-length(tmpind), length(timevals), timevals(tmpind(1)), timevals(tmpind(end)));
        verboseprintf(verbose, '         frequency could be computed with the requested accuracy\n');
    end
    timevals = timevals(tmpind);
end

% find closet points in data
% --------------------------
timeindices = round(eeg_lat2point(timevals, 1, srate, tlimits, 1E-3));
if length(timeindices) < length(unique(timeindices))
    timeindices = unique_bc(timeindices)
    verboseprintf(verbose, 'Warning: duplicate times, reduce the number of output times\n');
end
if length(unique(timeindices(2:end)-timeindices(1:end-1))) > 1
    verboseprintf(verbose, 'Finding closest points for time variable\n');
    verboseprintf(verbose, 'Time values for time/freq decomposition is not perfectly uniformly distributed\n');
else
    verboseprintf(verbose, 'Distribution of data point for time/freq decomposition is perfectly uniform\n');
end
timevals    = timevect(timeindices);

% DEPRECATED, FOR C INTERFACE
function nofunction()
% C PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = [ 'tmpcrossf' num2str(round(rand(1)*1000)) ];
f = fopen([ filename '.in'], 'w');
fwrite(f, tmpsaveall, 'int32');
fwrite(f, g.detret, 'int32');
fwrite(f, g.srate, 'int32');
fwrite(f, g.maxfreq, 'int32');
fwrite(f, g.padratio, 'int32');
fwrite(f, g.cycles, 'int32');
fwrite(f, g.winsize, 'int32');
fwrite(f, g.timesout, 'int32');
fwrite(f, g.subitc, 'int32');
fwrite(f, g.type, 'int32');
fwrite(f, trials, 'int32');
fwrite(f, g.naccu, 'int32');
fwrite(f, length(X), 'int32');
fwrite(f, X, 'double');
fwrite(f, Y, 'double');
fclose(f);

command = [ '!cppcrosff ' filename '.in ' filename '.out' ];
eval(command);

f = fopen([ filename '.out'], 'r');
size1 = fread(f, 'int32', 1);
size2 = fread(f, 'int32', 1);
Rreal = fread(f, 'double', [size1 size2]);
Rimg  = fread(f, 'double', [size1 size2]);
Coher.R = Rreal + j*Rimg;
Boot.Coherboot.R = [];
Boot.Rsignif = [];

function verboseprintf(verbose, varargin)
if strcmpi(verbose, 'on')
    fprintf(varargin{:});
end

