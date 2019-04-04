% pac() - compute phase-amplitude coupling (power of first input
%         correlation with phase of second). There is no graphical output
%         to this function.
%
% Usage:
%   >> pac(x,y,srate);
%   >> [coh,timesout,freqsout1,freqsout2,cohboot] ...
%                     = pac(x,y,srate,'key1', 'val1', 'key2', val2' ...);
% Inputs:
%    x       = [float array] 2-D data array of size (times,trials) or
%              3-D (1,times,trials)
%    y       = [float array] 2-D or 3-d data array
%    srate   = data sampling rate (Hz)
%
%    Most important optional inputs
%       'method'    = ['mod'|'corrsin'|'corrcos'|'latphase'] modulation
%                     method or correlation of amplitude with sine or cosine of 
%                     angle (see ref). 'laphase' compute the phase
%                     histogram at a specific time and requires the
%                     'powerlat' option to be set.
%       'freqs'     = [min max] frequency limits. Default [minfreq 50], 
%                     minfreq being determined by the number of data points, 
%                     cycles and sampling frequency. Use 0 for minimum frequency
%                     to compute default minfreq. You may also enter an 
%                     array of frequencies for the spectral decomposition
%                     (for FFT, closest computed frequency will be returned; use
%                     'padratio' to change FFT freq. resolution).
%       'freqs2'    = [float array] array of frequencies for the second
%                     argument. 'freqs' is used for the first argument. 
%                     By default it is the same as 'freqs'.
%       'wavelet'   = 0  -> Use FFTs (with constant window length) { Default } 
%                   = >0 -> Number of cycles in each analysis wavelet 
%                   = [cycles expfactor] -> if 0 < expfactor < 1,  the number 
%                     of wavelet cycles expands with frequency from cycles
%                     If expfactor = 1, no expansion; if = 0, constant
%                     window length (as in FFT)            {default wavelet: 0}
%                   = [cycles array] -> cycle for each frequency. Size of array
%                      must be the same as the number of frequencies 
%                     {default cycles: 0}
%       'wavelet2'  = same as 'wavelet' for the second argument. Default is
%                     same as cycles. Note that if the lowest frequency for X
%                     and Y are different and cycle is [cycles expfactor], it
%                     may result in discrepencies in the number of cycles at
%                     the same frequencies for X and Y.
%       'ntimesout' = Number of output times (int<frames-winframes). Enter a 
%                     negative value [-S] to subsample original time by S.
%       'timesout'  = Enter an array to obtain spectral decomposition at 
%                     specific time values (note: algorithm find closest time 
%                     point in data and this might result in an unevenly spaced
%                     time array). Overwrite 'ntimesout'. {def: automatic}
%       'powerlat'  = [float] latency in ms at which to compute phase
%                     histogram
%       'tlimits'   = [min max] time limits in ms.
%
%    Optional Detrending:
%       'detrend'   = ['on'|'off'], Linearly detrend each data epoch   {'off'}
%       'rmerp'     = ['on'|'off'], Remove epoch mean from data epochs {'off'}
%
%    Optional FFT/DFT Parameters:
%       'winsize'   = If cycles==0: data subwindow length (fastest, 2^n<frames);
%                     If cycles >0: *longest* window length to use. This
%                     determines the lowest output frequency. Note that this
%                     parameter is overwritten if the minimum frequency has been set
%                     manually and requires a longer time window {~frames/8}
%       'padratio'  = FFT-length/winframes (2^k)                    {2}
%                     Multiplies the number of output frequencies by dividing
%                     their spacing (standard FFT padding). When cycles~=0, 
%                     frequency spacing is divided by padratio.
%       'nfreqs'    = number of output frequencies. For FFT, closest computed
%                     frequency will be returned. Overwrite 'padratio' effects
%                     for wavelets. Default: use 'padratio'.
%       'freqscale' = ['log'|'linear'] frequency scale. Default is 'linear'.
%                     Note that for obtaining 'log' spaced freqs using FFT, 
%                     closest correspondant frequencies in the 'linear' space 
%                     are returned.
%       'subitc'    = ['on'|'off'] subtract stimulus locked Inter-Trial Coherence 
%                    (ITC) from x and y. This computes the  'intrinsic' coherence
%                     x and y not arising from common synchronization to 
%                     experimental events. See notes. {default: 'off'}
%       'itctype'   = ['coher'|'phasecoher'] For use with 'subitc', see timef()
%                     for more details {default: 'phasecoher'}.
%       'subwin'    = [min max] sub time window in ms (this windowing is
%                     performed after the spectral decomposition).
%
% Outputs: 
%        pac         = Matrix (nfreqs1,nfreqs2,timesout) of coherence (complex).
%                      Use 20*log(abs(crossfcoh)) to vizualize log spectral diffs. 
%        timesout    = Vector of output times (window centers) (ms).
%        freqsout1   = Vector of frequency bin centers for first argument (Hz).
%        freqsout2   = Vector of frequency bin centers for second argument (Hz).
%        alltfX      = single trial spectral decomposition of X
%        alltfY      = single trial spectral decomposition of Y
%
% Author: Arnaud Delorme, SCCN/INC, UCSD 2005-
%
% Ref: Testing for Nested Oscilations (2008) J Neuro Methods 174(1):50-61
%
% See also: timefreq(), crossf()

% Copyright (C) 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [crossfcoh, timesout1, freqs1, freqs2, crossfcohall, alltfX, alltfY] = pac(X, Y, srate, varargin);
    
if nargin < 1
    help pac; 
    return; 
end

% deal with 3-D inputs
% --------------------
if ndims(X) == 3, X = reshape(X, size(X,2), size(X,3)); end
if ndims(Y) == 3, Y = reshape(Y, size(Y,2), size(Y,3)); end
frame = size(X,2);

g = finputcheck(varargin, ...
                { 'alpha'         'real'     [0 0.2]                  [];
                  'baseboot'      'float'    []                       0;
                  'boottype'      'string'   {'times','trials','timestrials'}  'timestrials';
                  'detrend'       'string'   {'on','off'}              'off';
                  'freqs'         'real'     [0 Inf]                  [0 srate/2];
                  'freqs2'        'real'     [0 Inf]                  [];
                  'freqscale'     'string'   { 'linear','log' }       'linear';
                  'itctype'       'string'   {'phasecoher','phasecoher2','coher'}  'phasecoher';
                  'nfreqs'        'integer'  [0 Inf]                  [];
                  'lowmem'        'string'   {'on','off'}              'off';
                  'method'        'string'   { 'mod','corrsin','corrcos','latphase' }         'mod';
                  'naccu'         'integer'  [1 Inf]                   250;
                  'newfig'        'string'   {'on','off'}              'on';
                  'padratio'      'integer'  [1 Inf]                   2;
                  'rmerp'         'string'   {'on','off'}              'off';
                  'rboot'         'real'     []                        [];
                  'subitc'        'string'   {'on','off'}              'off';
                  'subwin'        'real'     []                        []; ...
                  'gammapowerlim' 'real'     []                        []; ...
                  'powerlim'      'real'     []                        []; ...
                  'powerlat'      'real'     []                        []; ...
                  'gammabase'     'real'     []                        []; ...
                  'timesout'      'real'     []                        []; ...
                  'ntimesout'     'integer'  []                        200; ...
                  'tlimits'       'real'     []                        [0 frame/srate];
                  'title'         'string'   []                        '';
                  'vert'          { 'real','cell' }  []                [];
                  'wavelet'       'real'     [0 Inf]                   0;
                  'wavelet2'      'real'     [0 Inf]                   [];
                  'winsize'       'integer'  [0 Inf]                   max(pow2(nextpow2(frame)-3),4) }, 'pac');

if ischar(g), error(g); end

% more defaults
% -------------
if isempty(g.wavelet2), g.wavelet2 = g.wavelet; end
if isempty(g.freqs2),   g.freqs2   = g.freqs;   end

% remove ERP if necessary
% -----------------------
X = squeeze(X);
Y = squeeze(Y);X = squeeze(X);
trials = size(X,2);
if strcmpi(g.rmerp, 'on')
    X = X - repmat(mean(X,2), [1 trials]);
    Y = Y - repmat(mean(Y,2), [1 trials]);
end

% perform timefreq decomposition
% ------------------------------
[alltfX freqs1 timesout1] = timefreq(X, srate, 'ntimesout',  g.ntimesout, 'timesout',  g.timesout,  'winsize',  g.winsize, ...
                                'tlimits', g.tlimits, 'detrend',   g.detrend,   'itctype',  g.itctype, ...
                                'subitc',  g.subitc,  'wavelet',   g.wavelet,   'padratio', g.padratio, ...
                                'freqs',   g.freqs,   'freqscale', g.freqscale, 'nfreqs',   g.nfreqs); 
[alltfY freqs2 timesout2] = timefreq(Y, srate, 'ntimesout',  g.ntimesout, 'timesout',  g.timesout,  'winsize',  g.winsize, ...
                                'tlimits', g.tlimits, 'detrend',   g.detrend,   'itctype',  g.itctype, ...
                                'subitc',  g.subitc,  'wavelet',   g.wavelet2,  'padratio', g.padratio, ...
                                'freqs',   g.freqs2,  'freqscale', g.freqscale, 'nfreqs',   g.nfreqs); 

% check time limits
% -----------------
if ~isempty(g.subwin)
    ind1      = find(timesout1 > g.subwin(1) & timesout1 < g.subwin(2));
    ind2      = find(timesout2 > g.subwin(1) & timesout2 < g.subwin(2));
    alltfX    = alltfX(:, ind1, :);
    alltfY    = alltfY(:, ind2, :);
    timesout1 = timesout1(ind1);
    timesout2 = timesout2(ind2);
end
if length(timesout1) ~= length(timesout2) || any( timesout1 ~= timesout2)
    disp('Warning: Time points are different for X and Y. Use ''timesout'' to specify common time points');
    [vals ind1 ind2 ] = intersect_bc(timesout1, timesout2);
    fprintf('Searching for common time points: %d found\n', length(vals));
    if length(vals) < 10, error('Less than 10 common data points'); end
    timesout1 = vals;
    timesout2 = vals;
    alltfX = alltfX(:, ind1, :);
    alltfY = alltfY(:, ind2, :);
end

% scan accross frequency and time
% -------------------------------
%if isempty(g.alpha)
%    disp('Warning: if significance mask is not applied, result might be slightly')
%    disp('different (since angle is not made uniform and amplitude interpolated)')
%end

cohboot =[];
if ~strcmpi(g.method, 'latphase')
    for find1 = 1:length(freqs1)
        for find2 = 1:length(freqs2)           
            for ti = 1:length(timesout1)

                % get data
                % --------
                tmpalltfx = squeeze(alltfX(find1,ti,:));            
                tmpalltfy = squeeze(alltfY(find2,ti,:));

                %if ~isempty(g.alpha)
                %    tmpalltfy = angle(tmpalltfy);
                %    tmpalltfx = abs(  tmpalltfx);
                %    [ tmp cohboot(find1,find2,ti,:) newamp newangle ] = ...
                %        bootcircle(tmpalltfx, tmpalltfy, 'naccu', g.naccu); 
                %    crossfcoh(find1,find2,ti) = sum ( newamp .* exp(j*newangle) );
                %else 
                tmpalltfy = angle(tmpalltfy);
                tmpalltfx = abs(  tmpalltfx);
                if strcmpi(g.method, 'mod')
                    crossfcoh(find1,find2,ti) = sum( tmpalltfx .* exp(j*tmpalltfy) );
                elseif strcmpi(g.method, 'corrsin')
                    tmp = corrcoef( sin(tmpalltfy), tmpalltfx);
                    crossfcoh(find1,find2,ti) = tmp(2);
                else
                    tmp = corrcoef( cos(tmpalltfy), tmpalltfx);
                    crossfcoh(find1,find2,ti) = tmp(2);
                end
            end
        end
    end
elseif 1
    % this option computes power at a given latency
    % then computes the same as above (vectors)
    
    %if isempty(g.powerlat)
    %    error('You need to specify a latency for the ''powerlat'' option');
    %end
        
    gammapower = mean(10*log10(alltfX(:,:,:).*conj(alltfX)),1); % average all frequencies for power
    if isempty(g.gammapowerlim)
        g.gammapowerlim = [ min(gammapower(:)) max(gammapower(:)) ];
    end
    fprintf('Gamma power limits: %3.2f to %3.2f\n', g.gammapowerlim(1), g.gammapowerlim(2)); 
    power = 10*log10(alltfY(:,:,:).*conj(alltfY));
    if isempty(g.powerlim)
        for freq = 1:size(power,1)
            g.powerlim(freq,:) = [ min(power(freq,:)) max(power(freq,:)) ];
        end
    end
    for freq = 1:size(power,1)
        fprintf('Freq %d power limits: %3.2f to %3.2f\n', freqs2(freq), g.powerlim(freq,1), g.powerlim(freq,2)); 
    end
            
    % power plot
    %figure; plot(timesout2/1000, (mean(power(9,:,:),3)-mean(power(9,:)))/50);
    %hold on; plot(linspace(0, length(Y)/srate, length(Y)), mean(Y'), 'g');

    % phase with power
    % figure; plot(timesout2/1000, (mean(phaseangle(9,:,:),3)-mean(phaseangle(9,:)))/50);
    % hold on; plot(timesout1/1000, (mean(gammapower,3)-mean(gammapower(:)))/100, 'r');
    %figure; plot((mean(phaseangle(9,:,:),3)-mean(phaseangle(9,:)))/50+j*(mean(gammapower,3)-mean(gammapower(:)))/100, '.');
    
    matsize               = 32;
    matcenter             = (matsize-1)/2+1;
    matrixfinalgammapower = zeros(size(alltfY,1),size(alltfX,3),matsize,matsize);
    matrixfinalcount      = zeros(size(alltfY,1),size(alltfX,3),matsize,matsize);    
    
    % get power indices
    if isempty(g.gammabase)
        g.gammabase = mean(gammapower(:));
    end
    fprintf('Gamma power average: %3.2f\n', g.gammabase); 
    gammapoweradd  = gammapower-g.gammabase;
    gammapower     = floor((gammapower-g.gammapowerlim(1))/(g.gammapowerlim(2)-g.gammapowerlim(1))*(matsize-2))+1;
    phaseangle     = angle(alltfY);
    posx           = zeros(size(power));
    posy           = zeros(size(power));
    for freq = 1:length(freqs2)
        fprintf('Processing frequency %3.2f\n', freqs2(freq));
        power(freq,:,:) = (power(freq,:,:)-g.powerlim(freq,1))/(g.powerlim(freq,2)-g.powerlim(freq,1))*(matsize-3)/2+1;
        complexval      = power(freq,:,:).*exp(j*phaseangle(freq,:,:));
        posx(freq,:,:)  = round(real(complexval)+matcenter);
        posy(freq,:,:)  = round(imag(complexval)+matcenter);
        for trial = 1:size(alltfX,3) % scan trials
            for time = 1:size(alltfX,2)
                %matrixfinal(freq,posx(freq,time,trial),posy(freq,time,trial),gammapower(1,time,trial)) = ...
                %    matrixfinal(freq,posx(freq,time,trial),posy(freq,time,trial),gammapower(1,time,trial))+1;
                matrixfinalgammapower(freq,trial,posx(freq,time,trial),posy(freq,time,trial)) = ...
                    matrixfinalgammapower(freq,trial,posx(freq,time,trial),posy(freq,time,trial))+gammapoweradd(1,time,trial);
                matrixfinalcount(freq,trial,posx(freq,time,trial),posy(freq,time,trial)) = ...
                    matrixfinalcount(freq,trial,posx(freq,time,trial),posy(freq,time,trial))+1;
            end
        end
        %matrixfinal(freq,:,:,:) = convn(squeeze(matrixfinal(freq,:,:,:)), gs, 'same');
        %tmpmat = posx(index,:)+(posy(index,:)-1)*64+(gammapower(:)-1)*64*64;
        matrixfinalcount(freq, find(matrixfinalcount(freq,:) == 0)) = 1;
        matrixfinalgammapower(freq,:,:,:) = matrixfinalgammapower(freq,:,:,:)./matrixfinalcount(freq,:,:,:);
    end
    
    % average and smooth
    matrixfinalgammapowermean = squeeze(mean(matrixfinalgammapower,2));
    for freq = 1:length(freqs2)
        matrixfinalgammapowermean(freq,:,:) = conv2(squeeze(matrixfinalgammapowermean(freq,:,:)), gauss2d(5,5), 'same');
    end
    %matrixfinalgammapower = matrixfinalgammapower/size(alltfX,3)/size(alltfX,2);
    
    %vect = linspace(-pi,pi,50);    
    %for f = 1:length(freqs2)
    %    crossfcoh(f,:) = hist(tmpalltfy(f,:), vect);
    %end
    
    % smoothing of output image
    % -------------------------
    %gs = gauss2d(6, 6, 6);
    %crossfcoh = convn(crossfcoh, gs, 'same');
    %freqs1    = freqs2;
    %timesout1 = linspace(-180, 180, size(crossfcoh,2));

    crossfcoh    = matrixfinalgammapowermean;
    crossfcohall = matrixfinalgammapower;
     
else
    % this option computes power at a given latency
    % then computes the same as above (vectors)
    
    %if isempty(g.powerlat)
    %    error('You need to specify a latency for the ''powerlat'' option');
    %end
        
    gammapower = mean(10*log10(alltfX(:,:,:).*conj(alltfX)),1); % average all frequencies for power
    if isempty(g.gammapowerlim)
        g.gammapowerlim = [ min(gammapower(:)) max(gammapower(:)) ];
    end
    power = 10*log10(alltfY(:,:,:).*conj(alltfY));
    if isempty(g.powerlim)
        for freq = 1:size(power,1)
            g.powerlim(freq,:) = [ min(power(freq,:)) max(power(freq,:)) ];
        end
    end

    % power plot
    %figure; plot(timesout2/1000, (mean(power(9,:,:),3)-mean(power(9,:)))/50);
    %hold on; plot(linspace(0, length(Y)/srate, length(Y)), mean(Y'), 'g');

    % phase with power
    % figure; plot(timesout2/1000, (mean(phaseangle(9,:,:),3)-mean(phaseangle(9,:)))/50);
    % hold on; plot(timesout1/1000, (mean(gammapower,3)-mean(gammapower(:)))/100, 'r');
    %figure; plot((mean(phaseangle(9,:,:),3)-mean(phaseangle(9,:)))/50+j*(mean(gammapower,3)-mean(gammapower(:)))/100, '.');
    
    matsize               = 64;
    matcenter             = (matsize-1)/2+1;
    matrixfinal           = zeros(size(alltfY,1),64,64,64);
    matrixfinalgammapower = zeros(size(alltfY,1),matsize,matsize);
    matrixfinalcount      = zeros(size(alltfY,1),matsize,matsize);    
    
    % get power indices
    gammapoweradd  = gammapower-mean(gammapower(:));
    gammapower     = floor((gammapower-g.gammapowerlim(1))/(g.gammapowerlim(2)-g.gammapowerlim(1))*(matsize-1))+1;
    phaseangle     = angle(alltfY);
    posx           = zeros(size(power));
    posy           = zeros(size(power));
    gs             = gauss3d(6, 6, 6);
    for freq = 1:size(alltfY)
        fprintf('Processing frequency %3.2f\n', freqs2(freq));
        power(freq,:,:) = (power(freq,:,:)-g.powerlim(freq,1))/(g.powerlim(freq,2)-g.powerlim(freq,1))*(matsize-2)/2;
        complexval      = power(freq,:,:).*exp(j*phaseangle(freq,:,:));
        posx(freq,:,:)  = round(real(complexval)+matcenter);
        posy(freq,:,:)  = round(imag(complexval)+matcenter);
        for trial = 1:size(alltfX,3) % scan trials
            for time = 1:size(alltfX,2)
                %matrixfinal(freq,posx(freq,time,trial),posy(freq,time,trial),gammapower(1,time,trial)) = ...
                %    matrixfinal(freq,posx(freq,time,trial),posy(freq,time,trial),gammapower(1,time,trial))+1;
                matrixfinalgammapower(freq,posx(freq,time,trial),posy(freq,time,trial)) = ...
                    matrixfinalgammapower(freq,posx(freq,time,trial),posy(freq,time,trial))+gammapoweradd(1,time,trial);
                matrixfinalcount(freq,posx(freq,time,trial),posy(freq,time,trial)) = ...
                    matrixfinalcount(freq,posx(freq,time,trial),posy(freq,time,trial))+1;
            end
        end
        %matrixfinal(freq,:,:,:) = convn(squeeze(matrixfinal(freq,:,:,:)), gs, 'same');
        %tmpmat = posx(index,:)+(posy(index,:)-1)*64+(gammapower(:)-1)*64*64;
        matrixfinalcount(freq, find(matrixfinalcount(freq,:) == 0)) = 1;
        matrixfinalgammapower(freq,:,:) = matrixfinalgammapower(freq,:,:)./matrixfinalcount(freq, :,:);
        matrixfinalgammapower(freq,:,:) = conv2(squeeze(matrixfinalgammapower(freq,:,:)), gauss2d(5,5), 'same');
    end
    %matrixfinalgammapower = matrixfinalgammapower/size(alltfX,3)/size(alltfX,2);
    
    %vect = linspace(-pi,pi,50);    
    %for f = 1:length(freqs2)
    %    crossfcoh(f,:) = hist(tmpalltfy(f,:), vect);
    %end
    
    % smoothing of output image
    % -------------------------
    %gs = gauss2d(6, 6, 6);
    %crossfcoh = convn(crossfcoh, gs, 'same');
    %freqs1    = freqs2;
    %timesout1 = linspace(-180, 180, size(crossfcoh,2));

    crossfcoh = matrixfinalgammapower;
    
end

% 7/31/2014 Ramon: crossfcohall sometimes does not exist depending on choice of input options
if ~exist('crossfcohall', 'var')
    crossfcohall = [];
end    


