% crossfreq() - compute cross-frequency coherences. Power of first input
%               correlation with phase of second.
%
% Usage:
%   >> crossfreq(x,y,srate);
%   >> [coh,timesout,freqsout1,freqsout2,cohboot] ...
%                     = crossfreq(x,y,srate,'key1', 'val1', 'key2', val2' ...);
% Inputs:
%    x       = [float array] 2-D data array of size (times,trials) or
%              3-D (1,times,trials)
%    y       = [float array] 2-D or 3-d data array
%    srate   = data sampling rate (Hz)
%
%    Most important optional inputs
%       'mode'      = ['amp_amp'|'amp_phase'|'phase_phase'] correlation mode
%                     is either amplitude-amplitude ('amp_amp'), amplitude
%                     and phase ('amp_phase') and phase-phase ('phase_phase').
%                     Default is 'amp_phase'. 
%       'method'    = ['mod'|'corrsin'|'corrcos'] modulation method ('mod')
%                     or correlation of amplitude with sine or cosine of 
%                     angle (see ref).
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
%       'lowmem'    = ['on'|'off'] compute frequency, by frequency to save
%                     memory. Default 'off'.
%
%    Optional Bootstrap Parameters:
%       'alpha'     = If non-0, compute two-tailed bootstrap significance prob. 
%                     level. Show non-signif. output values as green.   {0}
%       'naccu'     = Number of bootstrap replications to accumulate. Note that
%                     naccu might be automatically increate depending on the
%                     value for 'alpha' {250}
%       'baseboot'  = Bootstrap baseline subtract time window in ms. If only one 
%                     is entered, baseline is from beginning of data to this 
%                     value. Note: you must specify 'tlimits' for bootstrap {0}
%       'boottype'  = ['times'|'timestrials'|'trials'] Bootstrap type: Either
%                     shuffle windows ('times') or windows and trials ('timestrials')
%                     or trials only using a separate bootstrap for each time window
%                     ('trials'). Option 'times' is not recommended but requires less
%                     memory {default 'timestrials'}
%       'rboot'     = Input bootstrap coherence limits (e.g., from crossfreq())
%                     The bootstrap type should be identical to that used
%                     to obtain the input limits. {default: compute from data}
%
%    Optional Plotting Parameters:
%       'title'     = Optional figure title                              {none}
%       'vert'      = [times_vector] plot vertical dashed lines at specified 
%                     times in ms. Can also be a cell array specifying line aspect. 
%                     I.e. { { 0 'color' 'b' 'linewidth' 2 } {1000 'color' 'r' }} 
%                     would draw two lines, one blue thick line at latency 0 and one 
%                     thin red line at latency 1000.
%       'newfig'    = ['on'|'off'] Create new figure for difference plots {'on'}
%                     
% Outputs: 
%        crossfcoh   = Matrix (nfreqs1,nfreqs2,timesout) of coherence (complex).
%                      Use 20*log(abs(crossfcoh)) to vizualize log spectral diffs. 
%        timesout    = Vector of output times (window centers) (ms).
%        freqsout1   = Vector of frequency bin centers for first argument (Hz).
%        freqsout2   = Vector of frequency bin centers for second argument (Hz).
%        cohboot     = Matrix (nfreqs1,nfreqs2,2) of p-value coher signif. 
%                      values. if 'boottype' is 'trials',
%                      (nfreqs1,nfreqs2,timesout,2)
%        alltfX      = single trial spectral decomposition of X
%        alltfY      = single trial spectral decomposition of Y
%
% Author: Arnaud Delorme & Scott Makeig, SCCN/INC, UCSD 2003-
%
% Ref: Testing for Nested Oscilations (2008) J Neuro Methods 174(1):50-61
%
% See also: timefreq(), crossf()

% Copyright (C) 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [crossfcoh, timesout1, freqs1, freqs2, cohboot, alltfX, alltfY] = ...
        crossfreq(X, Y, srate, varargin);
    
if nargin < 1
    help crossfreq; 
    return; 
end;

% deal with 3-D inputs
% --------------------
if ndims(X) == 3, X = reshape(X, size(X,2), size(X,3)); end;
if ndims(Y) == 3, Y = reshape(Y, size(Y,2), size(Y,3)); end;
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
                  'mode'          'string'   { 'amp_amp','amp_phase','phase_phase' } 'amp_phase';
                  'method'        'string'   { 'mod','corrsin','corrcos' }         'mod';
                  'naccu'         'integer'  [1 Inf]                   250;
                  'newfig'        'string'   {'on','off'}              'on';
                  'padratio'      'integer'  [1 Inf]                   2;
                  'rmerp'         'string'   {'on','off'}              'off';
                  'rboot'         'real'     []                        [];
                  'subitc'        'string'   {'on','off'}              'off';
                  'subwin'        'real'     []                        []; ...
                  'timesout'      'real'     []                        []; ...
                  'ntimesout'     'integer'  []                        200; ...
                  'tlimits'       'real'     []                        [0 frame/srate];
                  'title'         'string'   []                        '';
                  'vert'          { 'real','cell' }  []                [];
                  'wavelet'       'real'     [0 Inf]                   0;
                  'wavelet2'      'real'     [0 Inf]                   [];
                  'winsize'       'integer'  [0 Inf]                   max(pow2(nextpow2(frame)-3),4) }, 'crossfreq');

if isstr(g), error(g); end;

% more defaults
% -------------
if isempty(g.wavelet2), g.wavelet2 = g.wavelet; end;
if isempty(g.freqs2),   g.freqs2   = g.freqs;   end;

% remove ERP if necessary
% -----------------------
X = squeeze(X);
Y = squeeze(Y);X = squeeze(X);
trials = size(X,2);
if strcmpi(g.rmerp, 'on')
    X = X - repmat(mean(X,2), [1 trials]);
    Y = Y - repmat(mean(Y,2), [1 trials]);
end;

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
end;
if length(timesout1) ~= length(timesout2) | any( timesout1 ~= timesout2)
    disp('Warning: Time points are different for X and Y. Use ''timesout'' to specify common time points');
    disp('Searching for common points');
    [vals ind1 ind2 ] = intersect(timesout1, timesout2);
    if length(vals) < 10, error('Less than 10 common data points'); end;
    timesout1 = vals;
    timesout2 = vals;
    alltfX = alltfX(:, ind1, :);
    alltfY = alltfY(:, ind2, :);
end;

% scan accross frequency and time
% -------------------------------
if isempty(g.alpha)
    disp('Warning: if significance mask is not applied, result might be slightly')
    disp('different (since angle is not made uniform and amplitude interpolated)')
end;

cohboot =[];
for find1 = 1:length(freqs1)
    for find2 = 1:length(freqs2)           
        for ti = 1:length(timesout1)
            
            % get data
            % --------
            tmpalltfx = squeeze(alltfX(find1,ti,:));            
            tmpalltfy = squeeze(alltfY(find2,ti,:));

            if ~isempty(g.alpha)
                tmpalltfy = angle(tmpalltfy);
                tmpalltfx = abs(  tmpalltfx);
                [ tmp cohboot(find1,find2,ti,:) newamp newangle ] = ...
                    bootcircle(tmpalltfx, tmpalltfy, 'naccu', g.naccu); 
                crossfcoh(find1,find2,ti) = sum ( newamp .* exp(j*newangle) );
            else 
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
                end;
            end;
        end;
    end;
end;

    

