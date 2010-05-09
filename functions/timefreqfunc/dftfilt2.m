% dftfilt2() - discrete complex wavelet filters
%
% Usage:
%   >> wavelet = dftfilt2( freqs, cycles, srate, cyclefact)
%
% Inputs:
%   freqs    - frequency array
%   cycles   - cycles array. If one value is given, all wavelets have
%              the same number of cycles. If two values are given, the
%              two values are used for the number of cycles at the lowest
%              frequency and at the highest frequency, with linear 
%              interpolation between these values for intermediate
%              frequencies
%   srate    - sampling rate (in Hz)
%
%   cycleinc - ['linear'|'log'] increase mode if [min max] cycles is
%              provided in 'cycle' parameter. {default: 'linear'}
%   type     - ['sinus'|'morlet'] wavelet type is a sinusoid with
%              cosine (real) and sine (imaginary) parts tapered by
%              a Hanning or Morlet function. 'morlet' is a typical Morlet 
%              wavelet (with p=2*pi and sigma=0.7) best matching the 
%              'sinus' Hanning taper) {default: 'morlet'}
% Output:
%   wavelet - cell array of wavelet filters
%
% Note: The length of the window is automatically computed from the 
%       number of cycles and is always made odd.
%
% Authors: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 3/28/2003

% Copyright (C) 3/28/2003 Arnaud Delorme 8, SCCN/INC/UCSD, arno@salk.edu
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

function wavelet = dftfilt2( freqs, cycles, srate, cycleinc, type);

    if nargin < 3
        error('3 arguments required');
    end;
    if nargin < 5
        type = 'morlet';
    end;

    % compute number of cycles at each frequency
    % ------------------------------------------
    if length(cycles) == 1
        cycles = cycles*ones(size(freqs));
    elseif length(cycles) == 2
        if nargin == 4 & strcmpi(cycleinc, 'log') % cycleinc
            cycles = linspace(log(cycles(1)), log(cycles(2)), length(freqs));
            cycles = exp(cycles);
        else
            cycles = linspace(cycles(1), cycles(2), length(freqs));
        end;
    end;
    
    % compute wavelet
    for index = 1:length(freqs)
        
        % number of cycles depend on window size 
        % number of cycles automatically reduced if smaller window
        % note: as the number of cycle changes, the frequency shifts a little
        %       this has to be fixed
        
        winlen = cycles(index)*srate/freqs(index);
        winlenint = floor(winlen);
        if mod(winlenint,2) == 1, winlenint = winlenint+1; end; 
        winval = linspace(winlenint/2, -winlenint/2, winlenint+1);        
        
        if strcmpi(type, 'sinus') % Hanning
            win = exp(2i*pi*freqs(index)*winval/srate);
            wavelet{index} = win .* hanning(length(winval))';

        else % Morlet
            t = freqs(index)*winval/srate;
            p = 2*pi;
            s = cycles(index)/5;
            wavelet{index} = exp(j*t*p)/sqrt(2*pi) .* ...
                (exp(-t.^2/(2*s^2))-sqrt(2)*exp(-t.^2/(s^2)-p^2*s^2/4));
        end;    
    end;
    
    
    return;
    
    % testing
    % -------
    wav1 = dftfilt2(5, 5, 256); wav1 = wav1{1};
    abs1 = linspace(-floor(length(wav1)),floor(length(wav1)), length(wav1));
    figure; plot(abs1, real(wav1), 'b');
    
    wav2 = dftfilt2(5, 3, 256); wav2 = wav2{1};
    abs2 = linspace(-floor(length(wav2)),floor(length(wav2)), length(wav2)); 
    hold on; plot(abs2, real(wav2), 'r');
    
    wav3 = dftfilt2(5, 1.4895990, 256); wav3 = wav3{1};
    abs3 = linspace(-floor(length(wav3)),floor(length(wav3)), length(wav3)); 
    hold on; plot(abs3, real(wav3), 'g');

    wav4 = dftfilt2(5, 8.73, 256); wav4 = wav4{1};
    abs4 = linspace(-floor(length(wav4)),floor(length(wav4)), length(wav4)); 
    hold on; plot(abs4, real(wav4), 'm');
    
    % more testing
    % ------------
    freqs = exp(linspace(0,log(10),10));
    win = dftfilt2(freqs, [3 3], 256, 'linear', 'sinus');
    win = dftfilt2(freqs, [3 3], 256, 'linear', 'morlet');
    
    freqs = [12.0008   13.2675   14.5341   15.8007   17.0674   18.3340   19.6007   20.8673   22.1339   23.4006   24.6672   25.9339   27.2005 28.4671   29.7338   31.0004   32.2670   33.5337   34.8003   36.0670   37.3336   38.6002   39.8669   41.1335   42.4002   43.6668 44.9334   46.2001   47.4667 ...
             48.7334   50.0000];
    
    win = dftfilt2(freqs, [3 12], 256, 'linear'); size(win)
    
    winsize = 0;
    for index = 1:length(win)
        winsize = max(winsize,length(win{index}));
    end;
    allwav = zeros(winsize, length(win));
    for index = 1:length(win)
        wav1 = win{index};
        abs1 = linspace(-(length(wav1)-1)/2,(length(wav1)-1)/2, length(wav1));
        allwav(abs1+(winsize+1)/2,index) = wav1(:);
    end;
    figure; imagesc(imag(allwav));


% syemtric hanning function
function w = hanning(n)
if ~rem(n,2)
   w = .5*(1 - cos(2*pi*(1:n/2)'/(n+1)));
   w = [w; w(end:-1:1)];
else
   w = .5*(1 - cos(2*pi*(1:(n+1)/2)'/(n+1)));
   w = [w; w(end-1:-1:1)];
end
