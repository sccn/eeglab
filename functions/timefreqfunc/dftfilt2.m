% dftfilt2() - discrete Fourier filter
%
% Usage:
%   >> wavelet = dftfilt2( freqs, cycles, srate, cyclefact)
%
% Inputs:
%   freqs   - frequency array
%   cycles  - cycles array. If one vale is given, all wavelet use
%             the same number of cycles. If 2 numbers are given, these
%             2 number are used for the number of cycles at the lowest
%             frequency and at the higher frequency. 
%   srate   - sampling rate
%   cycleinc - ['linear'|'log'] increase mode if [min max] cycles
%              provided in 'cycle' parameter. Default is 'log' but
%              use 'linear' if frequencies are 'log' spaced.
%
% Output:
%   wavelet - cell array of wavelet filter to apply onto the data
%             (if 1 frequency only, directelly return the wavelet not
%             a cell array).
%
% Note: The size of the window is automatically computed from the 
%       number of cycles ans is always odd.
%
% Authors: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 3/28/2003

% Copyright (C) 3/28/2003 Arnaud Delorme 8, SCCN/INC/UCSD, arno@sccn.ucsd.edu
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

function wavelet = dftfilt2( freqs, cycles, srate, cycleinc);

    if nargin < 3
        error('3 arguments required');
    end;
    
    % compute number of cycles at each frequency
    % ------------------------------------------
    if length(cycles) == 1
        cycles = cycles*ones(size(freqs));
    elseif length(cycles) == 2
        if nargin == 5 & strcmpi(cycleinc, 'linear') % cycleinc
            cycles = linspace(cycles(1), cycles(2), length(freqs));
        else
            cycles = linspace(log(cycles(1)), log(cycles(2)), length(freqs));
            cycles = exp(cycles);
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
        winval = linspace(-winlenint/2, winlenint/2, winlenint+1);        
        
        win = exp(2i*pi*freqs(index)*winval/srate);
        wavelet{index} = win .* hanning(length(winval))';
        
    end;

    if length(freqs) == 1
        wavelet  = wavelet{1};
    end;
    
    return;
    
    % testing
    % -------
    wav1 = dftfilt2(5, 5, 256);
    abs1 = linspace(-floor(length(wav1)),floor(length(wav1)), length(wav1));
    figure; plot(abs1, real(wav1), 'b');
    
    wav2 = dftfilt2(5, 3, 256);
    abs2 = linspace(-floor(length(wav2)),floor(length(wav2)), length(wav2)); 
    hold on; plot(abs2, real(wav2), 'r');
    
    wav3 = dftfilt2(5, 1.4895990, 256);
    abs3 = linspace(-floor(length(wav3)),floor(length(wav3)), length(wav3)); 
    hold on; plot(abs3, real(wav3), 'g');

    wav4 = dftfilt2(5, 8.73, 256);
    abs4 = linspace(-floor(length(wav4)),floor(length(wav4)), length(wav4)); 
    hold on; plot(abs4, real(wav4), 'm');
    
    