% compute_ERSP_times() - computes the widest possible ERSP/ITC time window,   
%        which depends on requested ERSP/ITC parameters such as epoch limits, 
%        frequency range, wavelet parameters, sampling rate and frequency 
%        resolution that are used by timef(). 
%        This helper function is called by pop_preclust() & std_ersp(). 
% Example:
%    [time_range, winsize] = compute_ersp_times(cycles,  ALLEEG(seti).srate, ...
%                              [ALLEEG(seti).xmin ALLEEG(seti).xmax]*1000, freq(1),padratio);
%
% Authors: Hilit Serby & Arnaud Delorme, SCCN, INC, UCSD, Feb 03, 2005

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, Feb 03, 2005, hilit@sccn.ucsd.edu
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

function [time_range, winsize] = compute_ERSP_times(cycles, srate, epoch_lim, lowfreq, padratio) 

if cycles == 0 %FFT option
    if ~exist('padratio')
        error('You must enter padratio value for FFT ERSP');
    end
    lowfreq = lowfreq*padratio;
    t = 1/lowfreq;%time window in sec
    winsize = t*srate;%time window in points
    %time window in points (must be power of 2) for FFT
    winsize =pow2(nextpow2(winsize));
    %winsize =2^round(log2(winsize)); 
else %wavelet
    t = cycles(1)/lowfreq; %time window in sec
    winsize  = round(t*srate); %time window in points
end

time_range(1) = epoch_lim(1) + .5*t*1000;
time_range(2) = epoch_lim(2) - .5*t*1000;
