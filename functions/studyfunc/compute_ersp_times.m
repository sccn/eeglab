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
