% compute_ERSP_times() - computes the widest possible ERSP/ITC time window,   
%        which depends on requested ERSP/ITC parameters such as epoch limits, 
%        frequency range, wavelet parameters, sampling rate and frequency 
%        resolution that are used by timef(). 
%        This helper function is called by pop_preclust() & cls_ersp(). 

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
    winsize  = t*srate; %time window in points
end

time_range(1) = epoch_lim(1) + .5*t*1000;
time_range(2) = epoch_lim(2) - .5*t*1000;