% compute_ERSP_times() - computes the maximum time window of the ERSP,  
%        based on the frequency range requested, and other input parameters 
%        that are used by timef(). 
%        This is a helper function called from pop_preclust() & cls_ersp(). 

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