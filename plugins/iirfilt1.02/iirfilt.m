% iirfilt() -  (high|low|band)-pass filter data using an elliptic IIR filter 
%              Multiple data channels and epochs are supported. Forward and
%              reverse filtering are used to avoid phase distortions.
% Usage:
%  >> [smoothdata] = iirfilt(data,srate,locutoff,hicutoff);
%  >> [smoothdata,filtwts] = iirfilt(data,srate,locutoff,hicutoff, ...
%                                    epochframes, trans_bw, revfilt,rp,rs)
% Inputs:
%   data        = (channels,frames*epochs) data to filter
%   srate       = data sampling rate (Hz)
%   locutoff    = low-edge frequency in pass band (Hz).  If 0, lowpass only.
%   hicutoff    = high-edge frequency in pass band (Hz). If 0, highpass only.
%   epochframes = frames per epoch (filter each epoch separately)
%                 {default|0: data is 1 epoch}
%   trans_bw    = width (in Hz) of transition interval between stop and 
%                 pass bands {default: 1 Hz, or if pass band f < 5 Hz, f/3 Hz}.
%                 Enter 0 to use default.
%   revfilt     = [0|1] reverse filter (i.e. bandpass filter to notch 
%                 or band-reject filter). {0}
%   rp          = ripple amplitude in dB in the pass band {default: 0.0025 dB, 
%                 or if pass band f < 5 Hz, 0.01 dB}
%   rs          = ripple amplitude in dB in the stop band {default: 40 dB,
%                 or if pass band f < 5 Hz, 30 dB}
%   causal      = ['on'|'off'] use causal filter. Default 'off'.
%
% Note: Requires the MATLAB Signal Processing Toolbox.
%
% Authors: Maksym Pozdin (mpozdin.ece04@gtalumni.org, IOL/ONRC,2004), 
%          with Arnaud Delorme and Scott Makeig (SCCN/INC/UCSD, La Jolla CA)
%
% See also: eegfilt(), eegfiltfft()

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

% Notes:
% Currently, with cutoff less then 5 Hz, both HPF and LPF use more 'relaxed' 
% filter parameters to achieve best results. 

function [smoothdata,b,a] = iirfilt(data,srate,locutoff,hicutoff,epochframes, trans_bw, revfilt, rp, rs, causal)

if nargin<4
    fprintf('');
    help iirfilt
    return
end

a = 0;
b = 0;
data = double(data);

if exist('ellipord') ~= 2 | exist('ellip') ~= 2
   error('*** ellip() requires the Matlab Signal Processing toolbox. ***');
end

[chans frames] = size(data);
if chans > 1 & frames == 1,
    help iirfilt
    error('input data should be a row vector.');
end
nyq            = srate*0.5;  % Nyquist frequency
MINFREQ = 5;

if locutoff>0 & hicutoff > 0 & locutoff > hicutoff,
    error('locutoff > hicutoff ???');
end
if locutoff < 0 | hicutoff < 0,
   error('locutoff | hicutoff < 0 ???');
end


if locutoff>nyq,
    error('locutoff cannot be > srate/2');
end

if hicutoff>=nyq
   hicutoff = 0; 
end
if nargin < 6
    trans_bw = 0;
end;

if trans_bw==0

    if locutoff < MINFREQ & locutoff > 0
        trans_bw=locutoff/3;
    elseif hicutoff < MINFREQ & hicutoff > 0
        trans_bw=hicutoff/3;
        %rp=0.01;
        %rs=30;
    else
        % Default transition bandwidth is 1 Hz
        trans_bw=1;
    end;
end

if nargin<7
   revfilt = 0;
end
if nargin<8 || isempty(rp)
   % Ripple in the passband
   rp=0.0025;
end
if nargin<9 || isempty(rs)
   % Ripple in the stopband
   rs=40;
end
if nargin<10
   causal = 'off';
end

if nargin<5
        epochframes = 0;
end
if epochframes ==0,
    epochframes = frames;    % default
end
epochs = fix(frames/epochframes);
if epochs*epochframes ~= frames,
    error('epochframes does not divide frames.\n');
end

if locutoff > 0 & hicutoff ==0 & revfilt==0              % highpass filter

    ws=(locutoff-trans_bw)/nyq;
    wp=(locutoff)/nyq;
    [N,wn] = ellipord(wp,ws,rp,rs);
    fprintf('HPF has been designed.\n');
    fprintf('HPF has cutoff of %1.1f Hz, transition bandwidth of %1.1f Hz and its order is %1.1f\n',locutoff, trans_bw,N);
    [b,a]=ellip(N,rp,rs,wn,'high');
    
elseif hicutoff > 0 & locutoff==0 & revfilt==0                %  lowpass filter

    ws=(hicutoff+trans_bw)/nyq;
    wp=(hicutoff)/nyq;
    [N,wn] = ellipord(wp,ws,rp,rs);
    fprintf('LPF has been designed.\n');
    fprintf('LPF has cutoff of %1.1f Hz, transition bandwidth of %1.1f Hz and its order is %1.1f\n',hicutoff, trans_bw,N);
    [b,a]=ellip(N,rp,rs,wn);
    
elseif hicutoff > 0 & locutoff> 0 & revfilt==0               %  bandpass filter

    %first part LPF
    ws=(hicutoff+trans_bw)/nyq;
    wp=(hicutoff)/nyq;
    [N,wn] = ellipord(wp,ws,rp,rs);
    fprintf('BPF has been designed.\n');
    fprintf('LPF has cutoff of %1.1f Hz, transition bandwidth of %1.1f Hz and its order is %1.1f\n',hicutoff, trans_bw,N);
    %[zl, pl, kl]=ellip(N,rp,rs,wn);
    %[lpf_sos,lpf_g] = zp2sos(zl,pl, kl);
    [bl,al]=ellip(N,rp,rs,wn);
   
    %second part HPF
    ws=(locutoff-trans_bw)/nyq;
    wp=(locutoff)/nyq;
    [N,wn] = ellipord(wp,ws,rp,rs);
    fprintf('HPF has cutoff of %1.1f Hz, transition bandwidth of %1.1f Hz and its order is %1.1f\n',locutoff, trans_bw,N);
    %[zh, ph, kh]=ellip(N,rp,rs,wn);
    %[hpf_sos,hpf_g] = zp2sos(zh,ph, kh);
    [bh,ah]=ellip(N,rp,rs,wn,'high');
    %help fvtool
    %b=conv(bh,bl);a=conv(ah,al);
    b.bl=bl;b.bh=bh; a.al=al;a.ah=ah;
    
elseif hicutoff > 0 & locutoff> 0 & revfilt==1               %  bandreject filter

    %first part LPF
    ws=(locutoff+trans_bw)/nyq;
    wp=(locutoff)/nyq;
    [N,wn] = ellipord(wp,ws,rp,rs);
    fprintf('BRF has been designed.\n');
    fprintf('LPF has cutoff of %1.1f Hz, transition bandwidth of %1.1f Hz and its order is %1.1f\n',locutoff, trans_bw,N);
    [bl,al]=ellip(N,rp,rs,wn);
   
    %second part HPF
    ws=(hicutoff-trans_bw)/nyq;
    wp=(hicutoff)/nyq;
    [N,wn] = ellipord(wp,ws,rp,rs);
    fprintf('HPF has cutoff of %1.1f Hz, transition bandwidth of %1.1f Hz and its order is %1.1f\n',hicutoff, trans_bw,N);
    [bh,ah]=ellip(N,rp,rs,wn,'high');
    b.bl=bl;b.bh=bh; a.al=al;a.ah=ah;
    
else       
    error('You must provide a non-0 low or high cut-off frequency');
end


smoothdata = zeros(chans,frames);
lastwarn('');
for e = 1:epochs                % filter each epoch, channel
    for c=1:chans
        if isstruct(a) & isstruct(b) & revfilt==0            %BPF - filter with LPF and HPF in series
            
            if strcmpi(causal, 'on')
                smoothdata1(c,(e-1)*epochframes+1:e*epochframes) = filter(b.bl,a.al,data(c,(e-1)*epochframes+1:e*epochframes));
                smoothdata( c,(e-1)*epochframes+1:e*epochframes) = filter(b.bh,a.ah,smoothdata1(c,(e-1)*epochframes+1:e*epochframes));
            else
                smoothdata1(c,(e-1)*epochframes+1:e*epochframes) = filtfilt(b.bl,a.al,data(c,(e-1)*epochframes+1:e*epochframes));
                smoothdata( c,(e-1)*epochframes+1:e*epochframes) = filtfilt(b.bh,a.ah,smoothdata1(c,(e-1)*epochframes+1:e*epochframes));
            end;
            
        elseif isstruct(a) & isstruct(b) & revfilt==1         %BRF - filter with LPF and HPF in parallel
            if strcmpi(causal, 'on')
                smoothdata1(c,(e-1)*epochframes+1:e*epochframes) = filter(b.bl,a.al,data(c,(e-1)*epochframes+1:e*epochframes));
                smoothdata2(c,(e-1)*epochframes+1:e*epochframes) = filter(b.bh,a.ah,data(c,(e-1)*epochframes+1:e*epochframes));
            else
                smoothdata1(c,(e-1)*epochframes+1:e*epochframes) = filtfilt(b.bl,a.al,data(c,(e-1)*epochframes+1:e*epochframes));
                smoothdata2(c,(e-1)*epochframes+1:e*epochframes) = filtfilt(b.bh,a.ah,data(c,(e-1)*epochframes+1:e*epochframes));
            end;
            smoothdata=smoothdata1+smoothdata2;               %combine final results
        else
            if strcmpi(causal, 'on')
                smoothdata(c,(e-1)*epochframes+1:e*epochframes) = filter(b,a,data(c,(e-1)*epochframes+1:e*epochframes));
            else
                smoothdata(c,(e-1)*epochframes+1:e*epochframes) = filtfilt(b,a,data(c,(e-1)*epochframes+1:e*epochframes));
            end;
        end;
        if epochs == 1
            if rem(c,20) ~= 0
                fprintf('.');
            else
                fprintf('%d',c);
            end
        end
    end
end
fprintf('\n');
[LASTMSG, LASTID] = lastwarn;
if ~isempty(LASTMSG)
    disp('Warning: the warning message (for example "matrix close to singular")');
    disp('         indicates that some of your data segment are too small to be');
    disp('         filtered or that you low edge frequency is too low.');
    disp('         Check the data by looking at it (raw data and data spectrum).');
    disp('         If necessary, reload data and modify the filter');
end;

%% To test my filter: draws Frequency response
% N=10001; n=1:N;
% xx=zeros(size(n));
% xx(ceil(length(xx)/2))=1;
% if isstruct(a) & isstruct(b)
%     yy1=filtfilt(b.bl,a.al,xx);
%     yy=filtfilt(b.bh, a.ah,yy1);
% else
%     yy=filtfilt(b,a,xx);
% end;
% YY=fft(yy,100000);
% ff=-nyq:2*nyq/length(YY):nyq-2*nyq/length(YY);
% figure;plot(ff,abs(fftshift(YY)));grid on;title('Frequency response');
