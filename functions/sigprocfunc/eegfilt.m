% eegfilt() -  (high|low|band)-pass filter data using two-way least-squares 
%              FIR filtering. Optionally uses the window method instead of 
%              least-squares. Multiple data channels and epochs supported.
%              Requires the MATLAB Signal Processing Toolbox.
% Usage:
%  >> [smoothdata] = eegfilt(data,srate,locutoff,hicutoff);
%  >> [smoothdata,filtwts] = eegfilt(data,srate,locutoff,hicutoff, ...
%                                    epochframes,filtorder,revfilt,firtype,causal);
% Inputs:
%   data        = (channels,frames*epochs) data to filter
%   srate       = data sampling rate (Hz)
%   locutoff    = low-edge frequency in pass band (Hz)  {0 -> lowpass}
%   hicutoff    = high-edge frequency in pass band (Hz) {0 -> highpass}
%   epochframes = frames per epoch (filter each epoch separately {def/0: data is 1 epoch}
%   filtorder   = length of the filter in points {default 3*fix(srate/locutoff)}
%   revfilt     = [0|1] reverse filter (i.e. bandpass filter to notch filter). {default 0}
%   firtype     = 'firls'|'fir1' {'firls'}
%   causal      = [0|1] use causal filter if set to 1 (default 0)
%
% Outputs:
%    smoothdata = smoothed data
%    filtwts    = filter coefficients [smoothdata <- filtfilt(filtwts,1,data)]
%
% See also: firls(), filtfilt()

% Author: Scott Makeig, Arnaud Delorme, Clemens Brunner SCCN/INC/UCSD, La Jolla, 1997

% Copyright (C) 4-22-97 from bandpass.m Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 05-08-97 fixed frequency bound computation -sm
% 10-22-97 added MINFREQ tests -sm
% 12-05-00 added error() calls -sm
% 01-25-02 reformated help & license, added links -ad
% 03-20-12 added firtype option -cb

function [smoothdata,filtwts] = eegfilt(data,srate,locutoff,hicutoff,epochframes,filtorder,revfilt,firtype,causal)

if nargin<4
    fprintf('');
    help eegfilt
    return
end

%if ~exist('firls')
%   error('*** eegfilt() requires the signal processing toolbox. ***');
%end

[chans frames] = size(data);
if chans > 1 & frames == 1,
    help eegfilt
    error('input data should be a row vector.');
end
nyq            = srate*0.5;  % Nyquist frequency
%MINFREQ = 0.1/nyq;
MINFREQ = 0;

minfac         = 3;    % this many (lo)cutoff-freq cycles in filter
min_filtorder  = 15;   % minimum filter length
trans          = 0.15; % fractional width of transition zones

if locutoff>0 & hicutoff > 0 & locutoff > hicutoff,
    error('locutoff > hicutoff ???\n');
end
if locutoff < 0 | hicutoff < 0,
    error('locutoff | hicutoff < 0 ???\n');
end

if locutoff>nyq,
    error('Low cutoff frequency cannot be > srate/2');
end

if hicutoff>nyq
    error('High cutoff frequency cannot be > srate/2');
end

if nargin<6
    filtorder = 0;
end
if nargin<7
    revfilt = 0;
end
if nargin<8
    firtype = 'firls';
end
if nargin<9
    causal = 0;
end

if strcmp(firtype, 'firls')
    warning('Using firls to estimate filter coefficients. We recommend that you use fir1 instead, which yields larger attenuation. In future, fir1 will be used by default!');
end

if isempty(filtorder) | filtorder==0,
    if locutoff>0,
        filtorder = minfac*fix(srate/locutoff);
    elseif hicutoff>0,
        filtorder = minfac*fix(srate/hicutoff);
    end
    
    if filtorder < min_filtorder
        filtorder = min_filtorder;
    end
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

if filtorder*3 > epochframes,   % Matlab filtfilt() restriction
    fprintf('eegfilt(): filter order is %d. ',filtorder);
    error('epochframes must be at least 3 times the filtorder.');
end
if (1+trans)*hicutoff/nyq > 1
    error('high cutoff frequency too close to Nyquist frequency');
end;

if locutoff > 0 & hicutoff > 0,    % bandpass filter
    if revfilt
        fprintf('eegfilt() - performing %d-point notch filtering.\n',filtorder);
    else
        fprintf('eegfilt() - performing %d-point bandpass filtering.\n',filtorder);
    end;
    fprintf('            If a message, ''Matrix is close to singular or badly scaled,'' appears,\n');
    fprintf('            then Matlab has failed to design a good filter. As a workaround, \n');
    fprintf('            for band-pass filtering, first highpass the data, then lowpass it.\n');
    if strcmp(firtype, 'firls')
        f=[MINFREQ (1-trans)*locutoff/nyq locutoff/nyq hicutoff/nyq (1+trans)*hicutoff/nyq 1];
        fprintf('eegfilt() - low transition band width is %1.1g Hz; high trans. band width, %1.1g Hz.\n',(f(3)-f(2))*srate/2, (f(5)-f(4))*srate/2);
        m=[0       0                      1            1            0                      0];
    elseif strcmp(firtype, 'fir1')
        filtwts = fir1(filtorder, [locutoff, hicutoff]./(srate/2));
    end
elseif locutoff > 0                % highpass filter
    if locutoff/nyq < MINFREQ
        error(sprintf('eegfilt() - highpass cutoff freq must be > %g Hz\n\n',MINFREQ*nyq));
    end
    fprintf('eegfilt() - performing %d-point highpass filtering.\n',filtorder);
    if strcmp(firtype, 'firls')
        f=[MINFREQ locutoff*(1-trans)/nyq locutoff/nyq 1];
        fprintf('eegfilt() - highpass transition band width is %1.1g Hz.\n',(f(3)-f(2))*srate/2);
        m=[   0             0                   1      1];
    elseif strcmp(firtype, 'fir1')
        filtwts = fir1(filtorder, locutoff./(srate/2), 'high');
    end
elseif hicutoff > 0                %  lowpass filter
    if hicutoff/nyq < MINFREQ
        error(sprintf('eegfilt() - lowpass cutoff freq must be > %g Hz',MINFREQ*nyq));
    end
    fprintf('eegfilt() - performing %d-point lowpass filtering.\n',filtorder);
    if strcmp(firtype, 'firls')
        f=[MINFREQ hicutoff/nyq hicutoff*(1+trans)/nyq 1];
        fprintf('eegfilt() - lowpass transition band width is %1.1g Hz.\n',(f(3)-f(2))*srate/2);
        m=[     1           1              0                 0];
    elseif strcmp(firtype, 'fir1')
        filtwts = fir1(filtorder, hicutoff./(srate/2));
    end
else
    error('You must provide a non-0 low or high cut-off frequency');
end
if revfilt
    if strcmp(firtype, 'fir1')
        error('Cannot reverse filter using ''fir1'' option');
    else
        m = ~m;
    end;
end;

if strcmp(firtype, 'firls')
    filtwts = firls(filtorder,f,m); % get FIR filter coefficients
end

smoothdata = zeros(chans,frames);
for e = 1:epochs                % filter each epoch, channel
    for c=1:chans
        try
            if causal
                 smoothdata(c,(e-1)*epochframes+1:e*epochframes) = filter(  filtwts,1,data(c,(e-1)*epochframes+1:e*epochframes));
            else smoothdata(c,(e-1)*epochframes+1:e*epochframes) = filtfilt(filtwts,1,data(c,(e-1)*epochframes+1:e*epochframes));
            end;
        catch,
            if causal
                 smoothdata(c,(e-1)*epochframes+1:e*epochframes) = filter(  filtwts,1,double(data(c,(e-1)*epochframes+1:e*epochframes)));
            else smoothdata(c,(e-1)*epochframes+1:e*epochframes) = filtfilt(filtwts,1,double(data(c,(e-1)*epochframes+1:e*epochframes)));
            end;
        end;
        if epochs == 1
            if rem(c,20) ~= 0, fprintf('.');
            else               fprintf('%d',c);
            end
        end
    end
end
fprintf('\n');
