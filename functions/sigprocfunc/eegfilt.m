% eegfilt() - (high|low|band)-pass filter EEG data using two-way least-squares 
%              FIR filtering. Multiple data channels and epochs supported.
%              (Requires the MATLAB Signal Processing Toolbox)
% Usage:
%  >> [smoothdata] = eegfilt(data,srate,locutoff,hicutoff);
%  >> [smoothdata,filtwts] = eegfilt(data,srate,locutoff,hicutoff, ...
%                                             epochframes,filtorder);
% Inputs:
%   data        = (channels,frames*epochs) data to filter
%   srate       = data sampling rate (Hz)
%   locutoff    = low edge frequency in pass band (Hz)  {0 -> lowpass}
%   hicutoff    = high edge frequency in pass band (Hz) {0 -> highpass}
%   epochframes = frames per epoch (filter each epoch separately {def/0: all}
%   filtorder   = length of the filter in points {default 3*fix(srate/locutoff)}
%   revfilt     = [0|1] reverse filter (i.e. bandpass filter to notch filter). {0}
%
% Outputs:
%    smoothdata = smoothed data
%    filtwts    = filter coefficients [smoothdata <- filtfilt(filtwts,1,data)]
%
% Author: Scott Makeig, Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 1997 

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

% $Log: not supported by cvs2svn $
% Revision 1.8  2003/01/24 00:23:33  arno
% print information about transition bands
%
% Revision 1.7  2003/01/23 23:53:25  arno
% change text
%
% Revision 1.6  2003/01/23 23:47:26  arno
% same
%
% Revision 1.5  2003/01/23 23:40:59  arno
% implementing notch filter
%
% Revision 1.4  2002/08/09 14:55:36  arno
% update transition band
%
% Revision 1.3  2002/08/09 02:06:27  arno
% updating transition
%
% Revision 1.2  2002/08/09 02:04:34  arno
% debugging
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 5-08-97 fixed frequency bound computation -sm
% 10-22-97 added MINFREQ tests -sm
% 12-05-00 added error() calls -sm
% 01-25-02 reformated help & license, added links -ad 

function [smoothdata,filtwts] = eegfilt(data,srate,locutoff,hicutoff,epochframes,filtorder, revfilt)

if nargin<4
    fprintf('');
    help eegfilt
    return
end

if ~exist('firls')
   error('*** eegfilt() requires the signal processing toolbox. ***');
end

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
    help eegfilt
    error('locutoff > hicutoff ???\n');
end
if locutoff < 0 | hicutoff < 0,
   help eegfilt
   error('locutoff | hicutoff < 0 ???\n');
end

if locutoff>nyq,
    help eegfilt
    error('locutoff cannot be > srate/2');
end

if hicutoff>=nyq
   hicutoff = 0; 
end

if nargin<6
   filtorder = 0;
end
if nargin<7
   revfilt = 0;
end

if isempty(filtorder) | filtorder==0,
   if locutoff>0,
     filtorder = minfac*fix(srate/locutoff);
   elseif hicutoff>0,
     filtorder = minfac*fix(srate/hicutoff);
   end
     
   if filtorder < min_filtorder
        filtorder == min_filtorder;
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
    fprintf('filter order is %d. ',filtorder);
    error('epochframes must be 3 times the filtorder.');
end
if (1+trans)*hicutoff/nyq > 1
	error('high cutoff frequency to close to Nyquist frequency');
end;

if locutoff > 0 & hicutoff > 0,    % bandpass filter
    if revfilt
         fprintf('eegfilt() - performing %d-point notch filtering.\n',filtorder);
    else fprintf('eegfilt() - performing %d-point bandpass filtering.\n',filtorder);
    end; 
    f=[MINFREQ (1-trans)*locutoff/nyq locutoff/nyq hicutoff/nyq (1+trans)*hicutoff/nyq 1]; 
    fprintf('eegfilter() - low transition band is %1.1g Hz; high one is %1.1g Hz.\n',(f(3)-f(2))*srate, (f(5)-f(4))*srate);
    m=[0       0                      1            1            0                      0]; 
elseif locutoff > 0                % highpass filter
 if locutoff/nyq < MINFREQ
    help eegfilt
    fprintf('eegfilt() - highpass cutoff freq must be > %g Hz\n\n',MINFREQ*nyq);
    return
 end
 fprintf('eegfilter() - performing %d-point highpass filtering.\n',filtorder);
 f=[MINFREQ locutoff*(1-trans)/nyq locutoff/nyq 1]; 
 fprintf('eegfilter() - highpass transition band is %1.1g Hz.\n',(f(3)-f(2))*srate);
 m=[   0             0                   1      1];
elseif hicutoff > 0                %  lowpass filter
 if hicutoff/nyq < MINFREQ
    help eegfilt
    fprintf('eegfilt() - lowpass cutoff freq must be > %g Hz\n\n',MINFREQ*nyq);
    return
 end
 fprintf('eegfilt() - performing %d-point lowpass filtering.\n',filtorder);
 f=[MINFREQ hicutoff/nyq hicutoff*(1+trans)/nyq 1]; 
 fprintf('eegfilter() - lowpass transition band is %1.1g Hz.\n',(f(3)-f(2))*srate);
 m=[     1           1              0                 0];
else
 help eegfilt
 return
end
if revfilt
    m = ~m;
end;

filtwts = firls(filtorder,f,m); % get FIR filter coefficients

smoothdata = zeros(chans,frames);
for e = 1:epochs                % filter each epoch, channel 
    for c=1:chans
      smoothdata(c,(e-1)*epochframes+1:e*epochframes) ...
        = filtfilt(filtwts,1,data(c,(e-1)*epochframes+1:e*epochframes));
    end
end

