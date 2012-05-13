% eegfiltsig() - (high|low|band)-pass filter using butterworth, Chebychev
%                or ellipsoid filter from the signal processing toolbox
%
% Usage:
%  >> [smoothdata] = eegfiltsig(data,srate,'key', 'val');
%
% Inputs:
%   data        = (channels,frames*epochs) data to filter
%   srate       = data sampling rate (Hz)
%
% Optional inputs
%   'locutoff'    = low-edge frequency in pass band (Hz) 
%   'hicutoff'    = high-edge frequency in pass band (Hz)
%   'lotrans'     = low-edge transition in Hz
%   'hitrans'     = high-edge transition in Hz
%   'filtmode'    = ['pass'|'high'|'low'|'stop'] filter mode. Default is
%                  'pass' for passband. If only 'locutoff' is defined 
%                  then the filter mode is automatically set to 'high'. If
%                  only 'hicutoff' is defined then the filter mode is 
%                  automatically set to 'low'.
%   'filttype'    = ['butter'|'cheby1'|'cheby2'|'ellip'] filter type. Default
%                  is 'butter' for Butterworth
%   'rippledb'    = [float] max ripple in the pass band in dB. Default 2.
%   'stopdb'      = [float] minimum attenuation in the stop band. Default 20.
%
% Outputs:
%    smoothdata = smoothed data
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2012 

% Copyright (C) Arnaud Delorme, SCCN/INC/UCSD, arno@sccn.ucsd.edu
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

function data = eegfiltsig(data, srate, varargin)

if nargin < 1
    help myfilter;
    return;
end;

opt = finputcheck(varargin, { 'rippledb'    'float'   []    2;
                              'stopdb'      'float'   []    20;
                              'locutoff'    'float'   []    [];
                              'lotrans'     'float'   []    [];
                              'hicutoff'    'float'   []    [];
                              'hitrans'     'float'   []    [];                              
                              'filtertype'  'string'  { 'butter' 'elipsoid' 'cheby1' 'cheby2' } 'butter';
                              'filtmode'    'string'  { 'high' 'low' 'pass' 'stop' } 'pass' });                              
if isstr(opt), error(opt), end;

if isempty(opt.hicutoff)
    opt.filtmode = 'high';
end;
if isempty(opt.locutoff)
    opt.filtmode = 'low';
end;
if isempty(opt.lotrans)
    opt.lotrans = opt.locutoff-max(opt.locutoff-3, opt.locutoff/2);
end;
if isempty(opt.hitrans)
    opt.hitrans = 1;
end;

switch(opt.filtertype)
    case 'butter',   ordfunc = @buttord;  filtfunc = @butter;
    case 'cheby1',   ordfunc = @cheb1ord; filtfunc = @cheby1;
    case 'cheby2',   ordfunc = @cheb2ord; filtfunc = @cheby2;
    case 'elipsoid', ordfunc = @ellipord; filtfunc = @ellip;
end;

Wp = [ opt.locutoff/srate opt.hicutoff/srate ];
Ws = [ (opt.locutoff-opt.lotrans)/srate (opt.hicutoff+opt.hitrans)/srate ];
[n,Wn] = ordfunc(Wp,Ws,opt.rippledb,opt.stopdb);

if strcmpi(opt.filtmode, 'pass'), [b,a] = filtfunc(n,Wn);
else                              [b,a] = filtfunc(n,Wn, opt.filtmode);
end;

fprintf('Filter order %d (ripple max %3.2f dB, stop band %3.2f dB)\n', n, opt.rippledb, opt.stopdb);
if strcmpi(opt.filtmode, 'high')
    fprintf('High pass transition band: %3.2f Hz - %3.2f Hz\n', Ws(1)*srate, Wp(1)*srate);
end;
if strcmpi(opt.filtmode, 'low')
    fprintf('Low pass transition band: %3.2f Hz - %3.2f Hz\n', Ws(end)*srate, Wp(end)*srate);
end;
freqz(b,a,2000,256);

fprintf('Filtering:');
for index = 1:size(data,1)
    fprintf('.');
    data(index,:) = filtfilt(b, a, data(index,:));
end;
fprintf('\n');
return;
    
    
    
Ws = [max(locutoff-3, locutoff/2)/srate (hicutoff+1)/srate ];
Wp = [ locutoff/srate hicutoff/srate ];
[n,Wn] = buttord(Wp,Ws,rippledb,stopdb);
[b,a] = butter(n,Wn, filtmode{:});
fprintf('Filter order %d (ripple max %3.2f dB, stop band %3.2f dB)\n', n, rippledb, stopdb);
if isempty(filtmode) || strcmpi(filtmode{1}, 'high')
    fprintf('High pass transition band: %3.2f Hz - %3.2f Hz\n', Ws(1)*srate, Wp(1)*srate);
end;
if isempty(filtmode) || strcmpi(filtmode{1}, 'low')
    fprintf('Low pass transition band: %3.2f Hz - %3.2f Hz\n', Ws(end)*srate, Wp(end)*srate);
end;
freqz(b,a,2000,256);

fprintf('Filtering:');
for index = 1:size(data,1)
    fprintf('.');
    data(index,:) = filtfilt(b, a, data(index,:));
end;
fprintf('\n');
