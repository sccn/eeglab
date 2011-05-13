% std_topomovie - make movie in the frequency domain
%
% Usage:    
%  >> [STUDY] = std_movie(STUDY, ALLEEG, key1, val1, key2, val2, ...);  
%
% Inputs:
%   STUDY      - STUDY structure comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - vector of EEG dataset structures for the dataset(s) in the STUDY, 
%                typically created using load_ALLEEG().  
%   'channels' - [numeric vector]  specific channel group to plot. By
%                default, the grand mean channel spectrum is plotted (using the 
%                same format as for the cluster component means described above)
%   'moviemode' - ['erp'|'spec'|'ersptime'] movie mode. Currently only
%                 'spec' is implemented.
%   'erspfreq'  - [min max] frequency range when making movie of ERSP. The
%                 ERSP values are averaged over the selected frequency
%                 range. Not implemented
%   'limitbeg'  - [min max] limits at the beginning of the movie
%   'limitend'  - [min max] limits at the end of the movie
%   'freqslim'  - [freqs] array of landmark frequencies to set color
%                 limits.
%   'movieparams' - [low inc high] lower limit, higher limit and increment. If
%                increment is omited, movie is generate at every possible
%                increment (max freq resolution or max time resolution)
%
% Authors: Arnaud Delorme, CERCO, August, 2006
%
% See also: std_specplot(), std_erppplot(), std_erspplot()

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function [STUDY, M] = std_movie(STUDY, ALLEEG, varargin);

if nargin < 2
    help std_specplot;
    return;
end;

[opt addparams ] = finputcheck( varargin, ...
                             { 'erspfreq'    'real'    [] [];
                               'movieparams' 'real'    [] [];
                               'channels'    'cell'    [] {};
                               'freqslim'    'real'    [] [];
                               'limitbeg'    'real'    [] [];
                               'limitend'    'real'    [] [];
                               'moviemode'   'string'  { 'ersptime','erp','spec' } 'spec' }, 'std_movie', 'ignore');

if isstr(opt), error(opt); end;
tmpchanlocs =  ALLEEG(1).chanlocs;
if isempty(opt.channels), opt.channels = { tmpchanlocs.labels }; end;
if ~strcmpi(opt.moviemode, 'spec'), error('Only spec has been implemented so far'); end;

% read data once
% --------------
STUDY = std_specplot(STUDY, ALLEEG, 'channels', opt.channels, 'topofreq', [10 11]); close;

% find first data channel with info
% ---------------------------------
for cind = 1:length(STUDY.changrp)
    if ~isempty(STUDY.changrp(cind).specdata), break; end;
end;

% generate movie
% --------------
if isempty(opt.movieparams), 
    opt.movieparams(1) = STUDY.changrp(cind).specfreqs(1);
    opt.movieparams(2) = STUDY.changrp(cind).specfreqs(end);
end;
if length(opt.movieparams) == 2
    opt.movieparams(3) = opt.movieparams(2);
    opt.movieparams(2) = STUDY.changrp(cind).specfreqs(2)-STUDY.changrp(cind).specfreqs(1);
end;
if length(opt.movieparams) == 3
    opt.movieparams = [opt.movieparams(1):opt.movieparams(2):opt.movieparams(3)];
end;

% find limits
% -----------
if isempty(opt.limitbeg)
    [STUDY specdata] = std_specplot(STUDY, ALLEEG, 'channels', opt.channels, 'topofreq', opt.movieparams(1)); close;
    opt.limitbeg = [ min(specdata{1}) max(specdata{1}) ];
    [STUDY specdata] = std_specplot(STUDY, ALLEEG, 'channels', opt.channels, 'topofreq', opt.movieparams(end)); close;
    opt.limitend = [ min(specdata{1}) max(specdata{1}) ];
end;
lowlims  = linspace(opt.limitbeg(1), opt.limitend(1), length(opt.movieparams));
highlims = linspace(opt.limitbeg(2), opt.limitend(2), length(opt.movieparams));

% limits at specific frequencies
% ------------------------------
if ~isempty(opt.freqslim)
    if opt.freqslim(1)   ~= opt.movieparams(1)  , opt.freqslim = [ opt.movieparams(1) opt.freqslim ]; end;
    if opt.freqslim(end) ~= opt.movieparams(end), opt.freqslim = [ opt.freqslim opt.movieparams(end) ]; end;

    for ind = 1:length(opt.freqslim)
        [tmp indf(ind)] = min(abs(opt.freqslim(ind) - opt.movieparams));
        [STUDY specdata] = std_specplot(STUDY, ALLEEG, 'channels', opt.channels, 'topofreq', opt.movieparams(indf(ind))); close;
        minimum(ind) = min(specdata{1});
        maximum(ind) = max(specdata{1});
    end;
    indf(1)  = 0;
    lowlims  = [ ];
    highlims = [ ];
    for ind = 2:length(opt.freqslim)
        lowlims  = [lowlims linspace(minimum(ind-1), minimum(ind), indf(ind)-indf(ind-1)) ];
        highlims = [highlims linspace(maximum(ind-1), maximum(ind), indf(ind)-indf(ind-1)) ];
    end;
end;

% make movie
% ----------
for ind = 1:length(opt.movieparams)
    STUDY = std_specplot(STUDY, ALLEEG, 'channels', opt.channels, 'topofreq', opt.movieparams(ind), 'caxis', [lowlims(ind) highlims(ind)]);
    pos = get(gcf, 'position');
    set(gcf, 'position', [pos(1) pos(2) pos(3)*2 pos(4)*2]);
    M(ind) = getframe(gcf);
    close;
end;
figure; axis off; movie(M);
