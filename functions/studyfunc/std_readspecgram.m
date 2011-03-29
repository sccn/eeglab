% std_readspecgram() - returns the stored mean power spectrogram for an ICA component 
%                  or a data channel in a specified dataset.  The spectrogram is 
%                  assumed to have been saved in a Matlab file, 
%                  "[dataset_name].datspecgram", in the same
%                  directory as the dataset file. If this file doesn't exist,
%                  use std_specgram() to create it.
% Usage:    
%  >> [spec, times, freqs] = std_readspecgram(ALLEEG, setindx, component, timerange, freqrange);  
%
% Inputs:
%   ALLEEG     - a vector of dataset EEG structures (may also be one dataset). 
%                Must contain the dataset of interest (the 'setindx' below).
%   setindx    - [integer] an index of an EEG dataset in the ALLEEG
%                structure for which to read a component spectrum.
%   component  - [integer] index of the component in the selected EEG dataset 
%                for which to return the spectrum
%   freqrange  - [min max in Hz] frequency range to return
%
%
% Outputs:
%   spec      - the log-power spectrum of the requested ICA component in the
%               specified dataset (in dB)
%   freqs     - vector of spectral frequencies (in Hz)
%
%  See also  std_spec(), pop_preclust(), std_preclust()
%
% Authors:  Arnaud Delorme, SCCN, INC, UCSD, February, 2008

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, October 11, 2004, hilit@sccn.ucsd.edu
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

function [X, t, f] = std_readspecgram(ALLEEG, abset, comp, timerange, freqrange, rmsubjmean);

if nargin < 4
    timerange = [];
end;
if nargin < 5
    freqrange = [];
end;
    
X = [];

if iscell(comp)
    % find channel indices list
    % -------------------------
    chanind  = [];
    tmpchanlocs = ALLEEG(abset).chanlocs;
    chanlabs = lower({ tmpchanlocs.labels });
    for index = 1:length(comp)
        tmp = strmatch(lower(comp{index}), chanlabs, 'exact');
        if isempty(tmp)
            error([ 'Channel ''' comp{index} ''' not found in dataset ' int2str(abset)]);
        else    
            chanind = [ chanind tmp ];
        end;
    end;
    filename = fullfile( ALLEEG(abset).filepath,[ ALLEEG(abset).filename(1:end-3) 'datspecgram']);
    prefix = 'chan';
    inds   = chanind;
elseif comp(1) < 0
    filename = fullfile( ALLEEG(abset).filepath,[ ALLEEG(abset).filename(1:end-3) 'datspecgram']);
    prefix = 'chan';
    inds   = -comp;
else
    filename = fullfile( ALLEEG(abset).filepath,[ ALLEEG(abset).filename(1:end-3) 'icaspecgram']);
    prefix = 'comp';
    inds   = comp;
end;

for k=1:length(inds)
    try,
        warning backtrace off;
        erpstruct = load( '-mat', filename, [ prefix int2str(inds(k)) ], 'freqs', 'times');
        warning backtrace on;
    catch
        error( [ 'Cannot read file ''' filename '''' ]);
    end;

    tmpdat    = getfield(erpstruct, [ prefix int2str(inds(k)) ]);
    if k == 1
        X = zeros([size(tmpdat) length(comp)]);
    end;
    X(:,:,k) = tmpdat;
    f        = getfield(erpstruct, 'freqs');
    t        = getfield(erpstruct, 'times');
end;

% select frequency range of interest
% ----------------------------------
if ~isempty(freqrange)
    maxind = max(find(f <= freqrange(end)));
    minind = min(find(f >= freqrange(1)));
else
    %if not, use whole spectrum
    maxind = length(f);
    minind = 1;
end
f = f(minind:maxind);
X = X(minind:maxind,:,:);

% select time range of interest
% -----------------------------
if ~isempty(timerange)
    maxind = max(find(t <= timerange(end)));
    minind = min(find(t >= timerange(1)));
else
    %if not, use whole spectrum
    maxind = length(t);
    minind = 1;
end
t = t(minind:maxind);
X = X(:,minind:maxind,:);

return;
