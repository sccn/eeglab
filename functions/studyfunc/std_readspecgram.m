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

function [X, t, f] = std_readspecgram(ALLEEG, abset, comp, timerange, freqrange, rmsubjmean);

if nargin < 4
    timerange = [];
end
if nargin < 5
    freqrange = [];
end
    
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
        end
    end
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
end

for k=1:length(inds)
    try,
        warning backtrace off;
        erpstruct = load( '-mat', filename, [ prefix int2str(inds(k)) ], 'freqs', 'times');
        warning backtrace on;
    catch
        error( [ 'Cannot read file ''' filename '''' ]);
    end

    tmpdat    = getfield(erpstruct, [ prefix int2str(inds(k)) ]);
    if k == 1
        X = zeros([size(tmpdat) length(comp)]);
    end
    X(:,:,k) = tmpdat;
    f        = getfield(erpstruct, 'freqs');
    t        = getfield(erpstruct, 'times');
end

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
