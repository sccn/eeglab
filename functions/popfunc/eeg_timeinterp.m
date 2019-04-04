% eeg_timeinterp() - perform spline interpolation of a portion
%                    of data based on prior and post activity. See
%                    eeg_interp() for interpolation of bad channels.
%
% Usage:
%   >> OUTEEG = eeg_timeinterp( INEEG, samples, 'key', 'val');
%
% Inputs:
%   INEEG       - input EEG structure
%   samplerange - [min max] range sample points in continuous 
%                 or epoched data. Only one sample point range 
%                 may be given at a time.
%
% Optional inputs:
%   'elecinds'  - indices of electrodes to interpolate (default all)
%   'epochinds' - indices of epochs for epoched data (default all)
%   'interpwin' - [integer] number of data points to use before and
%                 after the sample range. Default is 5 times the 
%                 interpolated sample range.
%   'epochcont' - ['on'|'off'] epochs are contiguous. Only works for
%                 interpolating the end of epochs (default 'off')
%
% Outputs:
%   OUTEEG         - EEGLAB data structure
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, 2007

% Copyright (C) Arnaud Delorme, SCCN/INC/UCSD, 2007
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

function EEG = eeg_timeinterp( EEG, samples, varargin);

    if nargin < 2
        help eeg_timeinterp;
        return;
    end
    
    opt = finputcheck(varargin, { 'epochinds'   'integer'   []       [1:EEG.trials]; ...
                                  'interpwin'   'integer'   []       5; ...
                                  'elecinds'    'integer'   []       [1:EEG.nbchan]; ...
                                  'epochcont'   'string'    { 'on';'off' }  'off' }, 'eeg_timeinterp');

    if ischar(opt), error(opt); end
    
    srange = samples(2)-samples(1);
    data   = EEG.data;
    pnts   = EEG.pnts;
    
    if strcmpi(opt.epochcont, 'on')
        data(:,end+1:end+srange*opt.interpwin,1:end-1) = data(:,1:srange*opt.interpwin,2:end);
        pnts = pnts + srange*opt.interpwin;
    end
    
    % determine region to interpolate
    % and region to use for interpolation
    % -----------------------------------
    samplesin  = [min(samples(1)-srange*opt.interpwin,1):samples(1)-1 samples(2)+1:min(samples(2)+srange*opt.interpwin, pnts)];
    samplesout = [samples(1):samples(2)];
    
    if length(opt.epochinds) > 1, fprintf('Trials:'); end
    for index = opt.epochinds
        for elec = opt.elecinds
            EEG.data(elec,samplesout,index) = spline( samplesin, data(elec, samplesin, index), samplesout);
        end
        if length(opt.epochinds) > 1, 
            fprintf('.'); 
            if mod(index,40) == 01, fprintf('\n'); end
        end
    end
    if length(opt.epochinds) > 1, fprintf('\n'); end
