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

function EEG = eeg_timeinterp( EEG, samples, varargin);

    if nargin < 2
        help eeg_timeinterp;
        return;
    end;
    
    opt = finputcheck(varargin, { 'epochinds'   'integer'   []       [1:EEG.trials]; ...
                                  'interpwin'   'integer'   []       5; ...
                                  'elecinds'    'integer'   []       [1:EEG.nbchan]; ...
                                  'epochcont'   'string'    { 'on';'off' }  'off' }, 'eeg_timeinterp');

    if isstr(opt), error(opt); end;
    
    srange = samples(2)-samples(1);
    data   = EEG.data;
    pnts   = EEG.pnts;
    
    if strcmpi(opt.epochcont, 'on')
        data(:,end+1:end+srange*opt.interpwin,1:end-1) = data(:,1:srange*opt.interpwin,2:end);
        pnts = pnts + srange*opt.interpwin;
    end;
    
    % determine region to interpolate
    % and region to use for interpolation
    % -----------------------------------
    samplesin  = [min(samples(1)-srange*opt.interpwin,1):samples(1)-1 samples(2)+1:min(samples(2)+srange*opt.interpwin, pnts)];
    samplesout = [samples(1):samples(2)];
    
    if length(opt.epochinds) > 1, fprintf('Trials:'); end;
    for index = opt.epochinds
        for elec = opt.elecinds
            EEG.data(elec,samplesout,index) = spline( samplesin, data(elec, samplesin, index), samplesout);
        end;
        if length(opt.epochinds) > 1, 
            fprintf('.'); 
            if mod(index,40) == 01, fprintf('\n'); end;
        end;
    end;
    if length(opt.epochinds) > 1, fprintf('\n'); end;
