% eeg_rereject() - reject the same data portion as the ones defined
%                  in a reference dataset. This function is useful when
%                  have filtered a dataset and want to reject again the
%                  same portions of data.
% Usage:
%   >> OUTEEG = eeg_rereject(EEG, cleanEEG)
%
% Inputs:
%   EEG       - input dataset
%   cleanEEG  - EEGLAB dataset with portion of data rejected (information 
%               stored in boundary event types).
%
% Optional inputs:
%   'blank'    - ['on'|'off'] replace by 0 instead of removing the data
%                portions. Default is 'off'. Note that this does not
%                add any boundary events.
%   'rmdiscont' - ['on'|'off'] remove discontinuities (by shifting the
%                data in the different discontinuous portions). Note that
%                this discards the boundary events.
%
% Outputs:
%   OUTEEG     - output dataset
%
% Author: Arnaud Delorme, 2020
%
% Note: uses the resample() function from the signal processing toolbox
%       if present. Otherwise use griddata interpolation method (it should be
%       reprogrammed using spline interpolation for speed up).
%
% See also: resample(), eeglab()

% Copyright (C) 2020 Arnaud Delorme
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

function [EEG, boundloc]= eeg_rereject(EEG, CLEANEDEEG, varargin)

boundLat = [];
if nargin < 2
    help eeg_rereject;
    return;
end

% check inputs
% ------------
g = finputcheck(varargin, { 'blank'      'string'  {'on','off' }             'off';
                            'rmdiscont'  'string'  {'on','off' }             'off' });

if isstr(g), error(g); end

boundEvents = eeg_findboundaries(CLEANEDEEG);
boundloc = [ CLEANEDEEG.event(boundEvents).latency ];
dur      = [ CLEANEDEEG.event(boundEvents).duration ];
cumdur   = cumsum(dur);
boundloc = boundloc + [0 cumdur(1:end-1) ];
boundloc = [ boundloc; boundloc+dur-1]';
if strcmpi(g.blank, 'off')
    EEG = eeg_eegrej(EEG, ceil(boundloc));
    boundEvents = eeg_findboundaries(EEG.event);
    boundLat    = [ EEG.event(boundEvents).latency ];
else
    boundloc = round(boundloc);
    for iBound = size(boundloc,1):-1:1
        EEG.data(:,boundloc(iBound,1):boundloc(iBound,2)) = 0;
    end
    boundEvents = [];
    boundlocTmp = boundloc;
    boundlocTmp(:,1) = boundlocTmp(:,1)-1; % to compensate round(lat+0.5) below
    boundLat = sort(boundlocTmp(:))';
end

% remake data continuous
if strcmpi(g.rmdiscont, 'on')
    for ind = length(boundLat):-1:1
        if boundLat(ind) > 1
            lat = min(round(boundLat(ind)+0.5), size(EEG.data,2));
            EEG.data(:,lat:end) = EEG.data(:,lat:end) - EEG.data(:,lat) + EEG.data(:,lat-1);
            EEG.data(:,lat) = []; % remove 1 sample
            if ~isempty(boundEvents)
                for ind2 = boundEvents(ind):length(EEG.event)
                    EEG.event(ind2).latency = EEG.event(ind2).latency-1; % shift events by 1 sample
                end
            end
        end
    end
    EEG.event(boundEvents) = [];
    EEG.pnts = size(EEG.data,2);
end
EEG = eeg_checkset(EEG);
