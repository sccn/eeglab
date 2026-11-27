% STD_FIXFFTFREQS - correct the frequency labels of spectrum files computed
%                   by STD_SPEC in 'fft' mode before EEGLAB issue 994 was
%                   fixed. These files labeled the kept FFT bins 2:floor(nfft/2)
%                   with linspace(0, srate/2, floor(nfft/2)) minus DC instead
%                   of k*srate/nfft. The spectral values are correct, so only
%                   the 'freqs' field is corrected and saved back to the file.
%
% Usage:
%   [ok, corrected] = std_fixfftfreqs(filename, EEG);
%
% Inputs:
%   filename  - [string or cell] spectrum file(s) (.datspec or .icaspec)
%   EEG       - EEG dataset(s) the file was computed from (only srate, pnts,
%               trials, xmin and times are used)
%
% Outputs:
%   ok        - false if a file computed in 'fft' mode has a frequency grid
%               that cannot be identified with certainty (it must then be
%               recomputed), true otherwise
%   corrected - true if at least one file was corrected

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2026
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

function [ok, corrected] = std_fixfftfreqs(filename, EEG)

ok        = true;
corrected = false;
if iscell(filename) % one file per session
    for iFile = 1:length(filename)
        [okFile, corrFile] = std_fixfftfreqs(filename{iFile}, EEG);
        ok        = ok && okFile;
        corrected = corrected || corrFile;
    end
    return;
end
if ~exist(filename, 'file'), return; end

warning('off', 'MATLAB:load:variableNotFound');
spec = load('-mat', filename, 'freqs', 'parameters');
warning('on', 'MATLAB:load:variableNotFound');
if ~isfield(spec, 'freqs') || isempty(spec.freqs) || ~isfield(spec, 'parameters') || ~strcmpi(getparam(spec.parameters, 'specmode', ''), 'fft')
    return;
end

f     = spec.freqs(:)';
n     = length(f);
srate = EEG(1).srate;
tol   = 1e-4*srate/(2*n); % small fraction of a frequency bin
if any([EEG.srate] ~= srate), ok = false; return; end

% correct grid: f = (1:n)*srate/nfft with n = floor(nfft/2)-1
nfft = round(srate/f(1));
if floor(nfft/2)-1 == n && max(abs(f - (1:n)*srate/nfft)) < tol, return; end

% legacy grid: f = (1:n)*srate/(2*n) with nfft = 2*n+2 or 2*n+3, which is
% obtained from the stored parameters and the datasets as in STD_SPEC
if max(abs(f - (1:n)*srate/(2*n))) >= tol, ok = false; return; end
nfft = getparam(spec.parameters, 'nfft', []);
if isempty(nfft)
    if all([EEG.trials] == 1) || strcmpi(getparam(spec.parameters, 'continuous', 'off'), 'on')
        epochlim = getparam(spec.parameters, 'epochlim', [0 1]);
        nfft = round(epochlim(2)*srate-1) - round(epochlim(1)*srate) + 1; % length of EEG_REGEPOCHS epochs (see EPOCH)
    elseif any([EEG.pnts] ~= EEG(1).pnts) || any([EEG.xmin] ~= EEG(1).xmin)
        ok = false; return;
    else
        timerange = getparam(spec.parameters, 'timerange', []);
        if isempty(timerange)
            nfft = EEG(1).pnts;
        else
            nfft = sum(EEG(1).times >= timerange(1) & EEG(1).times < timerange(2));
        end
    end
end
if floor(nfft/2)-1 ~= n, ok = false; return; end

spec = load('-mat', filename);
spec.freqs = (1:n)*srate/nfft;
std_savedat(filename, spec);
corrected = true;
fprintf('Corrected the frequency labels of FFT spectrum file %s (EEGLAB issue 994)\n', filename);

% last value of a parameter in a { 'key' val ... } cell array
% ------------------------------------------------------------
function val = getparam(params, key, default)
val = default;
ind = find(strcmpi(params(1:2:end), key), 1, 'last');
if ~isempty(ind), val = params{2*ind}; end
