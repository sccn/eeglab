% pop_eegfiltnew() - Filter data using Hamming windowed sinc FIR filter
%
% Usage:
%   >> [EEG, com, b] = pop_eegfiltnew(EEG); % pop-up window mode
%   >> [EEG, com, b] = pop_eegfiltnew(EEG, locutoff, hicutoff, filtorder,
%                                     revfilt, usefft, plotfreqz, minphase);
%
% Inputs:
%   EEG       - EEGLAB EEG structure
%   locutoff  - lower edge of the frequency pass band (Hz)
%               {[]/0 -> lowpass} 
%   hicutoff  - higher edge of the frequency pass band (Hz)
%               {[]/0 -> highpass}
%
% Optional inputs:
%   filtorder - filter order (filter length - 1). Mandatory even
%   revfilt   - [0|1] invert filter (from bandpass to notch filter)
%               {default 0 (bandpass)}
%   usefft    - ignored (backward compatibility only)
%   plotfreqz - [0|1] plot filter's frequency and phase response
%               {default 0} 
%   minphase  - scalar boolean minimum-phase converted causal filter
%               {default false}
%
% Outputs:
%   EEG       - filtered EEGLAB EEG structure
%   com       - history string
%   b         - filter coefficients
%
% Note:
%   pop_eegfiltnew is intended as a replacement for the deprecated
%   pop_eegfilt function. Required filter order/transition band width is
%   estimated with the following heuristic in default mode: transition band
%   width is 25% of the lower passband edge, but not lower than 2 Hz, where
%   possible (for bandpass, highpass, and bandstop) and distance from
%   passband edge to critical frequency (DC, Nyquist) otherwise. Window
%   type is hardcoded to Hamming. Migration to windowed sinc FIR filters
%   (pop_firws) is recommended. pop_firws allows user defined window type
%   and estimation of filter order by user defined transition band width.
%
% Author: Andreas Widmann, University of Leipzig, 2012
%
% See also:
%   firfilt, firws, windows

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2008 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

function [EEG, com, b] = pop_eegfiltnew(EEG, locutoff, hicutoff, filtorder, revfilt, usefft, plotfreqz, minphase)

com = '';

if nargin < 1
    help pop_eegfiltnew;
    return
end
if isempty(EEG.data)
    error('Cannot filter empty dataset.');
end

% GUI
if nargin < 2

    geometry = {[3, 1], [3, 1], [3, 1], 1, 1, 1, 1};
    geomvert = [1 1 1 2 1 1 1];

    uilist = {{'style', 'text', 'string', 'Lower edge of the frequency pass band (Hz)'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', 'Higher edge of the frequency pass band (Hz)'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', 'FIR Filter order (Mandatory even. Default is automatic*)'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', {'*See help text for a description of the default filter order heuristic.', 'Manual definition is recommended.'}} ...
              {'style', 'checkbox', 'string', 'Notch filter the data instead of pass band', 'value', 0} ...
              {'Style', 'checkbox', 'String', 'Use minimum-phase converted causal filter (non-linear!; beta)', 'Value', 0} ...
              {'style', 'checkbox', 'string', 'Plot frequency response', 'value', 1}};

    result = inputgui('geometry', geometry, 'geomvert', geomvert, 'uilist', uilist, 'title', 'Filter the data -- pop_eegfiltnew()', 'helpcom', 'pophelp(''pop_eegfiltnew'')');

    if isempty(result), return; end

    locutoff = str2num(result{1});
    hicutoff = str2num(result{2});
    filtorder = str2num(result{3});
    revfilt = result{4};
    minphase = result{5};
    plotfreqz = result{6};
    usefft = [];

else
    
    if nargin < 3
        hicutoff = [];
    end
    if nargin < 4
        filtorder = [];
    end
    if nargin < 5 || isempty(revfilt)
        revfilt = 0;
    end
    if nargin < 6
        usefft = [];
    elseif usefft == 1
        error('FFT filtering not supported. Argument is provided for backward compatibility in command line mode only.')
    end
    if nargin < 7 || isempty(plotfreqz)
        plotfreqz = 0;
    end
    if nargin < 8 || isempty(minphase)
        minphase = 0;
    end
    
end

% Constants
TRANSWIDTHRATIO = 0.25;
fNyquist = EEG.srate / 2;

% Check arguments
if locutoff == 0, locutoff = []; end
if hicutoff == 0, hicutoff = []; end
if isempty(hicutoff) % Convert highpass to inverted lowpass
    hicutoff = locutoff;
    locutoff = [];
    revfilt = ~revfilt;
end
edgeArray = sort([locutoff hicutoff]);

if isempty(edgeArray)
    error('Not enough input arguments.');
end
if any(edgeArray < 0 | edgeArray >= fNyquist)
    error('Cutoff frequency out of range');
end

if ~isempty(filtorder) && (filtorder < 2 || mod(filtorder, 2) ~= 0)
    error('Filter order must be a real, even, positive integer.')
end

% Max stop-band width
maxTBWArray = edgeArray; % Band-/highpass
if revfilt == 0 % Band-/lowpass
    maxTBWArray(end) = fNyquist - edgeArray(end);
elseif length(edgeArray) == 2 % Bandstop
    maxTBWArray = diff(edgeArray) / 2;
end
maxDf = min(maxTBWArray);

% Transition band width and filter order
if isempty(filtorder)

    % Default filter order heuristic
    if revfilt == 1 % Highpass and bandstop
        df = min([max([maxDf * TRANSWIDTHRATIO 2]) maxDf]);
    else % Lowpass and bandpass
        df = min([max([edgeArray(1) * TRANSWIDTHRATIO 2]) maxDf]);
    end

    filtorder = 3.3 / (df / EEG.srate); % Hamming window
    filtorder = ceil(filtorder / 2) * 2; % Filter order must be even.
    
else

    df = 3.3 / filtorder * EEG.srate; % Hamming window
    filtorderMin = ceil(3.3 ./ ((maxDf * 2) / EEG.srate) / 2) * 2;
    filtorderOpt = ceil(3.3 ./ (maxDf / EEG.srate) / 2) * 2;
    if filtorder < filtorderMin
        error('Filter order too low. Minimum required filter order is %d. For better results a minimum filter order of %d is recommended.', filtorderMin, filtorderOpt)
    elseif filtorder < filtorderOpt
        warning('firfilt:filterOrderLow', 'Transition band is wider than maximum stop-band width. For better results a minimum filter order of %d is recommended. Reported might deviate from effective -6dB cutoff frequency.', filtorderOpt)
    end

end

filterTypeArray = {'lowpass', 'bandpass'; 'highpass', 'bandstop (notch)'};
fprintf('pop_eegfiltnew() - performing %d point %s filtering.\n', filtorder + 1, filterTypeArray{revfilt + 1, length(edgeArray)})
fprintf('pop_eegfiltnew() - transition band width: %.4g Hz\n', df)
fprintf('pop_eegfiltnew() - passband edge(s): %s Hz\n', mat2str(edgeArray))

% Passband edge to cutoff (transition band center; -6 dB)
dfArray = {df, [-df, df]; -df, [df, -df]};
cutoffArray = edgeArray + dfArray{revfilt + 1, length(edgeArray)} / 2;
fprintf('pop_eegfiltnew() - cutoff frequency(ies) (-6 dB): %s Hz\n', mat2str(cutoffArray))

% Window
winArray = windows('hamming', filtorder + 1);

% Filter coefficients
if revfilt == 1
    filterTypeArray = {'high', 'stop'};
    b = firws(filtorder, cutoffArray / fNyquist, filterTypeArray{length(cutoffArray)}, winArray);
else
    b = firws(filtorder, cutoffArray / fNyquist, winArray);
end

if minphase
    disp('pop_eegfiltnew() - converting filter to minimum-phase (non-linear!)');
    b = minphaserceps(b);
end

% Plot frequency response
if plotfreqz
    freqz(b, 1, 8192, EEG.srate);
end

% Filter
if minphase
    disp('pop_eegfiltnew() - filtering the data (causal)');
    EEG = firfiltsplit(EEG, b, 1);
else
    disp('pop_eegfiltnew() - filtering the data (zero-phase)');
    EEG = firfilt(EEG, b);
end


% History string
com = sprintf('%s = pop_eegfiltnew(%s, %s);', inputname(1), inputname(1), vararg2str({locutoff, hicutoff, filtorder, revfilt, usefft, plotfreqz}));

end
