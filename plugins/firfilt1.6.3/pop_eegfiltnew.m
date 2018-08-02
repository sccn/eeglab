% pop_eegfiltnew() - Filter data using Hamming windowed sinc FIR filter
%
% Usage:
%   >> [EEG, com, b] = pop_eegfiltnew(EEG); % pop-up window mode
%   >> [EEG, com, b] = pop_eegfiltnew(EEG, g.locutoff, g.hicutoff, g.filtorder,
%                                     g.revfilt, usefft, g.plotfreqz, minphase);
%
% Inputs:
%   EEG       - EEGLAB EEG structure
%   g.locutoff  - lower edge of the frequency pass band (Hz)
%               {[]/0 -> lowpass} 
%   g.hicutoff  - higher edge of the frequency pass band (Hz)
%               {[]/0 -> highpass}
%
% Optional inputs:
%   g.filtorder - filter order (filter length - 1). Mandatory even
%   g.revfilt   - [0|1] invert filter (from bandpass to notch filter)
%               {default 0 (bandpass)}
%   usefft    - ignored (backward compatibility only)
%   g.plotfreqz - [0|1] plot filter's frequency and phase response
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

function [EEG, com, b] = pop_eegfiltnew(EEG, varargin)

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

    geometry = {[3, 1], [3, 1], [3, 1], 1, 1, 1, 1 [2 1.5 0.5] [2 1.5 0.5]  };
    geomvert = [1 1 1 2 1 1 1 1 1];

    cb_type = 'pop_chansel(get(gcbf, ''userdata''), ''field'', ''type'',   ''handle'', findobj(''parent'', gcbf, ''tag'', ''chantypes''));';
    cb_chan = 'pop_chansel(get(gcbf, ''userdata''), ''field'', ''labels'', ''handle'', findobj(''parent'', gcbf, ''tag'', ''channels''));';
    
    uilist = {{'style', 'text', 'string', 'Lower edge of the frequency pass band (Hz)'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', 'Higher edge of the frequency pass band (Hz)'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', 'FIR Filter order (Mandatory even. Default is automatic*)'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', {'*See help text for a description of the default filter order heuristic.', 'Manual definition is recommended.'}} ...
              {'style', 'checkbox', 'string', 'Notch filter the data instead of pass band', 'value', 0} ...
              {'Style', 'checkbox', 'String', 'Use minimum-phase converted causal filter (non-linear!; beta)', 'Value', 0} ...
              {'style', 'checkbox', 'string', 'Plot frequency response', 'value', 1} ...
              { 'style' 'text'       'string' 'Channel type(s)' } ...
              { 'style' 'edit'       'string' '' 'tag' 'chantypes'}  ...
              { 'style' 'pushbutton' 'string' '...'  'callback' cb_type } ...
              { 'style' 'text'       'string' 'OR channel labels or indices' } ...
              { 'style' 'edit'       'string' '' 'tag' 'channels' }  ...
              { 'style' 'pushbutton' 'string' '...' 'callback' cb_chan }              
              };

    % channel labels
    % --------------
    if ~isempty(EEG(1).chanlocs)
        tmpchanlocs = EEG(1).chanlocs;        
    else
        tmpchanlocs = [];
        for index = 1:EEG(1).nbchan
            tmpchanlocs(index).labels = int2str(index);
            tmpchanlocs(index).type = '';
        end
    end
    
    result = inputgui('geometry', geometry, 'geomvert', geomvert, 'uilist', uilist, 'title', 'Filter the data -- pop_eegfiltnew()', 'helpcom', 'pophelp(''pop_eegfiltnew'')', 'userdata', tmpchanlocs);

    if isempty(result), return; end
    options = {};
    if ~isempty(result{1}), options = { options{:} 'locutoff' str2num( result{1}) }; end
    if ~isempty(result{2}), options = { options{:} 'hicutoff' str2num( result{2}) }; end
    if ~isempty(result{3}), options = { options{:} 'filtorder' result{3} }; end
    if result{4}, options = { options{:} 'revfilt' result{4} }; end
    if result{5}, options = { options{:} 'minphase' result{5} }; end
    if result{6}, options = { options{:} 'plotfreqz' result{6} }; end
    if ~isempty(result{7} ), options = { options{:} 'chantype' parsetxt(result{7}) }; end
    if ~isempty(result{8}) && isempty( result{7} )
       [ chaninds, chanlist ] = eeg_decodechan(EEG(1).chanlocs, result{8});
       if isempty(chanlist), chanlist = chaninds; end
       options = { options{:}, 'channels' chanlist };
    end        
elseif ~ischar(varargin{1})
    % backward compatibility
    options = {};
    if nargin > 1, options = { options{:} 'locutoff'  varargin{1} }; end
    if nargin > 2, options = { options{:} 'hicutoff'  varargin{2} }; end
    if nargin > 3, options = { options{:} 'filtorder' varargin{3} }; end
    if nargin > 4, options = { options{:} 'revfilt'   varargin{4} }; end
    if nargin > 5, options = { options{:} 'minphase'  varargin{5} }; end
    if nargin > 6, options = { options{:} 'plotfreqz' varargin{6} }; end
else
    options = varargin;
end

% decode inputs
% -------------
fieldlist = { 'locutoff'           'real'       []            []; 
              'hicutoff'           'real'       []            [];
              'filtorder'          'integer'    []            [];
              'revfilt'            'integer'    [0 1]         0;
              'usefft'             'integer'    [0 1]         0;
              'minphase'           'integer'    [0 1]         0;
              'plotfreqz'          'integer'    [0 1]         0;
              'channels'      {'cell' 'string' 'integer' } []                {};
              'chantype'      {'cell' 'string'} []                {}  };
g = finputcheck( options, fieldlist, 'pop_eegfiltnew');
if ischar(g), error(g); end
if ~isempty(g.chantype) 
    g.channels = eeg_decodechan(EEG.chanlocs, g.chantype, 'type');
elseif ~isempty(g.channels) 
    g.channels = eeg_decodechan(EEG.chanlocs, g.channels);
else
    g.channels = [1:EEG(1).nbchan];
end
if g.usefft
    error('FFT filtering not supported. Argument is provided for backward compatibility in command line mode only.')
end

% Constants
TRANSWIDTHRATIO = 0.25;
fNyquist = EEG.srate / 2;

% Check arguments
if g.locutoff == 0, g.locutoff = []; end
if g.hicutoff == 0, g.hicutoff = []; end
if isempty(g.hicutoff) % Convert highpass to inverted lowpass
    g.hicutoff = g.locutoff;
    g.locutoff = [];
    g.revfilt = ~g.revfilt;
end
edgeArray = sort([g.locutoff g.hicutoff]);

if isempty(edgeArray)
    error('Not enough input arguments.');
end
if any(edgeArray < 0 | edgeArray >= fNyquist)
    error('Cutoff frequency out of range');
end

if ~isempty(g.filtorder) && (g.filtorder < 2 || mod(g.filtorder, 2) ~= 0)
    error('Filter order must be a real, even, positive integer.')
end

% Max stop-band width
maxTBWArray = edgeArray; % Band-/highpass
if g.revfilt == 0 % Band-/lowpass
    maxTBWArray(end) = fNyquist - edgeArray(end);
elseif length(edgeArray) == 2 % Bandstop
    maxTBWArray = diff(edgeArray) / 2;
end
maxDf = min(maxTBWArray);

% Transition band width and filter order
if isempty(g.filtorder)

    % Default filter order heuristic
    if g.revfilt == 1 % Highpass and bandstop
        df = min([max([maxDf * TRANSWIDTHRATIO 2]) maxDf]);
    else % Lowpass and bandpass
        df = min([max([edgeArray(1) * TRANSWIDTHRATIO 2]) maxDf]);
    end

    g.filtorder = 3.3 / (df / EEG.srate); % Hamming window
    g.filtorder = ceil(g.filtorder / 2) * 2; % Filter order must be even.
    
else

    df = 3.3 / g.filtorder * EEG.srate; % Hamming window
    g.filtorderMin = ceil(3.3 ./ ((maxDf * 2) / EEG.srate) / 2) * 2;
    g.filtorderOpt = ceil(3.3 ./ (maxDf / EEG.srate) / 2) * 2;
    if g.filtorder < g.filtorderMin
        error('Filter order too low. Minimum required filter order is %d. For better results a minimum filter order of %d is recommended.', g.filtorderMin, g.filtorderOpt)
    elseif g.filtorder < g.filtorderOpt
        warning('firfilt:filterOrderLow', 'Transition band is wider than maximum stop-band width. For better results a minimum filter order of %d is recommended. Reported might deviate from effective -6dB cutoff frequency.', g.filtorderOpt)
    end

end

filterTypeArray = {'lowpass', 'bandpass'; 'highpass', 'bandstop (notch)'};
fprintf('pop_eegfiltnew() - performing %d point %s filtering.\n', g.filtorder + 1, filterTypeArray{g.revfilt + 1, length(edgeArray)})
fprintf('pop_eegfiltnew() - transition band width: %.4g Hz\n', df)
fprintf('pop_eegfiltnew() - passband edge(s): %s Hz\n', mat2str(edgeArray))

% Passband edge to cutoff (transition band center; -6 dB)
dfArray = {df, [-df, df]; -df, [df, -df]};
cutoffArray = edgeArray + dfArray{g.revfilt + 1, length(edgeArray)} / 2;
fprintf('pop_eegfiltnew() - cutoff frequency(ies) (-6 dB): %s Hz\n', mat2str(cutoffArray))

% Window
winArray = windows('hamming', g.filtorder + 1);

% Filter coefficients
if g.revfilt == 1
    filterTypeArray = {'high', 'stop'};
    b = firws(g.filtorder, cutoffArray / fNyquist, filterTypeArray{length(cutoffArray)}, winArray);
else
    b = firws(g.filtorder, cutoffArray / fNyquist, winArray);
end

if g.minphase
    disp('pop_eegfiltnew() - converting filter to minimum-phase (non-linear!)');
    b = minphaserceps(b);
end

% Plot frequency response
if g.plotfreqz
    try
        freqz(b, 1, 8192, EEG.srate);
    catch
        warning( 'Plotting of frequency response requires signal processing toolbox.' )
    end
end

% Filter
if g.minphase
    disp('pop_eegfiltnew() - filtering the data (causal)');
    EEG = firfiltsplit(EEG, b, 1, g.channels);
else
    disp('pop_eegfiltnew() - filtering the data (zero-phase)');
    EEG = firfilt(EEG, b, [], g.channels);
end


% History string
com = sprintf('EEG = pop_eegfiltnew(EEG, %s);', vararg2str(options));

end
