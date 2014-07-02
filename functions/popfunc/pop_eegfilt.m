% pop_eegfilt() - interactively filter EEG dataset data using eegfilt()
%
% Usage:
%   >> EEGOUT = pop_eegfilt( EEG, locutoff, hicutoff, filtorder);
%
% Graphical interface:
%   "Lower edge ..." - [edit box] Lower edge of the frequency pass band (Hz) 
%                 Same as the 'locutoff' command line input.
%   "Higher edge ..." - [edit box] Higher edge of the frequency pass band (Hz) 
%                 Same as the 'hicutoff' command line input.
%   "Notch filter" - [edit box] provide the notch range, i.e. [45 55] 
%                 for 50 Hz). This option overwrites the low and high edge limits
%                 given above. Set the 'locutoff' and 'hicutoff' values to the
%                 values entered as parameters, and set 'revfilt to 1, to swap
%                 from bandpass to notch filtering.
%   "Filter length" - [edit box] Filter lenghth in point (default: see 
%                 >> help eegfilt). Same as 'filtorder' optional input.
%
% Inputs:
%   EEG       - input dataset
%   locutoff  - lower edge of the frequency pass band (Hz)  {0 -> lowpass}
%   hicutoff  - higher edge of the frequency pass band (Hz) {0 -> highpass}
%   filtorder - length of the filter in points {default 3*fix(srate/locutoff)}
%   revfilt   - [0|1] Reverse filter polarity (from bandpass to notch filter). 
%                     Default is 0 (bandpass).
%   usefft    - [0|1] 1 uses FFT filtering instead of FIR. Default is 0.
%   plotfreqz - [0|1] plot frequency response of filter. Default is 0.
%   firtype   - ['firls'|'fir1'] filter design method, default is 'firls'
%               from the command line
%   causal    - [0|1] 1 uses causal filtering. Default is 0.
%
% Outputs:
%   EEGOUT   - output dataset
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eegfilt(), eegfiltfft(), eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% 01-25-02 reformated help & license -ad 

function [EEG, com] = pop_eegfilt( EEG, locutoff, hicutoff, filtorder, revfilt, usefft, plotfreqz, firtype, causal)

com = '';
if nargin < 1
    help pop_eegfilt;
    return;
end;
if isempty(EEG(1).data)
    disp('Pop_eegfilt() error: cannot filter an empty dataset'); return;
end;

% warning
% -------
if exist('filtfilt') ~= 2
    disp('Warning: cannot find the signal processing toolbox');
    disp('         a simple fft/inverse fft filter will be used');
    usefft = 1;
end;

if nargin < 2
    % which set to save
    % -----------------
    uilist = { ...
        { 'style' 'text' 'string' 'Lower edge of the frequency pass band (Hz)' } ...
        { 'style' 'edit' 'string' '' } ...
        { 'style' 'text' 'string' 'Higher edge of the frequency pass band (Hz)' } ...
        { 'style' 'edit' 'string' '' } ...
        { 'style' 'text' 'string' 'FIR Filter order (default is automatic)' } ...
        { 'style' 'edit' 'string' '' } ...
        { 'style' 'checkbox' 'string' 'Notch filter the data instead of pass band' } ...
        { 'style' 'checkbox' 'string' 'Use (sharper) FFT linear filter instead of FIR filtering' 'value' 0 } ...
        { 'style' 'text' 'string' '(Use the option above if you do not have the Signal Processing Toolbox)' } ...
        { 'style' 'checkbox' 'string' 'Use causal filter (useful when performing causal analysis)' 'value' 0} ...
        { 'style' 'checkbox' 'string' 'Plot the filter frequency response' 'value' 0} ...
        { 'style' 'checkbox' 'string' 'Use fir1 (check, recommended) or firls (uncheck, legacy)' 'value' 1}};
    geometry = { [3 1] [3 1] [3 1] 1 1 1 1 1 1 };
    
    result = inputgui( 'geometry', geometry, 'uilist', uilist, 'title', 'Filter the data -- pop_eegfilt()', ...
        'helpcom', 'pophelp(''pop_eegfilt'')');
    
    if isempty(result), return; end;
    if isempty(result{1}), result{1} = '0'; end;
    if isempty(result{2}), result{2} = '0'; end;
    
    locutoff   	 = eval( result{1} );
    hicutoff 	 = eval( result{2} );
    if isempty( result{3} )
        filtorder = [];
    else filtorder    = eval( result{3} );
    end;
    revfilt = 0;
    if result{4},
        revfilt = 1;
        if locutoff == 0 | hicutoff == 0,
            error('Need both lower and higher edge for notch filter');
        end;
    end;
    if result{5}
         usefft = 1;
    else usefft = 0;
    end;
    if result{6}
         causal = 1;
    else causal = 0;
    end;
    plotfreqz = result{7};
    if locutoff == 0 & hicutoff == 0 return; end;
    if result{8}
        firtype = 'fir1';
    else
        firtype = 'firls';
    end
else
    if nargin < 3
        hicutoff = 0;
    end;
    if nargin < 4
        filtorder = [];
    end;
    if nargin < 5
        revfilt = 0;
    end;
    if nargin < 6
        usefft = 0;
    end;
    if nargin < 7
        plotfreqz = 0;
    end
    if nargin < 8
        firtype = 'firls';
    end
    if nargin < 8
        causal = 0;
    end
end;

if locutoff && hicutoff
    disp('WARNING: BANDPASS FILTERS SOMETIMES DO NOT WORK (MATLAB BUG)')
    disp('WARNING: PLOT SPECTRUM AFTER FILTERING TO ASSESS FILTER EFFICIENCY')
    disp('WARNING: IF FILTER FAILS, LOWPASS DATA THEN HIGHPASS DATA')
end;

% process multiple datasets
% -------------------------
if length(EEG) > 1
    [ EEG com ] = eeg_eval( 'pop_eegfilt', EEG, 'warning', 'on', 'params', ...
        { locutoff, hicutoff, filtorder, revfilt } );
    return;
end;

options = { EEG.srate, locutoff, hicutoff, 0 };
if ~isempty( filtorder )
    options = { options{:} filtorder };
else
    options = { options{:} 0 };
end;

options = {options{:} revfilt firtype causal};

if EEG.trials == 1
    if ~isempty(EEG.event) & isfield(EEG.event, 'type') & isstr(EEG.event(1).type)
        tmpevent = EEG.event;
        boundaries = strmatch('boundary', { tmpevent.type });
        if isempty(boundaries)
            if ~usefft
                [EEG.data, b] = eegfilt( EEG.data, options{:});
            else
                EEG.data = eegfiltfft( EEG.data, options{1:6});            % 7/30/2014 Ramon: {:} to {1:6}; 
            end;
        else
            options{4} = 0;
            disp('Pop_eegfilt:finding continuous data boundaries');
            tmplat = [ tmpevent.latency ];
            boundaries = tmplat(boundaries);
            boundaries = [0 floor(boundaries-0.49) EEG.pnts];
            try, warning off MATLAB:divideByZero
            catch, end;
            for n=1:length(boundaries)-1
                if boundaries(n)+1 < boundaries(n+1)
                    try
                        fprintf('Processing continuous data (%d:%d)\n',boundaries(n),boundaries(n+1));
                        if ~usefft
                            [EEG.data(:,boundaries(n)+1:boundaries(n+1)), b] = ...
                                eegfilt(EEG.data(:,boundaries(n)+1:boundaries(n+1)), options{:});
                        else
                            EEG.data(:,boundaries(n)+1:boundaries(n+1)) = ...
                                eegfiltfft(EEG.data(:,boundaries(n)+1:boundaries(n+1)), options{1:6});  % 7/30/2014 Ramon: {:} to {1:6}; 
                        end;
                    catch
                        fprintf('\nFilter error: continuous data portion too narrow (DC removed if highpass only)\n');
                        if locutoff ~= 0 & hicutoff == 0
                            tmprange = [boundaries(n)+1:boundaries(n+1)];
                            EEG.data(:,tmprange) = ...
                                EEG.data(:,tmprange) - repmat(mean(EEG.data(:,tmprange),2), [1 length(tmprange)]);
                        end;
                    end;
                end;
            end
            try, warning on MATLAB:divideByZero
            catch, end;
        end
    else
        if ~usefft
            [EEG.data, b] = eegfilt( EEG.data, options{:});
        else
            EEG.data = eegfiltfft( EEG.data, options{1:6});                % 7/30/2014 Ramon: {:} to {1:6}; 
        end;
    end;
else
    EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials);
    options{4} = EEG.pnts;
    if ~usefft
        [EEG.data, b] = eegfilt( EEG.data, options{:});
    else
        EEG.data = eegfiltfft( EEG.data, options{1:6});                    % 7/30/2014 Ramon: {:} to {1:6}; 
    end;
    % Note: reshape does not reserve new memory while EEG.data(:,:) does
end;
EEG.icaact = [];

if ~usefft & plotfreqz & exist('b') == 1
    freqz(b, 1, [], EEG.srate);
end

com = sprintf( '%s = pop_eegfilt( %s, %s, %s, [%s], [%s], %s, %s, ''%s'', %d);', inputname(1), inputname(1), ...
    num2str( locutoff), num2str( hicutoff), num2str( filtorder ), num2str( revfilt ), num2str(usefft), num2str(plotfreqz), firtype, causal);
return
