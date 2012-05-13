% pop_iirfilt() - interactively filter EEG dataset data using iirfilt()
%
% Usage:
%   >> EEGOUT = pop_iirfilt( EEG, locutoff, hicutoff);
%   >> EEGOUT = pop_iirfilt( EEG, locutoff, hicutoff, trans_bw, revfilt, causal);
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
%                 >> help pop_iirfilt). Same as 'trans_bw' optional input.
%
% Inputs:
%   EEG       - input dataset
%   locutoff  - lower edge of the frequency pass band (Hz)  {0 -> lowpass}
%   hicutoff  - higher edge of the frequency pass band (Hz) {0 -> highpass}
%   trans_bw  - length of the filter in points {default 3*fix(srate/locutoff)}
%   revfilt   - [0|1] Reverse filter polarity (from bandpass to notch filter).
%                     Default is 0 (bandpass).
%   causal    - [0|1] use causal filter.
%
% Outputs:
%   EEGOUT   - output dataset
%
% Authors: Maksym Pozdin (mpozdin.ece04@gtalumni.org, IOL/ONRC,2004), 
%          with Arnaud Delorme and Scott Makeig (SCCN/INC/UCSD, La Jolla CA)
%
% See also: iirfilt(), pop_eegfilt(), eegfilt(), eegfiltfft(), eeglab()

% Copyright (C) 2004 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [EEG, com] = pop_iirfilt( EEG, locutoff, hicutoff, trans_bw, revfilt, causal);

com = '';
if nargin < 1
	help pop_iirfilt;
	return;
end;	
if isempty(EEG(1).data)
    disp('pop_iirfilt() error: cannot filter an empty dataset'); return;
end;    
if nargin < 2
    % which set to save
    % -----------------
    uilist = { ...
        { 'style' 'text' 'string' 'Lower edge of the frequency pass band (Hz)' } ...
        { 'style' 'edit' 'string' '' } ...
        { 'style' 'text' 'string' 'Higher edge of the frequency pass band (Hz)' } ...
        { 'style' 'edit' 'string' '' } ...
        { 'style' 'text' 'string' 'Transition BW filter (default is automatic)' } ...
        { 'style' 'edit' 'string' '' } ...
        { 'style' 'checkbox' 'string' 'Notch filter the data instead of pass band' } ...
        { 'style' 'checkbox' 'string' 'Use causal filter (useful when performing causal analysis)' 'value' 0} ...
        { 'style' 'checkbox' 'string' 'Plot the filter frequency response' 'value' 0} };
    geometry = { [3 1] [3 1] [3 1] 1 1 1 };
    
    result = inputgui( 'geometry', geometry, 'uilist', uilist, 'title', 'Filter the data -- pop_iirfilt()', ...
        'helpcom', 'pophelp(''pop_iirfilt'')');
    
    if isempty(result), return; end;
    if isempty(result{1}), result{1} = '0'; end;
    if isempty(result{2}), result{2} = '0'; end;
    
    locutoff   	 = eval( result{1} );
    hicutoff 	 = eval( result{2} );
    if isempty( result{3} )
         trans_bw = [];
    else trans_bw = eval( result{3} );
    end;
    revfilt = 0;
    if result{4},
        revfilt = 1;
    end;
    if result{5}
         causal = 1;
    else causal = 0;
    end;
    plotfreqz = result{6};
else
    if nargin < 3
        hicutoff = 0;
    end;
    if nargin < 4
        trans_bw = [];
    end;
    if nargin < 5
        revfilt = 0;
    end;
    if nargin < 6
        causal = 0;
    end;
    plotfreqz = 0;
end;
 
options = { EEG.srate, locutoff, hicutoff, EEG.pnts };
if ~isempty( trans_bw )
	options = { options{:} trans_bw };
else 
	options = { options{:} 0 };
end;
if revfilt ~= 0
	options = { options{:} revfilt [] [] causal };
else
	options = { options{:} 0       [] [] causal };
end;

% warning
% -------
if exist('filtfilt') ~= 2
    disp('Warning: cannot find the signal processing toolbox');
    disp('         a simple fft/inverse fft filter will be used');
end;

% process multiple datasets
% -------------------------
if length(EEG) > 1
   [ EEG com ] = eeg_eval( 'pop_iirfilt', EEG, 'warning', 'on', 'params', ...
                           { locutoff, hicutoff, trans_bw, revfilt } );
   return;
end;

if EEG.trials == 1 
	if ~isempty(EEG.event) & isfield(EEG.event, 'type') & isstr(EEG.event(1).type)
        tmpevent = EEG.event;    
		boundaries = strmatch('boundary', {tmpevent.type});
		if isempty(boundaries)
            [EEG.data b a] = iirfilt( EEG.data, options{:}); 
		else
			options{4} = 0;
			disp('pop_iirfilt:finding continuous data boundaries');
			tmplat = cell2mat({tmpevent.latency});
            boundaries = tmplat(boundaries);
            boundaries = [0 round(boundaries-0.5) EEG.pnts];
            try, warning off MATLAB:divideByZero
            catch, end;
			for n=1:length(boundaries)-1
				try
                    fprintf('Processing continuous data (%d:%d)\n',boundaries(n),boundaries(n+1)); 
                    [ EEG.data(:,boundaries(n)+1:boundaries(n+1)) b a] = ...
                        iirfilt(EEG.data(:,boundaries(n)+1:boundaries(n+1)), options{:});
				catch
					fprintf('\nFilter error: continuous data portion too narrow (DC removed if highpass only)\n');
                    if locutoff ~= 0 & hicutoff == 0
                        tmprange = [boundaries(n)+1:boundaries(n+1)];
                        EEG.data(:,tmprange) = ...
                            EEG.data(:,tmprange) - repmat(mean(EEG.data(:,tmprange),2), [1 length(tmprange)]);
                    end;
				end;
			end
            try, warning on MATLAB:divideByZero
            catch, end;
		end
	else
        [EEG.data b a] = iirfilt( EEG.data, options{:});
	end;
else
    EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials);
    [EEG.data b a] = iirfilt( EEG.data, options{:});
	% Note: reshape does not reserve new memory while EEG.data(:,:) does
end;	

if plotfreqz && exist('b') == 1 && exist('a') == 1 && (~isequal(a,0) || ~isequal(b,0))
    freqz(b, a, [], EEG.srate);
elseif plotfreqz
    disp('Cannot plot frequency response of band pass or notch filter');
end

com = sprintf( '%s = pop_iirfilt( %s, %s, %s, [%s], %s, %s);', inputname(1), inputname(1), ...
			num2str( locutoff), num2str( hicutoff), num2str( trans_bw), num2str( revfilt), num2str( causal) );
