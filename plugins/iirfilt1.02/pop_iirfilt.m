% pop_iirfilt() - interactively filter EEG dataset data using iirfilt()
%
% Usage:
%   >> EEGOUT = pop_iirfilt( EEG, locutoff, hicutoff, trans_bw);
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
if isempty(EEG.data)
    disp('pop_iirfilt() error: cannot filter an empty dataset'); return;
end;    
if nargin < 2
	% which set to save
	% -----------------
   	promptstr = { 'Highpass: low edge of the frequency pass band (Hz) (0 -> lowpass)', ...
   				  'Lowpass: high edge of the frequency pass band (Hz) (0 -> highpass)', ...
   				  strvcat('Notch filter the data. Give the notch range, i.e. [45 55] for 50 Hz)', ...
                  '(this option overwrites the low and high edge limits given above)'), ...
                  'Transition BW filter (default: see >> help pop_iirfilt)' };
	inistr       = { '0', '0', '', '' };
   	result       = inputdlg2( promptstr, 'Filter the data -- pop_iirfilt()', 1,  inistr, 'pop_iirfilt');
	if size(result, 1) == 0 return; end;
	locutoff   	 = eval( result{1} );
	hicutoff 	 = eval( result{2} );
	if isempty( result{3} )
		 revfilt = 0;
	else 
        revfilt    = eval( [ '[' result{3} ']' ] );
        locutoff = revfilt(1);
        hicutoff = revfilt(2);
        revfilt = 1;
	end;
	if locutoff == 0 & hicutoff == 0 return; end;
	if isempty( result{4} )
		 trans_bw = [];
	else trans_bw    = eval( result{4} );
	end;
    causal = 'off';
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
        causal = 'off';
    end;
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
	options = { options{:} 0 [] [] causal };
end;

% warning
% -------
if exist('filtfilt') ~= 2
    disp('Warning: cannot find the signal processing toolbox');
    disp('         a simple fft/inverse fft filter will be used');
end;

if EEG.trials == 1 
	if ~isempty(EEG.event) & isfield(EEG.event, 'type') & isstr(EEG.event(1).type)
		boundaries = strmatch('boundary', {EEG.event.type});
		if isempty(boundaries)
            EEG.data = iirfilt( EEG.data, options{:}); 
		else
			options{4} = 0;
			disp('pop_iirfilt:finding continuous data boundaries');
			tmplat = cell2mat({EEG.event.latency});
            boundaries = tmplat(boundaries);
            boundaries = [0 round(boundaries-0.5) EEG.pnts];
            try, warning off MATLAB:divideByZero
            catch, end;
			for n=1:length(boundaries)-1
				try
                    fprintf('Processing continuous data (%d:%d)\n',boundaries(n),boundaries(n+1)); 
                    EEG.data(:,boundaries(n)+1:boundaries(n+1)) = ...
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
        EEG.data = iirfilt( EEG.data, options{:});
	end;
else
    EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials);
    EEG.data = iirfilt( EEG.data, options{:});
	% Note: reshape does not reserve new memory while EEG.data(:,:) does
end;	


com = sprintf( '%s = pop_iirfilt( %s, %s, %s, [%s], [%s]);', inputname(1), inputname(1), ...
			num2str( locutoff), num2str( hicutoff), num2str( trans_bw ), num2str( revfilt ) );
return;
