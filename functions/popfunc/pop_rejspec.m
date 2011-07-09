% pop_rejspec() - rejection of artifact in a dataset using 
%                 thresholding of frequencies in the data.
% Usage:
%   >>  pop_rejspec(INEEG, typerej); % pop-up interactive windo mode
%   >> [OUTEEG, Indices] = pop_rejspec( INEEG, typerej, elec_comp, ...
%         lowthresh, upthresh, startfreq, endfreq, superpose, reject);
%
% Pop-up window options:
%   "Electrode|Component" - [edit box] electrode or component number(s) to 
%                 take into consideration for rejection. Sets the 'elec_comp'
%                 parameter in the command line call (see below).
%   "Lower limits(s)" - [edit box] lower threshold limits(s) (in dB). 
%                 Sets the command line parameter 'lowthresh'. If more than
%                 one, apply to each electrode|component individually. If
%                 fewer than number of electrodes|components, apply the
%                 last values to all remaining electrodes|components.
%   "Upper limits(s)" - [edit box] upper threshold limit(s) in dB. 
%                 Sets the command line parameter 'upthresh'.
%   "Low frequency(s)" - [edit box] low-frequency limit(s) in Hz. 
%                 Sets the command line parameter 'startfreq'.
%   "High frequency(s)" - [edit box] high-frequency limit(s) in Hz. 
%                 Sets the command line parameter 'endfreq'.
%   "Display previous rejection marks?" - [edit box] either YES or NO. 
%                  Sets the command line input option 'superpose'.
%   "Reject marked trials?" - [edit box] either YES or NO. Sets the
%                 command line input option 'reject'.
%
% Command line inputs:
%   INEEG      - input dataset
%   typerej    - [1|0] data to reject on (0 = component activations; 1 = 
%              electrode data). {Default is 1}. 
%   elec_comp  - [e1 e2 ...] electrode|component number(s) to take into 
%              consideration during rejection
%   lowthresh  - lower threshold limit(s) in dB. Can be an array if 
%              several electrodes|components. If fewer values than number 
%              of electrodes|components, the last value is used for the 
%              remaining electrodes|components.
%   upthresh  - upper threshold limit(s) in dB (same syntax as lowthresh)
%   startfreq  - low frequency limit(s) in Hz (same syntax  as lowthresh)
%   endfreq    - high frequency limit(s) in Hz (same syntax  as lowthresh).
%              Options 'startfreq' and 'endfreq' define the frequncy range 
%              used during rejection.
%   superpose  - [0|1] 0 = Do not superpose rejection marks on previous
%              marks stored in the dataset. 1 = Show both previous and
%              current marks using different colors. {Default: 0}.
%   reject     - [0|1] 0 = Do not reject marked trials (but store the 
%              marks. 1 = Reject marked trials. {Default: 1}.
%
% Outputs:
%   OUTEEG     - output dataset with updated spectrograms
%   Indices    - index of rejected trials
%   Note: When eegplot() is called, modifications are applied to the current 
%   dataset at the end of the call to eegplot() (e.g., when the user presses 
%   the 'Reject' button).
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eegthresh(), eeglab(), eegplot(), pop_rejepoch()

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
% 03-07-02 added srate argument to eegplot call -ad
% 03-08-02 reworked spectrum to save space & add eeglab options -ad

function [EEG, Irej, com] = pop_rejspec( EEG, icacomp, varargin);
    %elecrange, negthresh, posthresh, ...
   	%startfreq, endfreq, superpose, reject);

Irej = [];
com  = '';
if nargin < 1
   help pop_rejspec;
   return;
end;  
if nargin < 2
   icacomp = 1;
end;  

if icacomp == 0
	if isempty( EEG.icasphere )
	    ButtonName=questdlg( 'Do you want to run ICA now ?', ...
                         'Confirmation', 'NO', 'YES', 'YES');
    	switch ButtonName,
        	case 'NO', disp('Operation cancelled'); return;   
        	case 'YES', [ EEG com ] = pop_runica(EEG);
    	end % switch
	end;
end;	
if exist('reject') ~= 1
    reject = 1;
end;

if nargin < 3

	% which set to save
	% -----------------
	promptstr   = { fastif(icacomp==0, 'Component number(s) (Ex: 2 4 5):', ...
                                           'Electrode number(s) (Ex: 2 4 5):'), ...
					'Lower limit(s) (dB):', ...
					'Upper limit(s) (dB):', ...
					'Low frequency(s) (Hz):', ...
					'High frequency(s) (Hz):', ...
               		'Display previous rejection marks? (YES or NO)', ...
         			'Reject marked trial(s)? (YES or NO)' };
	inistr      = { ['1:' int2str(EEG.nbchan)], ...
					'-30', ...
					'30', ...
					'15', ...
					'30', ...
               		'NO', ...
            		'NO' };

	result       = inputdlg2( promptstr, fastif(~icacomp, 'Reject by component spectra -- pop_rejspec()', ...
											   'Reject by data spectra -- pop_rejspec()'), 1,  inistr, 'pop_rejspec');
	size_result  = size( result );
	if size_result(1) == 0 return; end;
    options = {};
    options = { options{:} 'elecrange' eval( [ '[' result{1} ']' ]  ) };
    thresholdsLow  = eval( [ '[' result{2} ']' ]  );
    thresholdsHigh = eval( [ '[' result{3} ']' ]  );
    options = { options{:} 'threshold' [thresholdsLow(:) thresholdsHigh(:) ] };
    freqLimitsLow   = eval( [ '[' result{4} ']' ]  );
    freqLimitsHigh  = eval( [ '[' result{5} ']' ]  );
    options = { options{:} 'freqlimits' [freqLimitsLow(:) freqLimitsHigh(:) ] };
    
	switch lower(result{6}), case 'yes', superpose=1; otherwise, superpose=0; end;
	switch lower(result{7}), case 'yes', reject=1; otherwise, reject=0; end;
    options = { options{:} 'eegplotplotallrej' superpose };
    options = { options{:} 'eegplotreject'     reject };
else
    if isnumeric(varargin{3}) || ~isempty(str2num(varargin{3}))
        options = {};
        if isstr(varargin{1}), varargin{1} = str2num(varargin{1}); end;
        if isstr(varargin{2}), varargin{2} = str2num(varargin{2}); end;
        if isstr(varargin{3}), varargin{3} = str2num(varargin{3}); end;
        if isstr(varargin{4}), varargin{4} = str2num(varargin{4}); end;
        if isstr(varargin{5}), varargin{5} = str2num(varargin{5}); end;
        if nargin > 2, options = { options{:} 'elecrange'   varargin{1} }; end;
        if nargin > 3, options = { options{:} 'threshold'   [ varargin{2}; varargin{3}]' }; end;
        if nargin > 5, options = { options{:} 'freqlimits'  [ varargin{4}; varargin{5}]' }; end;
        if nargin > 7, options = { options{:} 'eegplotplotallrej' varargin{6}  }; end;
        if nargin > 8, options = { options{:} 'eegplotreject'     varargin{7}  }; end;
        if nargin > 9, options = { options{:} 'eegplotcom'        varargin{8}  }; end;
    else
        options = varargin;
    end;
end;

opt = finputcheck( options, { 'elecrange'     'integer'  []    [1:EEG.nbchan];
                              'threshold'     'real'     []    [-30 30];
                              'freqlimits'    'real'     []    [15 30];
                              'specdata'      'real'     []    EEG.specdata;
                              'eegplotcom'    'string'   []    '';
                              'method'        'string'   { 'fft';'multitaper' }    'multitaper';
                              'eegplotreject' 'integer'  []    0;
                              'eegplotplotallrej' 'integer'  []    0 }, 'pop_rejspec');
if isstr(opt), error(opt); end;

sizewin = 2^nextpow2(EEG.pnts);
if icacomp == 1
    [allspec, Irej, tmprejE, freqs ] = spectrumthresh( EEG.data, opt.specdata, ...
							opt.elecrange, EEG.srate, opt.threshold(:,1)', opt.threshold(:,2)', opt.freqlimits(:,1)', opt.freqlimits(:,2)', opt.method);
    rejE = zeros(EEG.nbchan, EEG.trials);
    rejE(opt.elecrange,Irej) = tmprejE;
else
    % test if ICA was computed
    % ------------------------
    icaacttmp = eeg_getdatact(EEG, 'component', [1:size(EEG.icaweights,1)]);
    [allspec, Irej, tmprejE, freqs ] = spectrumthresh( icaacttmp, EEG.specicaact, ...
                            opt.elecrange, EEG.srate, opt.threshold(:,1)', opt.threshold(:,2)', opt.freqlimits(:,1)', opt.freqlimits(:,2)', opt.method);
    rejE = zeros(EEG.nbchan, size(icaacttmp,1));
    rejE(opt.elecrange,Irej) = tmprejE;
end;

fprintf('%d channel selected\n', size(opt.elecrange(:), 1));
fprintf('%d/%d trials marked for rejection\n', length(Irej), EEG.trials);
rej = zeros( 1, EEG.trials);
rej(Irej) = 1;

if nargin < 3 || opt.eegplotplotallrej == 2
	nbpnts = size(allspec,2);
    if icacomp == 1 macrorej  = 'EEG.reject.rejfreq';
        			macrorejE = 'EEG.reject.rejfreqE';
    else			macrorej  = 'EEG.reject.icarejfreq';
        			macrorejE = 'EEG.reject.icarejfreqE';
    end;
	colrej = EEG.reject.rejfreqcol;
    
    elecrange  = opt.elecrange;
    superpose  = opt.eegplotplotallrej;
    reject     = opt.eegplotreject;
    topcommand = opt.eegplotcom;
	eeg_rejmacro; % script macro for generating command and old rejection arrays
    
	if icacomp == 1
		eegplot(EEG.data(opt.elecrange,:,:), 'winlength', 5, 'position', [100 550 800 500], ...
			'limits', [EEG.xmin EEG.xmax]*1000, 'xgrid', 'off', 'tag', 'childEEG' );
	else
		eegplot(icaacttmp(opt.elecrange,:,:), 'winlength', 5, 'position', [100 550 800 500], 'limits', ...
				[EEG.xmin EEG.xmax]*1000 , 'xgrid', 'off', 'tag', 'childEEG' );
	end;	
	eegplot( allspec(elecrange,:,:), 'srate', EEG.srate, 'freqlimits', [1 EEG.srate/2], 'command', ...
			 command, 'children', findobj('tag', 'childEEG'), 'position', [100 50 800 500], eegplotoptions{:}); 
end;
if ~isempty(rej)
	if icacomp	== 1
		EEG.reject.rejfreq = rej;
		EEG.reject.rejfreqE = rejE;
	else
		EEG.reject.icarejfreq = rej;
		EEG.reject.icarejfreqE = rejE;
	end;
    if reject
        EEG = pop_rejepoch(EEG, rej, 0);
    end;
    Irej = find(rej);
end;

% store variables
% ---------------
if icacomp == 1, EEG.specdata = allspec;
else,            EEG.specicaact = allspec;
end;
    
com = [com sprintf('%s = pop_rejspec( %s, %s);', inputname(1), ...
   inputname(1), vararg2str({icacomp, 'elecrange', opt.elecrange, 'threshold', opt.threshold, 'freqlimits', opt.freqlimits, ...
     'eegplotcom', opt.eegplotcom, 'eegplotplotallrej' opt.eegplotplotallrej 'eegplotreject' opt.eegplotreject })) ]; 

return;

% compute spectrum and reject artifacts
% -------------------------------------
function [specdata, Irej, Erej, freqs ] = spectrumthresh( data, specdata, elecrange, srate, negthresh, posthresh, startfreq, endfreq, method);

	% compute the fft if necessary - old version
    if isempty(specdata)
        if strcmpi(method, 'fft')
            sizewin = size(data,2);
            freqs = srate*[1, sizewin]/sizewin/2;
            specdata = fft( data-repmat(mean(data,2), [1 size(data,2) 1]), sizewin, 2);
            specdata = specdata( :, 2:sizewin/2+1, :);
            specdata = 10*log10(abs( specdata ).^2);
            specdata  = specdata - repmat( mean(specdata,3), [1 1 size(data,3)]);
        else
            if ~exist('pmtm')
                error('The signal processing toolbox needs to be installed');
            end;
            [tmp freqs] = pmtm( data(1,:,1), [],[],srate); % just to get the frequencies 	

            fprintf('Computing spectrum (using slepian tapers; done only once):\n');

            for index = 1:size(data,1)
                fprintf('%d ', index);    
                for indextrials = 1:size(data,3)
                    [ tmpspec(index,:,indextrials) freqs] = pmtm( data(index,:,indextrials) , [],[],srate);
                end;
            end;
            tmpspec  = 10*log(tmpspec);
            tmpspec  = tmpspec - repmat( mean(tmpspec,3), [1 1 size(data,3)]);
            specdata = tmpspec;
        end;
    else
        if strcmpi(method, 'fft')
            sizewin = size(data,2);
            freqs = srate*[1, sizewin]/sizewin/2;
        else
            [tmp freqs] = pmtm( data(1,:,1), [],[],srate); % just to get the frequencies 	
        end;
    end;
    
	% perform the rejection
	% ---------------------	
	[I1 Irej NS Erej] = eegthresh( specdata(elecrange, :, :), size(specdata,2), 1:length(elecrange), negthresh, posthresh, ...
										 [freqs(1) freqs(end)], startfreq, min(freqs(end), endfreq));
	fprintf('\n');    
	
