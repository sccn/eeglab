% pop_resample() - resample dataset (pop up window).
%
% Usage:
%   >> [OUTEEG] = pop_resample( INEEG ); % pop up interactive window
%   >> [OUTEEG] = pop_resample( INEEG, freq);
%   >> [OUTEEG] = pop_resample( INEEG, freq, fc, df);
%
% Graphical interface:
%   The edit box entitled "New sampling rate" contains the frequency in
%   Hz for resampling the data. Entering a value in this box  is the same 
%   as providing it in the 'freq' input from the command line.
%
% Inputs:
%   INEEG      - input dataset
%   freq       - frequency to resample (Hz)  
%
% Optional inputs:
%   fc         - anti-aliasing filter cutoff (pi rad / sample)
%                {default 0.9}
%   df         - anti-aliasing filter transition band width (pi rad /
%                sample) {default 0.2}
%
% Outputs:
%   OUTEEG     - output dataset
%
% Author: Arnaud Delorme, CNL/Salk Institute, 2001
%
% Note: uses the resample() function from the signal processing toolbox
%       if present. Otherwise use griddata interpolation method (it should be
%       reprogrammed using spline interpolation for speed up).
%
% See also: resample(), eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% 01-25-02 reformated help & license -ad 
% 03-08-02 remove ica activity resampling (now set to []) -ad
% 03-08-02 debug call to function help -ad
% 04-05-02 recompute event latencies -ad

function [EEG, command] = pop_resample( EEG, freq, fc, df); 

command = '';
if nargin < 1
    help pop_resample;
    return;
end    
if isempty(EEG(1).data)
    disp('Pop_resample error: cannot resample empty dataset'); return;
end

if ~plugin_askinstall('Firfilt', 'pop_firwsord'), return; end

if nargin < 2 

	% popup window parameters
	% -----------------------
	promptstr    = {['New sampling rate']};
	inistr       = { num2str(EEG(1).srate) };
	result       = inputdlg2( promptstr, 'Resample current dataset -- pop_resample()', 1,  inistr, 'pop_resample');
	if length(result) == 0 return; end
	freq         = eval( result{1} );

end

% Default cutoff frequency (pi rad / smp)
if nargin < 3 || isempty(fc)
    fc = 0.9;
end
% Default transition band width (pi rad / smp)
if nargin < 4 || isempty(df)
    df = 0.2;
end
% fc in range?
if fc < 0 || fc > 1
    error('Anti-aliasing filter cutoff freqeuncy out of range.')
end

% process multiple datasets
% -------------------------
if length(EEG) > 1
    if nargin < 2
        [ EEG command ] = eeg_eval( 'pop_resample', EEG, 'warning', 'on', 'params', { freq } );
    else
        [ EEG command ] = eeg_eval( 'pop_resample', EEG, 'params', { freq } );
    end
    return;
end

% finding the best ratio
% [p,q] = rat(freq/EEG.srate, 0.0001) % not used right now 
[p,q] = rat(freq/EEG.srate, 1e-12); % used! AW

% set variable
% ------------
EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
oldpnts  = EEG.pnts;

% resample for multiple channels
% -------------------------
if isfield(EEG, 'event') && isfield(EEG.event, 'type') && ischar(EEG.event(1).type)
    tmpevent = EEG.event;
    bounds = strmatch('boundary', { tmpevent.type });
    if ~isempty(bounds),
        disp('Data break detected and taken into account for resampling');
        bounds = [ tmpevent(bounds).latency ];
        bounds(bounds <= 0 | bounds > size(EEG.data,2)) = []; % Remove out of range boundaries
        bounds(mod(bounds, 1) ~= 0) = round(bounds(mod(bounds, 1) ~= 0) + 0.5); % Round non-integer boundary latencies
    end
    bounds = [1 bounds size(EEG.data, 2) + 1]; % Add initial and final boundary event
    bounds = unique(bounds); % Sort (!) and remove doublets
else 
    bounds = [1 size(EEG.data,2) + 1]; % [1:size(EEG.data,2):size(EEG.data,2)*size(EEG.data,3)+1];
end

eeglab_options;
if option_donotusetoolboxes
    usesigproc = 0;
elseif exist('resample') == 2
     usesigproc = 1;
else usesigproc = 0;
    disp('Signal Processing Toolbox absent: using custom interpolation instead of resample() function.');
    disp('This method uses cubic spline interpolation after anti-aliasing (see >> help spline)');    
end

fprintf('resampling data %3.4f Hz\n', EEG.srate*p/q);
eeglab_options;
for index1 = 1:size(EEG.data,1)
    fprintf('%d ', index1);	
    sigtmp = reshape(EEG.data(index1,:, :), oldpnts, EEG.trials);
    
    if index1 == 1
        tmpres = [];
        indices = [1];
        for ind = 1:length(bounds)-1
            tmpres  = [ tmpres; myresample( double( sigtmp(bounds(ind):bounds(ind+1)-1,:)), p, q, usesigproc, fc, df ) ];
            indices = [ indices size(tmpres,1)+1 ];
        end
        if size(tmpres,1) == 1, EEG.pnts  = size(tmpres,2);
        else                    EEG.pnts  = size(tmpres,1);
        end
        if option_memmapdata == 1
             tmpeeglab = mmo([], [EEG.nbchan, EEG.pnts, EEG.trials]);
        else tmpeeglab = zeros(EEG.nbchan, EEG.pnts, EEG.trials);
        end
    else
        for ind = 1:length(bounds)-1
            tmpres(indices(ind):indices(ind+1)-1,:) = myresample( double( sigtmp(bounds(ind):bounds(ind+1)-1,:) ), p, q, usesigproc, fc, df);
        end
    end; 
    tmpeeglab(index1,:, :) = tmpres;
end
fprintf('\n');	
EEG.srate   = EEG.srate*p/q;
EEG.data = tmpeeglab;
EEG.pnts    = size(EEG.data,2);
EEG.xmax    = EEG.xmin + (EEG.pnts-1)/EEG.srate; % cko: recompute xmax, since we may have removed a few of the trailing samples
EEG.times   = linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts);

% recompute all event latencies
% -----------------------------
if isfield(EEG.event, 'latency')
    fprintf('resampling event latencies...\n');
    
    if EEG.trials > 1 % Epoched data; not recommended
        
        warning( 'Resampling of epoched data is not recommended (due to anti-aliasing filtering)! Note: For epoched datasets recomputing urevent latencies is not supported. The urevent structure will be cleared.' )
        
        for iEvt = 1:length( EEG.event )
            
            EEG.event( iEvt ).latency = ( EEG.event( iEvt ).latency - ( EEG.event(iEvt).epoch - 1 ) * oldpnts - 1 ) * p / q + ( EEG.event(iEvt).epoch - 1 ) * EEG.pnts + 1;
            
            % % Recompute event latencies
            if isfield(EEG.event, 'duration') && ~isempty(EEG.event(iEvt).duration)
                EEG.event( iEvt ).duration = EEG.event( iEvt ).duration * p/q;
            end
        end
        
        EEG.urevent = [];
        
    else % Continuous data

        for iEvt = 1:length(EEG.event)

            % From >> help resample: Y is P/Q times the length of X (or the
            % ceiling of this if P/Q is not an integer).
            % That is, recomputing event latency by pnts / oldpnts will give
            % inaccurate results in case of multiple segments and rounded segment
            % length. Error is accumulated and can lead to several samples offset.
            % Blocker for boundary events.
            % Old version EEG.event(index1).latency = EEG.event(index1).latency * EEG.pnts /oldpnts;

            % Recompute event latencies relative to segment onset
            if strcmpi(EEG.event(iEvt).type, 'boundary') && mod(EEG.event(iEvt).latency, 1) == 0.5 % Workaround to keep EEGLAB style boundary events at -0.5 latency relative to DC event; actually incorrect
                iBnd = sum(EEG.event(iEvt).latency + 0.5 >= bounds);
                EEG.event(iEvt).latency = indices(iBnd) - 0.5;
            else
                iBnd = sum(EEG.event(iEvt).latency >= bounds);
                EEG.event(iEvt).latency = (EEG.event(iEvt).latency - bounds(iBnd)) * p / q + indices(iBnd);
            end
            
            % Recompute event duration relative to segment onset
            if isfield(EEG.event, 'duration') && ~isempty(EEG.event(iEvt).duration)
                EEG.event( iEvt ).duration = EEG.event( iEvt ).duration * p/q;
            end

        end

        if isfield(EEG, 'urevent') && isfield(EEG.urevent, 'latency')
            try
                for iUrevt = 1:length(EEG.urevent)
                    % Recompute urevent latencies relative to segment onset
                    if strcmpi(EEG.urevent(iUrevt).type, 'boundary') && mod(EEG.urevent(iUrevt).latency, 1) == 0.5 % Workaround to keep EEGLAB style boundary events at -0.5 latency relative to DC event; actually incorrect
                        iBnd = sum(EEG.urevent(iUrevt).latency + 0.5 >= bounds);
                        EEG.urevent(iUrevt).latency = indices(iBnd) - 0.5;
                    else
                        iBnd = sum(EEG.urevent(iUrevt).latency >= bounds);
                        EEG.urevent(iUrevt).latency = (EEG.urevent(iUrevt).latency - bounds(iBnd)) * p / q + indices(iBnd);
                    end

                end
            catch
                disp('pop_resample warning: ''urevent'' problem, reinitializing urevents');
                EEG = rmfield(EEG, 'urevent');
            end
        end
    
    end
    
    EEG = eeg_checkset(EEG, 'eventconsistency');
end

% resample for multiple channels ica
EEG.icaact = [];

% store dataset
fprintf('resampling finished\n');

EEG.setname = [EEG.setname ' resampled'];

command = sprintf('EEG = pop_resample( EEG, %d);', freq);
return;

% resample if resample is not present
% -----------------------------------
function tmpeeglab = myresample(data, p, q, usesigproc, fc, df)
    
    if length(data) < 2
        tmpeeglab = data;
        return;
    end
    %if size(data,2) == 1, data = data'; end
    if usesigproc
        % padding to avoid artifacts at the beginning and at the end
        % Andreas Widmann May 5, 2011
        
        %The pop_resample command introduces substantial artifacts at beginning and end
        %of data when raw data show DC offset (e.g. as in DC recorded continuous files)
        %when MATLAB Signal Processing Toolbox is present (and MATLAB resample.m command
        %is used).
        %Even if this artifact is short, it is a filtered DC offset and will be carried
        %into data, e.g. by later highpass filtering to a substantial amount (easily up
        %to several seconds).
        %The problem can be solved by padding the data at beginning and end by a DC
        %constant before resampling.

%         N = 10; % Resample default
%         nPad = ceil((max(p, q) * N) / q) * q; % # datapoints to pad, round to integer multiple of q for unpadding
%         tmpeeglab = resample([data(ones(1, nPad), :); data; data(end * ones(1, nPad), :)], pnts, new_pnts);

        % Conservative custom anti-aliasing FIR filter, see bug 1757
        nyq = 1 / max([p q]);
        fc = fc * nyq; % Anti-aliasing filter cutoff frequency
        df = df * nyq; % Anti-aliasing filter transition band width
        m = pop_firwsord('kaiser', 2, df, 0.002); % Anti-aliasing filter kernel
        b = firws(m, fc, windows('kaiser', m + 1, 5)); % Anti-aliasing filter kernel
        b = p * b; % Normalize filter kernel to inserted zeros
%         figure; freqz(b, 1, 2^14, q * 1000) % Debugging only! Sampling rate hardcoded as it is unknown in this context. Manually adjust for debugging!

        % Padding, see bug 1017
        nPad = ceil((m / 2) / q) * q; % Datapoints to pad, round to integer multiple of q for unpadding
        startPad = repmat(data(1, :), [nPad 1]);
        endPad = repmat(data(end, :), [nPad 1]);

        % Resampling
        tmpeeglab = resample([startPad; data; endPad], p, q, b);

        % Remove padding
        nPad = nPad * p / q; % # datapoints to unpad
        tmpeeglab = tmpeeglab(nPad + 1:end - nPad, :); % Remove padded data

    else % No Signal Processing toolbox
        
        % anti-alias filter
        % -----------------
%        data         = eegfiltfft(data', 256, 0, 128*pnts/new_pnts); % Downsample from 256 to 128 times the ratio of freq. 
%                                                                      % Code was verified by Andreas Widdman March  2014

%                                                                      % No! Only cutoff frequency for downsampling was confirmed.
%                                                                      % Upsampling doesn't work and FFT filtering introduces artifacts.
%                                                                      % Also see bug 1757. Replaced May 05, 2015, AW

        if p < q, nyq = p / q; else nyq = q / p; end
        fc = fc * nyq; % Anti-aliasing filter cutoff frequency
        df = df * nyq; % Anti-aliasing filter transition band width
        m = pop_firwsord('kaiser', 2, df, 0.002); % Anti-aliasing filter kernel
        b = firws(m, fc, windows('kaiser', m + 1, 5)); % Anti-aliasing filter kernel
%         figure; freqz(b, 1, 2^14, 1000) % Debugging only! Sampling rate hardcoded as it is unknown in this context. Manually adjust for debugging!

        if p < q % Downsampling, anti-aliasing filter
            data = firfiltdcpadded(b, data, 0);
        end

        % spline interpolation
        % --------------------
%         X            = [1:length(data)];
%         nbnewpoints  = length(data)*p/q;
%         nbnewpoints2 = ceil(nbnewpoints);
%         lastpointval = length(data)/nbnewpoints*nbnewpoints2;        
%         XX = linspace( 1, lastpointval, nbnewpoints2);

        % New time axis scaling, May 06, 2015, AW
        X = 0:length(data) - 1;
        newpnts  = ceil(length(data) * p / q);
        XX = (0:newpnts - 1) / (p / q);

        cs = spline( X, data);
        tmpeeglab = ppval(cs, XX)';

        if p > q % Upsampling, anti-imaging filter
            tmpeeglab = firfiltdcpadded(b, tmpeeglab, 0);
        end

    end
    
% firfiltdcpadded() - Pad data with DC constant and filter
%
% Usage:
%   >> data = firfiltdcpadded(data, b, causal);
%
% Inputs:
%   data      - raw data
%   b         - vector of filter coefficients
%   causal    - boolean perform causal filtering {default 0}
%
% Outputs:
%   data      - smoothed data
%
% Note:
%   firfiltdcpadded always operates (pads, filters) along first dimension.
%   Not memory optimized.
%
% Author: Andreas Widmann, University of Leipzig, 2013
%
% See also:
%   firfiltsplit

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2013 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

function [ data ] = firfiltdcpadded(b, data, causal)

% Defaults
if nargin < 3 || isempty(causal)
    causal = 0;
end

% Check arguments
if nargin < 2
    error('Not enough input arguments.');
end

% Filter's group delay
if mod(length(b), 2) ~= 1
    error('Filter order is not even.');
end
groupDelay = (length(b) - 1) / 2;
b = double(b); % Filter with double precision

% Pad data with DC constant
if causal
    startPad = repmat(data(1, :), [2 * groupDelay 1]);
    endPad = [];
else
    startPad = repmat(data(1, :), [groupDelay 1]);
    endPad = repmat(data(end, :), [groupDelay 1]);
end

% Filter data
data = filter(b, 1, double([startPad; data; endPad])); % Pad and filter with double precision

% Remove padded data
data = data(2 * groupDelay + 1:end, :);
