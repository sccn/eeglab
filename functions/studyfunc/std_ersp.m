% std_ersp() - Computes ERSP and/or ITC information for ICA components of a dataset, 
%              saves results into float files and and places pointers to the files 
%              in the dataset EEG structures. When the ERSP/ITC float files already exist, 
%              the function loads the ERSP/ITC information from it, unless the requested 
%              flag specifies differently: If so, a query window pops up. 
%
%              Re options to specify component numbers, the desired frequency range, 
%              timewindow, resolution, confidence level and the wavelet parameters 
%              (number of cycles): see >> help timef and >> timef details 
%
%              The function returns the masked (as per the requested alpha) ERSP or ITC
%              for the selected ICA components in the requested frequency range and 
%              time window (the two are dependent). Frequencies are equally log spaced.
%
%              Two float files are saved, one for ERSP and one for ITC information:
%              These contain the unprocessed ERSP|ITC image and its significant levels.
%              (See filename info below).
%
%              If the ERSPs/ITCs were previously saved in these files and the same set of
%              ERSP/ITC parameters are used, the values are not recomputed.
%              Vectors of frequencies and latencies for the ERSP/ITC images are returned 
%              separately, as well as the EEG.etc sub-structure modified with pointers 
%              to the output float files and some information about them. 
% Usage:  
%              >> [EEG_etc X times logfreqs ] = std_ersp(EEG, components,  ...
%                                                    freqrange, timewindow,  ...
%                                                         cycles, padratio, alpha,...
%                                                              type, powbase);
%
%              % If they do not exist, computes the ICA component ERSPs/ITCs 
%              % for a dataset. Saves the computed images in dataset name files 
%              % with extensions .icaersp and .icaitc
%              % Also updates the EEG structure in the Matlab environment 
%              % and saves the modified dataset to disk.
% Inputs:
%
%   EEG        - an EEG data structure. 
%   components - [numeric vector] of the EEG structure for which an ERSP   
%                 and ITC will be computed {default|[]: all components}
%   freqrange  - [minHz maxHz] the frequency range to compute the ERSP/ITC.
%   cycles     - If 0 -> Use FFTs (with constant window length) {default: 0}
%                 If >0 -> Number of cycles in each analysis wavelet 
%                 If [wavecycles factor] -> wavelet cycles increase with frequency 
%                 beginning at wavecyles (0 < factor < 1; factor = 1 -> no increase,
%                 standard wavelets; factor=0 -> fixed epoch length, as in FFT.
%   padratio   - FFT-length/winframes (2^k)                          
%                 Multiplies the number of output frequencies by dividing their spacing. 
%                 When cycles==0, frequency spacing is (low_freq/padratio).
%   alpha      - If in (0, 1), compute two-tailed permutation-based prob. thresholds
%                 and use these to mask the ERSP/ITC image (in output X)
%   type       - ['ersp'|'itc'] both ERSP and ITC images are computed and saved to disk, 
%                 but only one is returned (output X) {default: 'ersp'}
%   powbase    - [vector] baseline spectrum (power not dB, output from timef() run with 
%                 same parameters as above) to use in the timef() computation 
%                 {default|[] -> pre-time 0}
% Outputs:
%
%   EEG_etc   - the EEG dataset .etc sub-structure (i.e. EEG.etc), which is
%               modified with the output filenames and some information about
%               the float files that hold the dataset ERSP and/or ITC  
%               information. If an ERSP/ITC file already exists and wasn't
%               recomputed and this output will be empty. 
%   X         - the masked log ERSP/ITC of the requested ICA components in the 
%               selected frequency and time range. 
%   times     - a vector of time points for which the ERSPs/ITCs were computed. 
%   logfreqs  - a vector of (equally log spaced) frequencies (in Hz) at which the 
%               log ERSP/ITC was evaluated. 
%
% Files written or modified:     [dataset_filename].icaersp   <-- saved component ERSPs
%                                [dataset_filename].icaitc    <-- saved component ITCs
%                                [dataset_filename].set       <-- re-saved dataset
% Example: 
%            % create ERSP and ITC images on disk for all comps from dataset EEG
%            % use three-cycle wavelets (at 3 Hz) to >3-cycle wavelets at 50 Hz
%            % use probability masking at p < 0.01, padratio 4. See >> timef details
%            % returns log-freq spaced, probability-masked Xersp
%            >> [EEG_etc, Xersp, times, logfreqs] = std_ersp(EEG, ...
%                                                    [1:size(EEG.icawinv,2)],...
%                                                     [3 50], [3 0.5], 4, 0.01, 'ersp');
%
% See also: timef(), std_itc(), std_erp(), std_spec(), std_map(), std_preclust()
%
% Authors:  Hilit Serby, Arnaud Delorme, SCCN, INC, UCSD, January, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, October 11, 2004, hilit@sccn.ucsd.edu
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

% $Log: not supported by cvs2svn $
% Revision 1.18  2006/03/08 22:50:40  arno
% fix last change
%
% Revision 1.17  2006/03/08 22:49:10  arno
% detect the presence of file or not
%
% Revision 1.16  2006/03/08 20:29:08  arno
% rename func
%
% Revision 1.15  2006/03/08 19:43:37  scott
% fixed powbase definition -sm & ad
%
% Revision 1.14  2006/03/08 03:02:27  scott
% expand powbase to a matrix -sm
%
% Revision 1.13  2006/03/08 02:58:08  scott
% debug3
%
% Revision 1.12  2006/03/08 02:54:41  scott
% debug 2
%
% Revision 1.11  2006/03/08 02:38:14  scott
% debug
%
% Revision 1.10  2006/03/07 23:03:47  scott
% added optional posbase argument for timef(). NOTE: probability masking now
% allowed and computed
% but the masked version is not output by timef (I think - must clarify timef help)
% and therefore not used here... -sm
%
% Revision 1.9  2006/03/07 19:18:44  arno
% header
%
% Revision 1.8  2006/03/07 03:58:32  scott
% same -sm
%
% Revision 1.7  2006/03/07 03:32:22  scott
% fixing default|specified comps calculation -sm
%
% Revision 1.6  2006/03/07 02:30:54  scott
% worked on help msg; formatted, made accurate (only 2 files are now output, not 4).
% made the function accept the components argument (was always computing ersp/itc
% for ALL components, not just those asked for.  -sm
%
% Revision 1.5  2006/03/06 23:18:08  arno
% resave
%
% Revision 1.4  2006/03/03 23:39:05  arno
% put log
%

function [EEG_etc X times freqs] = std_ersp(EEG, comp, freqrange, timewindow, cycles, padratio, alpha, type, powbase)

EEG_etc = [];
if ~exist('type')
    type = 'ersp';
end

if isfield(EEG,'icaweights')
    numc = size(EEG.icaweights,1); % number of ICA comps
else
    error('EEG.icaweights not found');
end
if ~exist('comp') | isempty(comp)
   comps = 1:numc;
else
   comps = comp;
end

if exist('powbase') 
 if isempty(powbase) | isnan(powbase)
  powbase = NaN*ones(length(comps),1);  % default for timef()
 end
else
  powbase = NaN*ones(length(comps),1);  % default for timef()
end

if size(powbase,1) ~= length(comps)
   error('powbase should be of size (ncomps,nfreqs)');
end

% Check if ERSP information found in datasets and if fits requested parameters 
if isfield(EEG,'etc')
     if isfield(EEG.etc, [ 'ica' type ])
         if exist( fullfile(EEG.filepath, getfield(EEG.etc, [ 'ica' type ])))
             params = EEG.etc.icaerspparams;
             if sum(params.cycles ~= cycles)                   ...
                     | sum(params.freqrange ~= freqrange)   ...
                     | (padratio ~= params.padratio) ...
                     | (alpha~= params.alpha) ...
                 % if not as requested parameters, recompute ERSP/ITC
                 % i.e., continue
             else
                 return; % no need to compute ERSP
             end
         end;
     end
 end

% No ERSP/ITC information available
if isstr(EEG.data)
    TMP = eeg_checkset( EEG, 'loaddata' );  % load EEG.data and EEG.icaact
else
    TMP = EEG;
end
if isempty(TMP.icaact)                      % make icaact if necessary
    TMP.icaact = (TMP.icaweights*TMP.icasphere)* ...
                 reshape(TMP.data  , [ size(TMP.data,1)   size(TMP.data,2)*size(TMP.data,3) ]);
    TMP.icaact = reshape(TMP.icaact, [ size(TMP.icaact,1) size(TMP.data,2)*size(TMP.data,3) ]);
end;

% Compute ERSP parameters 
[time_range, winsize] = compute_ersp_times(cycles,  EEG.srate, ...
                                 [EEG.xmin EEG.xmax]*1000 , freqrange(1), padratio); 
if time_range(1) >= time_range(2)
    error(['std_ersp: parameters given for ' upper(type) ' calculation result in an invalid time range. Aborting. Please increase the lower frequency bound or change other parameters to resolve the problem. See >> timef details'] )
end

for k = 1:length(comps)  % for each (specified) component

    % Compute ERSP & ITC
    % [ersp,itc,powbase,times,freqs,erspboot,itcboot] = timef( TMP.icaact(comps(k), :) , ...
    %       EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, cycles ,'type', ...
    %          'phasecoher',  'plotersp', 'off', 'plotitc', 'off', 'powbase', powbase, ...
    %             'alpha',alpha,'padratio',padratio, 'plotphase','off','winsize',winsize);

    [ersp,itc,powbase,times,freqs,erspboot,itcboot] = timef( TMP.icaact(comps(k), :) , ...
          EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, cycles ,'type', ...
             'phasecoher',  'plotersp', 'off', 'plotitc', 'off', 'powbase', powbase(k,:), ...
                'padratio',padratio, 'plotphase','off','winsize',winsize);
   
    % Change frequency axis from linear scale to log scale (frequency values left in dB)

    [logfreqs,logersp] = logimagesc(times,freqs,ersp,'plot','off');  
    try
        logeboot(1,:) = interp1(log(freqs),erspboot(1,:),logfreqs','linear');
        logeboot(2,:) = interp1(log(freqs),erspboot(2,:),logfreqs','linear');
    catch
        logeboot = zeros(2, length(logfreqs));
    end;
    logbase = interp1(log(freqs),powbase,logfreqs','linear');

    [logfreqs,logitc] = logimagesc(times,freqs,itc,'plot','off'); 
    
    try
        logiboot = interp1(log(freqs),itcboot(1,:),logfreqs','linear');
    catch
        logiboot = zeros(1, length(logfreqs));
    end;

    if k == 1 % create data matrices
        all_ersp = zeros(length(freqs),(length(times)+3)*length(comps));
        all_itc = zeros(length(freqs),(length(times)+1)*length(comps)); % Save ITC info as well.
    end

    % store log-frequency spaced ersp and itc
    all_ersp(:,1+(k-1)*(length(times)+3):k*(length(times)+3) ) = [logersp logeboot' logbase'];
    all_itc(:,1+(k-1)*(length(times)+1):k*(length(times)+1) )  = [logitc  logiboot'];
end

% Save ERSP into float file
floatwrite(all_ersp, fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'icaersp']));
floatwrite(all_itc,  fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'icaitc']));

% Update the info in the dataset
EEG.etc.icaersp = [ EEG.filename(1:end-3) 'icaersp'];
EEG.etc.icaitc  = [ EEG.filename(1:end-3) 'icaitc'];

% Save ERSP parameters in the dataset
EEG.etc.icaerspparams.times      = times;
EEG.etc.icaerspparams.freqs      = freqs;
EEG.etc.icaerspparams.logfreqs   = logfreqs;
EEG.etc.icaerspparams.cycles     = cycles;
EEG.etc.icaerspparams.alpha      = alpha;
EEG.etc.icaerspparams.padratio   = padratio;
EEG.etc.icaerspparams.timewindow = timewindow;
EEG.etc.icaerspparams.freqrange  = freqrange;
EEG.etc.icaitcparams             = EEG.etc.icaerspparams;

try
    EEG.saved = 'no';
    EEG = pop_saveset( EEG, 'savemode', 'resave');
catch,
    error([ 'std_ersp: problem saving results into path ' EEG.filepath])
end
EEG_etc = EEG.etc; % return updated EEG.etc sub-structure
