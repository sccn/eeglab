%
% cls_ersp() - Computes ERSP and/or ITC information for ICA components of a dataset, 
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
%              If the ERSPs/ITCs were previously saved in these files and a new set of
%              ERSP/ITC parameters are used, there is an overwrite flag to allow this.
%              Depending on the flag, cls_ersp can open a pop-up window to ask the user
%              what to do, use the existing information, or overwrite it. Vectors of
%              frequencies and latencies for the ERSP/ITC images are returned separately, 
%              as well as the EEG.etc sub-structure modified with pointers to the float 
%              files and some information about them. 
% Usage:  
%              >> [EEG_etc X times freqs overwrite] = cls_ersp(EEG, components,        ...
%                                                              freqrange, timewindow,  ...
%                                                              cycles, padratio, alpha,...
%                                                              overwrite, type);
%
%              % If they do not exist, this computes the ICA component ERSPs/ITCs 
%              % for a dataset. Saves the computed images in dataset name files 
%              % with extensions .icaersp and .icaitc
%              % Also updates the EEG structure in the Matlab environment 
%              % and saves the modified dataset to disk.
% Inputs:
%
%   EEG        - an EEG data structure. 
%   components - [numeric vector] of the EEG structure for which an ERSP   
%                 and ITC will be computed {default: compute for all components}
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
%   overwrite  - [0|1|2] flag to indicate what to do with pre-existing ERSP/ITC info:
%                 0 - there is no ERSP/ITC  information saved in the dataset 
%                 1 - overwrite the saved ERSP/ITC of this dataset if exist
%                 2 - use exsiting ERSP/ITC dataset info, requested time range 
%                     and padratio (might be smaller then saved in dataset). 
%                 {default: 0}
%   type       - ['ersp'|'itc'] both ERSP and ITC images are computed and saved to disk, 
%                 but only one is returned (output X) {default: 'ersp'}
% Outputs:
%
%   EEG_etc   - the EEG dataset .etc sub-structure (i.e. EEG.etc), which is
%               modified with the pointers and some information about
%               the floating files that hold the dataset ERSP, logERSP, ITC  
%               and logITC information. If the ERSP/ITC file already exists 
%               and wasn't overwritten this output will be empty. 
%   X         - the masked log ERSP/ITC of the requested ICA components in the selected 
%               frequency and time range. 
%   times     - a time vector of the time points in which the ERSPs/ITCs were computed. 
%   freqs     - a frequency vector (log scaled) of the points in which the 
%               log ERSP/ITC was computed. 
%   overwrite - same as input option, only modified with what the user
%               asked for the rest of the datasets in STUDY (from the
%               pop-up menu, or from command line).
%
% Files written or modified:     [dataset_filename].icaersp   <-- saved component ERSPs
%                                [dataset_filename].icaitc    <-- saved component ITCs
%                                [dataset_filename].set       <-- re-saved dataset
%  Example: 
%            >> k = 1; 
%            >> [EEG_etc, X, overwrite] = cls_ersp(ALLEEG(k), [1:size(ALLEEG(k).icawinv,2)],...
%                                                     [2 50], [], [3 0.5], 4, 0.01, 0);
%
%  See also: timef(), cls_itc(), cls_erp(), cls_spec(), cls_scalp(), eeg_preclust(), 
%            eeg_createdata()
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
% Revision 1.5  2006/03/06 23:18:08  arno
% resave
%
% Revision 1.4  2006/03/03 23:39:05  arno
% put log
%

function [EEG_etc] = cls_ersp(EEG, comp, freqrange, timewindow, cycles, padratio, alpha, type)

EEG_etc = [];
if ~exist('type')
    type = 'ersp';
end
if isempty(comp)
   numc = size(EEG.icaweights,1); % number of ICA comps
   comp = 1:numc;
end

% Check if ERSP information found in datasets and if fits requested parameters 
if isfield(EEG,'etc')
     if isfield(EEG.etc, [ 'ica' type])
         params = EEG.etc.icaerspparams;
         if sum(params.cycles ~= cycles)                   ...
                    | sum(params.freqrange ~= freqrange)   ...
                           | (padratio ~= params.padratio) ...
                                  | (alpha~= params.alpha) ...
             % if not as requested parameters, recompute ERSP/ITC
         else
             return; % no need to compute ERSP
         end
     end
 end

% No ERSP/ITC information available
if isstr(EEG.data)
    TMP = eeg_checkset( EEG, 'loaddata' );  % load EEG.data and EEG.icaact
else
    TMP = EEG;
end
if isempty(TMP.icaact)
    TMP.icaact = (TMP.icaweights*TMP.icasphere)* ...
                 reshape(TMP.data  , [ size(TMP.data,1)   size(TMP.data,2)*size(TMP.data,3) ]);
    TMP.icaact = reshape(TMP.icaact, [ size(TMP.icaact,1) size(TMP.data,2)*size(TMP.data,3) ]);
end;

% Compute ERSP for all components
[time_range, winsize] = compute_ersp_times(cycles,  EEG.srate, ...
                                 [EEG.xmin EEG.xmax]*1000 , freqrange(1), padratio); 
if time_range(1) >= time_range(2)
    error(['cls_ersp: parameters given for ' upper(type) ' calculation result in an invalid time range. Aborting. Please change the lower frequency bound or other parameters to resolve the problem.'] )
end

for k = comps  % for each (specified) component

    % Compute ERSP & ITC
    [ersp,itc,powbase,times,freqs,erspboot,itcboot] = timef( TMP.icaact(k, :) , ...
          EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, cycles ,'type', ...
             'phasecoher',  'plotersp', 'off', 'plotitc', 'off', ...
                'alpha',alpha,'padratio',padratio, 'plotphase','off','winsize',winsize);
   
    % Change frequency axis from linear scale to log scale (frequency values left in dB)

    [logfreqs,logersp] = logimagesc(times,freqs,ersp,'plot','off');  
    logeboot(1,:) = interp1(log(freqs),erspboot(1,:),logfreqs','linear');
    logeboot(2,:) = interp1(log(freqs),erspboot(2,:),logfreqs','linear');
    logbase = interp1(log(freqs),powbase,logfreqs','linear');

    [logfreqs,logitc] = logimagesc(times,freqs,itc,'plot','off'); 
    logiboot = interp1(log(freqs),itcboot(1,:),logfreqs','linear');

    if k == comps(1)
        all_ersp = zeros(length(freqs),(length(times)+3)*numc);
        all_itc = zeros(length(freqs),(length(times)+1)*numc); % Save ITC info as well.
    end
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
    error([ 'cls_ersp: problem saving results into path ' EEG.filepath])
end
EEG_etc = EEG.etc; % return updated EEG.etc sub-structure
