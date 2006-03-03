% Usage:  
%   >> [EEG_etc] = cls_ersp(EEG, components, freqrange, timewindow,  cycles, padratio, alpha, overwrite, type)
%   If doesn't exist computes the ICA ERSPs / ITCs of a dataset. Updates the EEG structure in the 
%   Matlab environment and on the disk too!
%
% cls_ersp() - This function computes the ERSP & ITC information 
% of a dataset ICA components, saves that into a float files and 
% and saves pointers to the files in the EEG structure.
% If a float file with the information already exists, the function 
% loads the ERSP or ITC  information from it, unless the requested 
% parameters (input varaibles) were different, if so a window
% will pop-up to ask the user what to do.
% There is an option to specify only certain components, select a desired
% frequency range, timewindow, resolution, confidence level and the wavelet 
% parameters (number of cycles), see timef for details. 
% The type is which of the two (ERSP, ITC) to return, default is ersp. 
% The function returns the masked  (applied with the confidence level) log
% scaled (frequencies are in a log scale, so low frequencies are more
% pronounced than higher frequencies) ERSP or ITC, of the selected ICA components 
% in the requested frequency range and time window (the two are dependent).
% Four float files are saved, two for the ERSP and two for the ITC information.
% The two float files of ERSP / ITC are - one with the unprocessed ERSP /
% ITC and its separate significant level, and the other with the masked log
% scale ERSP / ITC.
% If the ERSPs / ITCs were already computed before, and a different set of
% parameters are desired, there is an overwrite variable to do so,
% which depending on the variable can open a pop-up window to ask the user
% what to do, use the existing information or overwrite it with new computation. 
% The frequency samples and time samples vectors are returned too,
% as well as the EEG sub-structure etc (i.e. EEG.etc), which is modified 
% with the pointers to the floating files and some information about them. 
%
%
% Inputs:
%   EEG            - an EEG data structure. 
%   components - [numeric vector] of the EEG structure for which an ERSP   
%                     and ITC will be computed.
%   freqrange - [minHz maxHz] the frequency range to compute the ERSP / ITC.
%   cycles      - If 0 -> Use FFTs (with constant window length) {0}
%                     If >0 -> Number of cycles in each analysis wavelet 
%                     If [wavecycles factor] -> wavelet cycles increase with frequency 
%                     beginning at wavecyles (0<factor<1; factor=1 -> no increase,
%                     standard wavelets; factor=0 -> fixed epoch length, as in FFT.
%   padratio  - FFT-length/winframes (2^k)                          
%                     Multiplies the number of output frequencies by
%                     dividing their spacing. When cycles==0, frequency
%                     spacing is (low_freq/padratio).
%   alpha        - If non-0, compute two-tailed bootstrap significance
%                     prob. level
%   overwrite  - [0|1|2] flag which indicates what to do with existing ERSP / ITC  information 
%                     0 - there is no ERSP / ITC  information saved in dataset 
%                     1- overwrite the saved ERSP / ITC of this dataset if exist
%                     2 - use exsiting ERSP / ITC dataset info, requested time range 
%                     and padratio might be smaller then saved in dataset. (default - 0).
%   type        - ['ersp'|'itc'] both ersp and itc are computed and saved, but only one 
%                     is returned (output X). default: 'ersp'. 
%
% Outputs:
%   EEG_etc    - the EEG dataset etc structure (i.e. EEG.etc), which is
%                      modified with the pointers and some information about
%                      the floating files that hold the dataset ERSP, logERSP, ITC  and logITC information.
%                      If the ERSP / ITC file already exists and wasn't overwritten this output will be empty. 
%   X              - the masked log ERSP / ITC of the requested ICA components in the selected 
%                     frequency and time range. 
%   times      - a time vector of the time points in which the ERSPs / ITCs were computed. 
%   freqs      - a frequency vector (log scaled) of the points in which the log ERSP / ITC was computed. 
%   overwrite - same as input option, only modified with what the user
%                     asked for the rest of the datasets in STUDY (from the
%                     pop-up menu, or from command line).
%
%  example: k = 1; [EEG_etc, X, overwrite] = cls_ersp(ALLEEG(k), [1:size(ALLEEG(k).icawinv,2)],[2 50], [], [3 0.5], 4, 0.01, 0);
%
%  See also  timef, cls_itc, cls_erp, cls_spec, cls_scalp, eeg_preclust, eeg_createdata         
%
% Authors:  Hilit Serby, SCCN, INC, UCSD, January, 2005

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

function [EEG_etc] = cls_ersp(EEG, comp, freqrange, timewindow, cycles, padratio, alpha, type)
EEG_etc = [];
if ~exist('type')
    type = 'ersp';
end

% Check if ersp information found in datasets and if fits requested parameters 
if isfield(EEG,'etc')
     if isfield(EEG.etc, [ 'ica' type])
         params = EEG.etc.icaerspparams;
         if sum(params.cycles ~= cycles) | sum(params.freqrange ~= freqrange)  | (padratio ~= params.padratio) | (alpha~= params.alpha) ...
                     % if not as requested parameters recompute ERSP / ITC
         else
             return; % no need to compute ERSP
         end
     end
 end

% No ERSP / ITC information available
if isstr(EEG.data)
    TMP = eeg_checkset( EEG, 'loaddata' ); %load EEG.data and EEG.icaact
else
    TMP = EEG;
end
if isempty(TMP.icaact)
    TMP.icaact = (TMP.icaweights*TMP.icasphere)* ...
                 reshape(TMP.data  , [ size(TMP.data,1)   size(TMP.data,2)*size(TMP.data,3) ]);
    TMP.icaact = reshape(TMP.icaact, [ size(TMP.icaact,1) size(TMP.data,2)*size(TMP.data,3) ]);
end;

%compute ERSP for all components
[time_range, winsize] = compute_ersp_times(cycles,  EEG.srate, [EEG.xmin EEG.xmax]*1000 , freqrange(1), padratio); 
if time_range(1) >= time_range(2)
    error(['cls_ersp: the parameters given for ' upper(type) ' calculation result in invalid time range. Aborting. Please change lower frequency bound or other parameters to resolve the problem.'] )
end

numc = size(EEG.icaweights,1); %number of ICA comp
for k = 1:numc
    % Compute ERSP & ITC
    [ersp,itc,powbase,times,freqs,erspboot,itcboot] = timef( TMP.icaact(k, :) , EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, cycles ,'type', ...
        'phasecoher',  'plotersp', 'off', 'plotitc', 'off', 'alpha',alpha,'padratio',padratio, 'plotphase','off','winsize',winsize);
   
    % Change frequency axis from linear scale to log scale (frequency values left in dB)
    [logfreqs,logersp] = logimagesc(times,freqs,ersp,'plot','off');  
    logeboot(1,:) = interp1(log(freqs),erspboot(1,:),logfreqs','linear');
    logeboot(2,:) = interp1(log(freqs),erspboot(2,:),logfreqs','linear');
    logbase = interp1(log(freqs),powbase,logfreqs','linear');
    [logfreqs,logitc] = logimagesc(times,freqs,itc,'plot','off'); 
    logiboot = interp1(log(freqs),itcboot(1,:),logfreqs','linear');

    if k == 1
        all_ersp = zeros(length(freqs),(length(times)+3)*numc);
        all_itc = zeros(length(freqs),(length(times)+1)*numc); %save itc info as well.
    end
    all_ersp(:,1+(k-1)*(length(times)+3):k*(length(times)+3) ) = [logersp logeboot' logbase'];
    all_itc(:,1+(k-1)*(length(times)+1):k*(length(times)+1) ) = [logitc logiboot'];
end

%save ERSP in file
floatwrite(all_ersp, fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'icaersp']));
floatwrite(all_itc,  fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'icaitc']));
%update the info in the dataset
EEG.etc.icaersp = [ EEG.filename(1:end-3) 'icaersp'];
EEG.etc.icaitc  = [ EEG.filename(1:end-3) 'icaitc'];
%save ersp parameters in the dataset
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
    EEG = pop_saveset( EEG, 'savemode', 'resave');
catch,
    error([ 'cls_ersp: problems saving into path ' EEG.filepath])
end
EEG_etc = EEG.etc;
