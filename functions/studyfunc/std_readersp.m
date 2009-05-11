% std_readersp() - Returns the equal-log-frequency spaced mean event-related spectral 
%                  perturbation(s) (ERSP(s)) for a requested ICA component. The ERSP 
%                  is assumed to have been saved in a Matlab file, 
%                          [dataset_name].icaersp 
%                  in the same folder as the dataset file.
%                  If this file does not exist, use std_ersp() to create it, else 
%                  use a pre-clustering function: pop_preclust() or std_preclust(), 
%                  that calls it. Interpretation of the ERSP requires some input 
%                  variables used to compute it: frequency range, window width, 
%                  resolution, probability threshold, and wavelet type (FFT or 
%                  wavelet_cycles). See >> timef help and >> timef details
% Usage:    
%     >>  [logersp, logfreqs, timevals, params, baseersp] = std_readersp(ALLEEG, setindex, component, ...
%                                                            time_range, freq_range);  
% Inputs:
%   ALLEEG     - vector of EEG datasets (can also be one EEG dataset). Must contain 
%                the dataset of interest (see 'setind' below).
%   setindex   - [integer] index of the EEG dataset in the ALLEEG structure for 
%                which to read a component log-frequency ERSP.
%   component  - [integer] component index in the selected EEG dataset for which 
%                to read the ERSP 
%   time_range - [min max ms] ERSP time (latency) range of interest
%   freq_range - [min max Hz] ERSP frequency range of interest
%
% Outputs:
%   logersp    - the log-frequency ERSP for the requested ICA component 
%                in the specified dataset. Dimensions: (equal log-spaced) 
%                frequencies by epoch latencies (unit: dB diff from baseline)
%   logfreqs   - vector of equal-log-spaced ERSP frequencies (Hz) 
%   timevals      - vector of ERSP times (latencies) (s)
%   params     - structure of timef() parameters saved with the ERSP
%   baseersp   - condition-specific baseline. Third dimesnsion corresponds to 
%                conditions.
%   
% See also:  std_ersp(), pop_preclust(), std_preclust(), timef()
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, February, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, October 11, 2004, arno@sccn.ucsd.edu
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
% Revision 1.32  2007/10/25 21:18:28  nima
% output parameters updated.
%
% Revision 1.31  2007/08/06 19:32:18  arno
% fix last changes
%
% Revision 1.30  2007/08/06 19:29:30  arno
% multiple condition baseline now in log space
%
% Revision 1.29  2007/05/19 01:15:09  toby
% error message incorrectly claiming precomputed data not found corrected
%
% Revision 1.28  2007/05/18 16:23:29  peter
% typo
%
% Revision 1.27  2007/05/18 16:17:54  peter
% typo
%
% Revision 1.26  2007/05/18 16:15:46  peter
% typo
%
% Revision 1.25  2007/05/18 16:07:09  peter
% added test for reading ERSP component vars.
%
% Revision 1.24  2007/01/26 18:06:35  arno
% nothing
%
% Revision 1.23  2006/11/10 02:28:46  arno
% bootstrap concatenation
%
% Revision 1.22  2006/11/03 19:45:06  arno
% do not mask ERSP by default
%
% Revision 1.21  2006/09/12 18:56:48  arno
% channel compatibility and many more features
%
% Revision 1.20  2006/05/13 17:46:22  arno
% transposing baseline to prevent crash
%
% Revision 1.19  2006/05/03 18:20:40  arno
% allowing to read data channels
%
% Revision 1.18  2006/03/29 17:47:10  scott
% help msg
%
% Revision 1.17  2006/03/28 14:58:10  arno
% reading ersp channel
%
% Revision 1.16  2006/03/22 00:46:32  scott
% help msg format only
%
% Revision 1.15  2006/03/14 03:02:51  scott
% help msg
%
% Revision 1.14  2006/03/14 02:23:18  scott
% help msg
%
% Revision 1.13  2006/03/14 01:59:38  scott
% help msg
%
% Revision 1.12  2006/03/13 19:09:03  arno
% no time range and freq range
%
% Revision 1.11  2006/03/11 07:28:49  arno
% header info
%
% Revision 1.10  2006/03/11 07:21:08  arno
% header
%
% Revision 1.9  2006/03/10 17:44:25  arno
% typo
%
% Revision 1.8  2006/03/10 15:49:19  arno
% fix reading ERSP
%
% Revision 1.7  2006/03/10 03:25:32  scott
% help msg -- ARNO, please check  -sm
%
% Revision 1.6  2006/03/10 00:39:39  arno
% error msg
%
% Revision 1.5  2006/03/09 23:29:34  arno
% implement new ERSP from Matlab and different structure ec...
%
% Revision 1.4  2006/03/09 19:37:06  arno
% header
%
% Revision 1.3  2006/03/08 20:34:24  arno
% rename func
%
% Revision 1.2  2006/03/07 22:16:23  arno
% use fullfile
%

function [logersp, logfreqs, timevals, params, baseersp] = std_readersp(ALLEEG, abset, comp, timewindow, freqrange);

if nargin < 4
    timewindow = [];
end;
if nargin < 5
    freqrange = [];
end;

% multiple entry
% --------------
if (length(comp) > 1 & comp(1) > 0) | length(comp) > length(abset) % recursive call if multiple components
    for index = 1:length(comp)
        [tmpersp, logfreqs, timevals, params, tmpbase] = std_readersp(ALLEEG, abset, comp(index), timewindow, freqrange);
        logersp(index,:,:,:) = tmpersp;
        baseersp(index,:,:)  = tmpbase;
    end;
    return;
end;

for k = 1: length(abset)    
    
    if comp < 0
        filename = fullfile( ALLEEG(abset(k)).filepath,[ ALLEEG(abset(k)).filename(1:end-3) 'datersp']);
        comptmp  = -comp(k);
        prefix   = 'chan';
    else    
        filename = fullfile( ALLEEG(abset(k)).filepath,[ ALLEEG(abset(k)).filename(1:end-3) 'icaersp']);
        prefix   = 'comp';
        comptmp  = comp;
    end;
    try
        tmpersp   = load( '-mat', filename, 'parameters', 'times', 'freqs');
    catch
        error( [ 'Cannot read file ''' filename '''' ]);
    end;
    
    tmpersp.parameters = removedup(tmpersp.parameters);
    params    = struct(tmpersp.parameters{:});
    params.times = tmpersp.times;
    params.freqs = tmpersp.freqs;
    tlen         = length(tmpersp.times);
    flen         = length(tmpersp.freqs);
    if isempty(comp)
        logersp   = [];
        logfreqs  = [];
        timevals  = [];
        return;
    end;
    if isinf(comp) % only read time and freqs
        abset = [];
        erspallboot{1} = [];
        logersp = [];
    else
        tmpersp2   = load( '-mat', filename, ...
                         [ prefix int2str(comptmp) '_ersp'], ...
                         [ prefix int2str(comptmp) '_erspbase'], ...
                         [ prefix int2str(comptmp) '_erspboot']);
        if ~isfield(tmpersp2,[ prefix int2str(comptmp) '_ersp']) |   ...
           ~isfield(tmpersp2,[ prefix int2str(comptmp) '_erspbase']) |   ...
           ~isfield(tmpersp2,[ prefix int2str(comptmp) '_erspboot']) 
          fprintf('\nERSP data for component %d not found!\n       In file %s\n',comptmp,filename); 
        end
        erspall{k}     = double(getfield(tmpersp2, [ prefix int2str(comptmp) '_ersp']));
        erspallboot{k} = double(getfield(tmpersp2, [ prefix int2str(comptmp) '_erspboot']));
        erspallbase{k} = double(getfield(tmpersp2, [ prefix int2str(comptmp) '_erspbase']));
    end;
end

% compute average baseline across conditions
% ------------------------------------------
if length(abset) > 1 
    % mean baseline for requested component across conditions
    % -------------------------------------------------------
	ave_baseline = zeros(size(erspallbase{1})); 
	for cond = 1:length(abset)
        ave_baseline = ave_baseline + erspallbase{cond}/length(abset);
	end    
    
    % apply mean baseline
    % -------------------
    for cond = 1:length(abset)
        % add back former baseline to the ERSP and subtract new baseline
        if any(erspallbase{cond}) > 100, 
            error([ 'You must recompute time-frequency decomposition' 10 ...
                    'changes greater than 100dB detected in ERSP which means' 10 ...
                    'that you have most likely computed these ERSP using a' 10 ...
                    'version of EEGLAB < 6.00' ]);
        end;
        
        erspall{cond} = erspall{cond} + repmat( erspallbase{cond}(:),[1 tlen]);
        erspall{cond} = erspall{cond} - repmat( ave_baseline(:)     ,[1 tlen]);
        
        % same for bootstrap array
        if ~isempty(erspallboot{cond})
            erspallboot{cond} = erspallboot{cond} + repmat(erspallbase{cond}',[1 2]);
            erspallboot{cond} = erspallboot{cond} - repmat(ave_baseline'     ,[1 2]);  
        end;
    end
end

% select plotting or clustering time/freq range
% ---------------------------------------------
if ~isempty(timewindow)
    if timewindow(1) > tmpersp.times(1) | timewindow(end) < tmpersp.times(end)
        maxind = max(find(tmpersp.times <= timewindow(end)));
        minind = min(find(tmpersp.times >= timewindow(1)));
    else
        minind = 1;
        maxind = tlen;
    end
else
    minind = 1;
    maxind = tlen;
end
if ~isempty(freqrange)
    if freqrange(1) > exp(1)^tmpersp.freqs(1) | freqrange(end) < exp(1)^tmpersp.freqs(end)
        fmaxind = max(find(tmpersp.freqs <= freqrange(end)));
        fminind = min(find(tmpersp.freqs >= freqrange(1)));
    else
        fminind = 1;
        fmaxind = flen;
    end
else
    fminind = 1;
    fmaxind = flen;
end

% Mask ERSP
% ---------
%if ~isempty(erspallboot{1})
%    for cond  = 1:length(abset)
%        minersp= repmat(erspallboot{cond}(1,:)',1,size(erspall{1},2));
%        maxersp= repmat(erspallboot{cond}(2,:)',1,size(erspall{1},2));
%        erspall{cond}(find(erspall{cond}<maxersp & erspall{cond}>minersp)) = 0;
%    end
%end;

% return parameters
% ----------------
for cond  = 1:length(abset)
    ersp = erspall{cond}(fminind:fmaxind,minind:maxind);
    logersp(:,:,cond) = ersp;
    baseersp(:,cond)  = erspallbase{cond}(fminind:fmaxind)';
end;
logfreqs = tmpersp.freqs(fminind:fmaxind);
timevals = tmpersp.times(minind:maxind);

% remove duplicates in the list of parameters
% -------------------------------------------
function cella = removedup(cella)
    [tmp indices] = unique(cella(1:2:end));
    if length(tmp) ~= length(cella)/2
        %fprintf('Warning: duplicate ''key'', ''val'' parameter(s), keeping the last one(s)\n');
    end;
    cella = cella(sort(union(indices*2-1, indices*2)));
