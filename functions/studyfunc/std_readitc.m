%
% std_readitc() - returns the log-frequency inter-trial coherence (ITC) for a 
%                 specified ICA component. The component ITCs for the dataset 
%                 are assumed to have been saved in a Matlab file, 
%                 [dataset_)name].icaitc, in the same directory as the dataset.
%                 If no such file exists, use std_ersp() to create it, else 
%                 the pre-clustering functions that call it: pop_preclust, 
%                 std_preclust().  The input variables used to compute the
%                 ITC are returnsd: frequency_range, time_range, resolution, 
%                 probability_threshold, and wavelet type (FFT | wavelet 
%                 cycles). See timef() for details. 
% Usage:    
%   >> [logersp, logfreqs, times] = std_readitc(ALLEEG, setindx, component, ...
%                                                       time_range, freq_range);  
% Inputs:
%   ALLEEG     - EEG dataset vector (can also be an EEG set). 
%                Must contain the dataset of interest (see 'setindx' below).
%   setindx    -  [integer] index of the EEG dataset in ALLEEG for which 
%                 to return the log-frequency ITC.
%   component  - [integer] component index in the selected EEG dataset for 
%                which to return the ITC. 
%   time_range - [min max in ms] ITC time window 
%   freq_range - [min max in Hz] ITC frequency range 
%
% Outputs:
%   logitc     - the equal log-spaced frequency ITC for the requested ICA 
%                component in the specified dataset. Its dimensions are 
%                (equal log-spaced) frequencies by times. 
%   logfreqs   - vector of ITC (log-spaced) frequencies, in Hz
%   times      - vector of ITC times (latencies), in ms.
%   params     - full structure of ITC parameters saved
%
%  See also  std_ersp(), std_readersp(), pop_preclust(), eeg_preclust(), 
%               eeg_createdata()
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
% Revision 1.11  2006/05/03 17:55:54  arno
% allow to read data channels
%
% Revision 1.10  2006/03/14 02:24:41  scott
% help msg
%
% Revision 1.9  2006/03/14 02:11:46  scott
% help msg
%
% Revision 1.8  2006/03/11 07:28:50  arno
% header info
%
% Revision 1.7  2006/03/10 15:49:07  arno
% fix reading ITC
%
% Revision 1.6  2006/03/10 00:39:37  arno
% error msg
%
% Revision 1.5  2006/03/10 00:01:52  arno
% remove old ITC code
%
% Revision 1.4  2006/03/09 23:40:10  arno
% read new ITC format
%
% Revision 1.3  2006/03/08 20:35:03  arno
% rename func
%
% Revision 1.2  2006/03/07 22:17:34  arno
% use fullfile
%

function [logitc, logfreqs, timevals, params] = std_readitc(ALLEEG, abset, comp, timewindow, freqrange);

if nargin < 4
    timewindow = [];
end;
if nargin < 5
    freqrange = [];
end;

% multiple entry
% --------------
if length(comp) > 1
    for index = 1:length(comp)
        [tmpitc, logfreqs, timevals, params] = std_readitc(ALLEEG, abset, comp(index), timewindow, freqrange);
        logitc(index,:,:,:) = tmpitc;
    end;
    return;
end;

for k = 1: length(abset)    
    
    if comp < 0
        filename = fullfile( ALLEEG(abset(k)).filepath,[ ALLEEG(abset(k)).filename(1:end-3) 'datitc']);
        comp   = -comp;
        prefix = 'chan';
    else    
        filename = fullfile( ALLEEG(abset(k)).filepath,[ ALLEEG(abset(k)).filename(1:end-3) 'icaitc']);
        prefix = 'comp';
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
    if isempty(comp)
        logitc    = [];
        logfreqs  = [];
        timevals  = [];
        return;
    end;
    tmpitc   = load( '-mat', filename, 'parameters', 'times', 'freqs', ...
                     [ prefix int2str(comp) '_itc'], ...
                     [ prefix int2str(comp) '_itcboot']);
    
    tlen      = length(tmpitc.times);
    flen      = length(tmpitc.freqs);
    itcall{k}     = double(getfield(tmpitc, [ prefix int2str(comp) '_itc']));
    itcallboot{k} = double(getfield(tmpitc, [ prefix int2str(comp) '_itcboot']));

end

% select plotting or clustering time/freq range
% ---------------------------------------------
if ~isempty(timewindow)
    if timewindow(1) > tmpitc.times(1) | timewindow(end) < tmpitc.times(end)
        maxind = max(find(tmpitc.times <= timewindow(end)));
        minind = min(find(tmpitc.times >= timewindow(1)));
    else
        minind = 1;
        maxind = tlen;
    end
else
    minind = 1;
    maxind = tlen;
end
if ~isempty(freqrange)
    if freqrange(1) > exp(1)^tmpitc.freqs(1) | freqrange(end) < exp(1)^tmpitc.freqs(end)
        fmaxind = max(find(tmpitc.freqs <= freqrange(end)));
        fminind = min(find(tmpitc.freqs >= freqrange(1)));
    else
        fminind = 1;
        fmaxind = flen;
    end
else
    fminind = 1;
    fmaxind = flen;
end

% Mask ITC
% ---------
if ~isempty(itcallboot{1})
    for cond  = 1:length(abset)
        %maxitc= repmat(itcallboot{cond}',1,size(itcall{1},2));
        %itcall{cond}(find(itcall{cond}<maxitc)) = 0;
    end
end;

% return parameters
% ----------------
for cond  = 1:length(abset)
    itc = itcall{cond}(fminind:fmaxind,minind:maxind);
    logitc(:,:,cond) = itc;
end;
logfreqs = tmpitc.freqs(fminind:fmaxind);
timevals = tmpitc.times(minind:maxind);

% remove duplicates in the list of parameters
% -------------------------------------------
function cella = removedup(cella)
    [tmp indices] = unique(cella(1:2:end));
    if length(tmp) ~= length(cella)/2
        fprintf('Warning: duplicate ''key'', ''val'' parameter(s), keeping the last one(s)\n');
    end;
    cella = cella(sort(union(indices*2-1, indices*2)));
