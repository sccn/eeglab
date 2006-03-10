% std_readersp() - Given the ALLEEG structure, a specific EEG dataset index, 
%                  and a specific component, return the mean log-frequency scaled 
%                  ERSP for a specified ICA component. The ERSP of the dataset ICA 
%                  component is assumed to be saved in a float file. The EEG dataset 
%                  includes a pointer to this file. If such a float file doesn't 
%                  exist, you can use the std_ersp() function to create it, or use 
%                  pop_preclust(), eeg_preclust(), and eeg_createdata(). The ERSP 
%                  information is specific to the input variables used to compute it: 
%                  frequency range, timewindow, resolution, probability threshold, 
%                  and wavelet type (FFT|wavelet_cycles), see >> timef details 
%                  Along with the ERSP of the selected ICA component the function 
%                  returns the log frequency vector of the ERSP samples. 
% Usage:    
%       >> [logersp, logfreqs] = std_readersp(ALLEEG, setindex, component);  
%
% Inputs:
%   ALLEEG     - an EEGLAB data structure, which holds EEG sets (can also be one EEG
%                dataset). ALLEEG must contain the dataset of interest (the setind).
%   setindex   - [integer] an index of an EEG dataset in the ALLEEG structure, for 
%                which to get the component log scaled ERSP.
%   component  - [integer] a component index in the selected EEG dataset for which 
%                an ERSP will be returned. 
% Outputs:
%   logersp    - the log scaled ERSP of the requested ICA component in the
%                      selected dataset. This is frequencies (log scaled)
%                      by time, time - frequency decomposition of the ICA
%                      activations.
%   logfreqs   - a vector of the frequencies points in log scale. 
%
%  See also  std_ersp, pop_preclust, eeg_preclust, eeg_createdata           
%
% Authors: Arnaud Delorme, Hilit Serby, SCCN, INC, UCSD, February, 2005

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

function [logersp, logfreqs, timevals, params] = std_readersp(ALLEEG, abset, comp, timewindow, freqrange);

for k = 1: length(abset)    
    
    filename = fullfile( ALLEEG(abset(k)).filepath,[ ALLEEG(abset(k)).filename(1:end-3) 'icaersp']);
    try
        tmpersp   = load( '-mat', filename, 'parameters', 'times', 'freqs');
    catch
        error( [ 'Cannot read file ''' filename '''' ]);
    end;
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
    tmp2ersp   = load( '-mat', filename, ...
                     [ 'comp' int2str(comp) '_ersp'], ...
                     [ 'comp' int2str(comp) '_erspbase'], ...
                     [ 'comp' int2str(comp) '_erspboot']);
    
    erspall{k}     = double(getfield(tmpersp2, [ 'comp' int2str(comp) '_ersp']));
    erspallboot{k} = double(getfield(tmpersp2, [ 'comp' int2str(comp) '_erspboot']));
    erspallbase{k} = double(getfield(tmpersp2, [ 'comp' int2str(comp) '_erspbase']));

end

% compute average baseline across conditions
% ------------------------------------------
if length(abset) > 1 
    % mean baseline for requested component across conditions
    % -------------------------------------------------------
	ave_baseline = 0; 
	for cond = 1:length(abset)
        ave_baseline = ave_baseline + erspallbase{cond}/length(abset);
	end    
    
    % apply mean baseline
    % -------------------
    for cond = 1:length(abset)
        % add back former baseline to the ERSP and subtract new baseline
        erspall{cond} = erspall{cond} + 10*log10(repmat( erspallbase{cond}',[1 tlen]));
        erspall{cond} = erspall{cond} - 10*log10(repmat( ave_baseline'     ,[1 tlen]));
        
        % same for bootstrap array
        if ~isempty(erspallboot{cond})
            erspallboot{cond} = erspallboot{cond} + 10*log10(repmat(erspallbase{cond}',[1 2]))';
            erspallboot{cond} = erspallboot{cond} - 10*log10(repmat(ave_baseline'     ,[1 2]))';  
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
if ~isempty(erspallboot{1})
    for cond  = 1:length(abset)
        minersp= repmat(erspallboot{cond}(1,:)',1,size(erspall{1},2));
        maxersp= repmat(erspallboot{cond}(2,:)',1,size(erspall{1},2));
        erspall{cond}(find(erspall{cond}<maxersp & erspall{cond}>minersp)) = 0;
    end
end;

% return parameters
% ----------------
for cond  = 1:length(abset)
    ersp = erspall{cond}(fminind:fmaxind,minind:maxind);
    logersp(:,:,cond) = ersp;
end;
logfreqs = tmpersp.freqs(fminind:fmaxind);
timevals = tmpersp.times(minind:maxind);

