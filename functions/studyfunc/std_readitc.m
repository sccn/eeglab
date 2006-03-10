% std_readitc() - Given the ALLEEG structure, a specific EEG dataset index, 
% and a specific component, the function returns the (frequency) log scaled 
% ITC of that ICA component. 
% The ITC of the dataset ICA components is assumed to be saved in a float 
% file, the EEG dataset include a pointer to this file. If such a float file doesn't exist,
% you can use the std_ersp() function to create it, or use the pre - clustering functions
% that call it: pop_preclust, eeg_preclust & eeg_createdata. The ITC information  
% is specific to the input variables used to compute it, the frequency range, 
% timewindow, resolution, confidence level and the wavelet type (FFT / wavelet
%  cycles), see timef for details. 
% Along with the ITC of the selected ICA component the function returns  
% the log frequency vector of the ITC samples. 
%
% Usage:    
%   >> [logersp, logfreqs] = std_readitc(ALLEEG, abset, component);  
%   This functions returns the log scaled ITC of an ICA component. 
%   The information is loaded from a float file, which a pointer 
%   to is saved in the EEG dataset. The float file was created
%   by the pre - clustering function std_ersp. 
%
% Inputs:
%   ALLEEG     - an EEGLAB data structure, which holds EEG sets (can also be one EEG set). 
%                      ALLEEG must contain the dataset of interest (the setind).
%   setind         -  [integer] an index of an EEG dataset in the ALLEEG
%                      structure, for which to get the component log scaled ITC.
%   component - [integer] a component index in the selected EEG dataset for which 
%                      an ITC will be returned. 
%
% Outputs:
%   logitc       - the log scaled ITC of the requested ICA component in the
%                      selected dataset. This is frequencies (log scaled)
%                      by time, time - frequency decomposition inter-trial coherence 
%                      of the ICA activations.
%   logfreqs   - a vector of the frequencies points in log scale. 
%
%  See also  std_ersp, std_readersp, pop_preclust, eeg_preclust, eeg_createdata           
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

for k = 1: length(abset)    
    
    filename = fullfile( ALLEEG(abset(k)).filepath,[ ALLEEG(abset(k)).filename(1:end-3) 'icaitc']);
    try
        tmpersp   = load( '-mat', filename, 'parameters', 'times', 'freqs');
    catch
        error( [ 'Cannot read file ''' filename '''' ]);
    end;
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
                     [ 'comp' int2str(comp) '_itc'], ...
                     [ 'comp' int2str(comp) '_itcboot']);
    
    tlen      = length(tmpitc.times);
    flen      = length(tmpitc.freqs);
    itcall{k}     = double(getfield(tmpitc, [ 'comp' int2str(comp) '_itc']));
    itcallboot{k} = double(getfield(tmpitc, [ 'comp' int2str(comp) '_itcboot']));

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
