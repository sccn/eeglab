% std_readitc()  - load ITC measures for data channels or 
%                  for all components of a specified cluster.
% Usage:
%         >> [STUDY, itcdata, times, freqs] = ...
%                   std_readitc(STUDY, ALLEEG, varargin);
% Inputs:
%       STUDY - studyset structure containing some or all files in ALLEEG
%      ALLEEG - vector of loaded EEG datasets
%
% Optional inputs:
%  'channels'  - [cell] list of channels to import {default: all}
%  'clusters'  - [integer] list of clusters to import {[]|default: all but
%                the parent cluster (1) and any 'NotClust' clusters}
%  'freqrange' - [min max] frequency range {default: whole measure range}
%  'subject'    - [string] select a specific subject {default:all}
%  'component'  - [integer] select a specific component in a cluster
%                 {default:all}
%  'singletrials' - ['on'|'off'] load single trials data (if available)
%
% Output:
%  STUDY    - updated studyset structure
%  itcdata  - [cell array] ITC data (the cell array size is 
%             condition x groups)
%  times    - [float array] array of time points
%  freqs    - [float array] array of frequencies
%
% Author: Arnaud Delorme, CERCO, 2006-

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

function [STUDY, erspdata, alltimes, allfreqs, erspbase] = std_readersp(STUDY, ALLEEG, varargin);

[STUDY, erspdata, alltimes, allfreqs] = std_readersp(STUDY, ALLEEG, 'infotype','itc', varargin{:});
return;

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
        %fprintf('Warning: duplicate ''key'', ''val'' parameter(s), keeping the last one(s)\n');
    end;
    cella = cella(sort(union(indices*2-1, indices*2)));
