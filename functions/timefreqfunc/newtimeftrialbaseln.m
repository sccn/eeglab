% newtimeftrialbaseln() - Remove baseline power values for single trials in 
%        newtimef. This only remove single trial baseline (not the average
%        baseline). 
%
% Usage:
%   >>  tf = newtimefbaseln(tf, tvals, 'key', val);
%
% Inputs:
%   tf       - [3-D or 4-D array] single-trial spectral estimates
%           [freqs x times x trials] or [channels x freqs x times x trials)
%           The function may also process cell arrays of such inputs.
%   tvals    - [array] time values
%
% Optional inputs: 'baseline', 'basenorm' and 'trialbase'. Same definition 
%  as for newtimef. If trialbase is 'off' this function does nothing.
%
% Outputs:
%   tf        - Baseline correct power (same size as input)
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, August 2016

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2016, arno@sccn.ucsd.edu
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

function PP = newtimeftrialbaseln(PPori, timesout, varargin)

if nargin < 3
    help newtimefbaseln;
    return;
end;

[g timefreqopts ] = finputcheck(varargin, ...
    {'basenorm'      'string'    {'on','off'} 'off';
    'baseline'      'real'      []          0;
    'trialbase'     'string'    {'on','off','full'} 'off';
    'verbose'       'string'    {'on','off'} 'on';
    }, 'newtimeftrialbaseln', 'ignore');
if ischar(g) error(g); return; end;
PP = PPori; if ~iscell(PP), PP = { PP }; end;

% ---------------
% baseline length
% ---------------
if size(g.baseline,2) == 2
    baseln = [];
    for index = 1:size(g.baseline,1)
        tmptime   = find(timesout >= g.baseline(index,1) & timesout <= g.baseline(index,2));
        baseln = union_bc(baseln, tmptime);
    end;
    if length(baseln)==0
        error( [ 'There are no sample points found in the default baseline.' 10 ...
            'This may happen even though data time limits overlap with' 10 ...
            'the baseline period (because of the time-freq. window width).' 10 ...
            'Either disable the baseline, change the baseline limits.' ] );
    end
else
    if ~isempty(find(timesout < g.baseline))
         baseln = find(timesout < g.baseline); % subtract means of pre-0 (centered) windows
    else baseln = 1:length(timesout); % use all times as baseline
    end
end;

for ind = 1:length(PP(:))
    
    P = PP{ind};
    
    % -----------------------------------------
    % remove baseline on a trial by trial basis
    % -----------------------------------------
    if strcmpi(g.trialbase, 'on'), tmpbase = baseln;
    else                           tmpbase = 1:size(P,2); % full baseline
    end;
    if ~strcmpi(g.trialbase, 'off')
        if ndims(P) == 4
            mbase = mean(P(:,:,tmpbase,:),3);
            if strcmpi(g.basenorm, 'on')
                mstd = std(P(:,:,tmpbase,:),[],3);
                P = bsxfun(@rdivide, bsxfun(@minus, P, mbase), mstd);
            else P = bsxfun(@rdivide, P, mbase);
            end;
        else
            mbase = mean(P(:,tmpbase,:),2);
            if strcmpi(g.basenorm, 'on')
                mstd = std(P(:,tmpbase,:),[],2);
                P = (P-repmat(mbase,[1 size(P,2) 1]))./repmat(mstd,[1 size(P,2) 1]); % convert to log then back to normal
            else
                P = P./repmat(mbase,[1 size(P,2) 1]);
                %P = 10 .^ (log10(P) - repmat(log10(mbase),[1 size(P,2) 1])); % same as above
            end;
        end;
    end;
    
    PP{ind} = P;
end;
if ~iscell(PPori) PP = PP{1}; end;
