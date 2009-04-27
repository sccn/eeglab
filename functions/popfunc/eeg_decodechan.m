% eeg_decodechan() - given an input EEG dataset structure, output a new EEG data structure 
%                retaining and/or excluding specified time/latency, data point, channel, 
%                and/or epoch range(s).
% Usage:
%   >> chanlist = eeg_decodechan(chanlocs, chanlist);
%
% Inputs:
%   chanlocs  - channel location structure
%   chanlist  - list of channels 'cz pz fz' etc or { 'cz' 'pz' 'fz' }
%
% Outputs:
%   chanlist  - cell array with a list of channel
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, 2009-
% 
% see also: eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [ chaninds chanlist ] = eeg_decodechan(chanlocs, chanstr);

if nargin < 2
    help eeg_decodechan;
    return;
end;

if isstr(chanstr)
    % convert chanstr
    % ---------------
    chanstr  = [ ' ' chanstr ' ' ];
    chanlist = {};
    sp = find(chanstr == ' ');
    for i = 1:length(sp)-1
        c = chanstr(sp(i)+1:sp(i+1)-1);
        if ~isempty(c)
            chanlist{end+1} = c;
        end;
    end;
else
    chanlist = chanstr;
end;

% convert to values
% -----------------
chanval  = [];
if iscell(chanlist)
    for ind = 1:length(chanlist)

        valtmp = str2num(chanlist{ind});
        if ~isempty(valtmp)
             chanval(end+1) = valtmp;
        else chanval(end+1) = 0;
        end;
    end;
else
    chanval = chanlist;
end;

% convert to numerical
% --------------------
if all(chanval) > 0
    chaninds = chanval;
    chanlist = chanval;
else
    chaninds = [];
    alllabs  = lower({ chanlocs.labels });
    chanlist = lower(chanlist);
    for ind = 1:length(chanlist)
        indmatch = strmatch(chanlist{ind}, alllabs, 'exact');
        if ~isempty(indmatch)
            for tmpi = 1:length(indmatch)
                chaninds(end+1) = indmatch(tmpi);
            end;
        else error([ 'Channel ''' chanlist{ind} ''' not found' ]);
        end;
    end;
end;
chaninds = sort(chaninds);
