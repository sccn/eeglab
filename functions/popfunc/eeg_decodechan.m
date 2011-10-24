% eeg_decodechan() - given an input EEG dataset structure, output a new EEG data structure 
%                retaining and/or excluding specified time/latency, data point, channel, 
%                and/or epoch range(s).
% Usage:
%   >> [chaninds chanlist] = eeg_decodechan(chanlocs, chanlist);
%
% Inputs:
%   chanlocs  - channel location structure
%   chanlist  - list of channels, numerical indices [1 2 3 ...] or string 
%               'cz pz fz' or cell array { 'cz' 'pz' 'fz' }
%
% Outputs:
%   chaninds  - integer array with the list of channel indices
%   chanlist  - cell array with a list of channel labels
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, 2009-
% 
% see also: eeglab()

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

function [ chaninds chanlist ] = eeg_decodechan(chanlocs, chanstr);

if nargin < 2
    help eeg_decodechan;
    return;
end;

if isempty(chanlocs) && isstr(chanstr)
    chaninds = str2num(chanstr);
    chanlist = chaninds;
    return;
end;
    
if isstr(chanstr)
    % convert chanstr
    % ---------------
    chanstr(find(chanstr == ']')) = [];
    chanstr(find(chanstr == '[')) = [];
    chanlistnum = [];
    chanstr  = [ ' ' chanstr ' ' ];
    chanlist = {};
    sp = find(chanstr == ' ');
    for i = 1:length(sp)-1
        c = chanstr(sp(i)+1:sp(i+1)-1);
        if ~isempty(c)
            chanlist{end+1} = c;
            if isnan(str2double(chanlocs(1).labels)) % channel labels are not numerical
                if ~isnan(str2double(c))
                    chanlistnum(end+1) = str2double(c);
                end;
            end;
        end;
    end;
    if length(chanlistnum) == length(chanlist)
        chanlist = chanlistnum;
    end;
else
    chanlist = chanstr;
end;

% convert to values
% -----------------
chanval = 0;
if isnumeric(chanlist)
    chanval = chanlist;
end;
% chanval  = [];
% if iscell(chanlist)
%     for ind = 1:length(chanlist)
% 
%         valtmp = str2double(chanlist{ind});
%         if ~isnan(valtmp)
%              chanval(end+1) = valtmp;
%         else chanval(end+1) = 0;
%         end;
%     end;
% else
%     chanval = chanlist;
% end;

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
        indmatch = find(strcmp(alllabs,chanlist{ind})); %#ok<STCI>
        if ~isempty(indmatch)
            for tmpi = 1:length(indmatch)
                chaninds(end+1) = indmatch(tmpi);
            end;
        else
            try,
                eval([ 'chaninds = ' chanlist{ind} ';' ]);
                if isempty(chaninds)
                     error([ 'Channel ''' chanlist{ind} ''' not found' ]);
                else 
                end;
            catch
                error([ 'Channel ''' chanlist{ind} ''' not found' ]);
            end;
        end;
    end;
end;
chaninds = sort(chaninds);
if ~isempty(chanlocs)
    chanlist = { chanlocs(chaninds).labels };
else
    chanlist = {};
end;
