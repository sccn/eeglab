% eeg_decodechan() - given an input EEG dataset structure, output a new EEG data structure 
%                retaining and/or excluding specified time/latency, data point, channel, 
%                and/or epoch range(s).
% Usage:
%   >> [chaninds chanlist] = eeg_decodechan(chanlocs, chanlist);
%
% Inputs:
%   chanlocs  - channel location structure
%   chanlist  - list of channels, numerical indices [1 2 3 ...] or string 
%               'cz pz fz' or cell array { 'cz' 'pz' 'fz' }. Can also be 
%               a list of channel types (see below)
%   field     - ['labels'|'type'] channel field to match. Default is
%               'labels'
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
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [ chaninds, chanlist ] = eeg_decodechan(chanlocs, chanstr, field)

if nargin < 2
    help eeg_decodechan;
    return;
end
if nargin < 3
    field = 'labels';
end

if isempty(chanlocs) && ischar(chanstr)
    chaninds = str2num(chanstr);
    chanlist = chaninds;
    return;
end
    
if ischar(chanstr)
    % convert chanstr
    % ---------------
    chanstr(find(chanstr == ']')) = [];
    chanstr(find(chanstr == '[')) = [];
    chanlistnum = [];
    chanstr  = [ ' ' chanstr ' ' ];
    try
        chanlist = eval( [ '{' chanstr '}' ] );
        chanlistnum = cellfun(@str2double, chanlist);
        if isnumeric(chanlist{1}) chanlist = [ chanlist{:} ]; end
    catch
        chanlist = {};
        sp = find(chanstr == ' ');
        for i = 1:length(sp)-1
            c = chanstr(sp(i)+1:sp(i+1)-1);
            if ~isempty(c)
                chanlist{end+1} = c;
                if isnan(str2double(chanlocs(1).(field))) % channel labels are not numerical
                    if ~isnan(str2double(c))
                        chanlistnum(end+1) = str2double(c);
                    end
                end
            end
        end
        if length(chanlistnum) == length(chanlist)
            chanlist = chanlistnum;
        end
    end
else
    chanlist = chanstr;
end

% convert to values
% -----------------
chanval = 0;
if isnumeric(chanlist)
    chanval = chanlist;
end
% chanval  = [];
% if iscell(chanlist)
%     for ind = 1:length(chanlist)
% 
%         valtmp = str2double(chanlist{ind});
%         if ~isnan(valtmp)
%              chanval(end+1) = valtmp;
%         else chanval(end+1) = 0;
%         end
%     end
% else
%     chanval = chanlist;
% end

% convert to numerical
% --------------------
if all(chanval) > 0
    chaninds = chanval;
    chanlist = chanval;
else
    chaninds = [];
    alllabs  = lower({ chanlocs.(field) });
    chanlist = lower(chanlist);
    for ind = 1:length(chanlist)
        indmatch = find(strcmp(alllabs,chanlist{ind})); %#ok<STCI>
        if ~isempty(indmatch)
            for tmpi = 1:length(indmatch)
                chaninds(end+1) = indmatch(tmpi);
            end
        else
            try,
                eval([ 'chaninds = ' chanlist{ind} ';' ]);
                if isempty(chaninds)
                     error([ 'Channel ''' chanlist{ind} ''' not found' ]);
                else 
                end
            catch
                error([ 'Channel ''' chanlist{ind} ''' not found' ]);
            end
        end
    end
end
chaninds = sort(chaninds);
if ~isempty(chanlocs)
    chanlist = { chanlocs(chaninds).(field) };
else
    chanlist = {};
end
