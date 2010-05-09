% biosig2eeglabevent() - convert biosig events to EEGLAB event structure
%
% Usage:
%   >> eeglabevent = biosig2eeglabevent( biosigevent, interval )
%
% Inputs:
%   biosigevent    - BioSig event structure
%   interval       - Period to extract events for, in frames.
%                    Default [] is all.
%
% Outputs:
%   eeglabevent    - EEGLAB event structure
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2006-

% Copyright (C) 13 2006- Arnaud Delorme, Salk Institute, arno@salk.edu
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

function event = biosig2eeglabevent(EVENT, interval)

if nargin < 2
    interval = [];
end;

event = [];
disp('Importing data events...');

% If the interval variable is empty, import all events.
if isempty(interval)
    if isfield(EVENT, 'Teeg')
        event = EVENT.Teeg;
    end
    if isfield(EVENT, 'TYP')
        for index = 1:length( EVENT.TYP )
            event(index).type = EVENT.TYP(index);
        end
    end
    if isfield(EVENT, 'POS')
        for index = 1:length( EVENT.POS )
            event(index).latency = EVENT.POS(index);
        end
    end
    if isfield(EVENT, 'DUR')
        if any( [ EVENT.DUR ] )
            for index = 1:length( EVENT.DUR )
                event(index).duration = EVENT.DUR(index);
            end
        end
    end
    if isfield(EVENT, 'CHN')
        if any( [ EVENT.CHN ] )
            for index = 1:length( EVENT.CHN )
                event(index).chanindex = EVENT.CHN(index);
            end
        end
    end
% If a subinterval was specified, select only events that fall in that range, and
% edit duration field if it exceeds that range.
elseif isfield(EVENT,'POS')
    count = 1;
    for index = 1:length(EVENT.POS)
        pos_tmp = EVENT.POS(index) - interval(1) + 1;
        if pos_tmp > 0 & EVENT.POS(index) <= interval(2)
            event(count).latency = pos_tmp;
            if isfield(EVENT, 'TYP'), event(count).type = EVENT.TYP(index); end
            if isfield(EVENT, 'CHN'), event(count).chanindex = EVENT.CHN(index); end
            if isfield(EVENT, 'DUR')
                event(count).duration = min(EVENT.DUR(index), interval(2) - EVENT.POS(index));
            end
            count = count + 1;
        end
    end
end
