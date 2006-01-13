% biosig2eeglabevent() - convert biosig events to EEGLAB event structure
%
% Usage:
%   >> eeglabevent = biosig2eeglabevent( biosigevent )
%
% Inputs:
%   biosigevent    - BioSig event structure
%
% Outputs:
%   eeglabevent    - EEGLAB event structure
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2006-

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $

function event = biosig2eeglabevent(EVENT)

disp('Importing data events...');
    if isfield(EVENT, 'Teeg')
        event = EVENT.Teeg;
    end
    if isfield(EVENT, 'TYP')
        for index = 1:length( EVENT.TYP )
            event(index).type = EVENT.TYP(index);
        end;
    end;
    if isfield(EVENT, 'POS')
        for index = 1:length( EVENT.POS )
            event(index).latency = EVENT.POS(index);
        end;
    end;
    if isfield(EVENT, 'DUR')
        if any( [ EVENT.DUR ] )
            for index = 1:length( EVENT.DUR )
                event(index).duration = EVENT.DUR(index);
            end;
        end;
    end;
    if isfield(EVENT, 'CHN')
        if any( [ EVENT.CHN ] )
            for index = 1:length( EVENT.CHN )
                event(index).chanindex = EVENT.CHN(index);
            end;
        end;
    end;            
