% eeg_urlatency() - find the real latency of an event in the continuous 
%                   data.
%
% Usage:
%   >> lat_out = eeg_urlatency( event, lat_in );
%
% Inputs:
%   event      - event structure. If this structure contain boundary
%                events, the length of these events is added to restore
%                the original latency from the relative latency in 'lat_in'
%   lat_in     - relative latency.
%
% Outputs:
%   lat_out     - output latency
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, April, 15, 2004
%
% See also: eeg_eegrej()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
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
% Revision 1.1  2004/04/20 01:11:39  arno
% Initial revision
%

function latout = eeg_urlatency( events, latin );
    
    if nargin < 2
        help eeg_urlatency;
        return;
    end;
    
    boundevents = { events.type };
    indbound = strmatch('boundary', boundevents);
    
    latout = latin;
    for index = indbound'
        if events(index).latency < latin
            latout = latout + events(index).length;
        end;
    end;
    