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
% Revision 1.5  2004/06/11 22:50:24  arno
% debug duration field
%
% Revision 1.4  2004/05/06 21:54:10  arno
% debug last
%
% Revision 1.3  2004/05/06 21:39:36  arno
% test event format (string?)
%
% Revision 1.2  2004/04/20 02:08:48  arno
% debug
%
% Revision 1.1  2004/04/20 01:11:39  arno
% Initial revision
%

function latout = eeg_urlatency( events, latin );
    
    if nargin < 2
        help eeg_urlatency;
        return;
    end;
    
    boundevents = { events.type };
    latout      = latin;
    if ~isempty(boundevents) & isstr(boundevents{1})
        indbound = strmatch('boundary', boundevents);
        
        if isfield(events, 'duration') & ~isempty(indbound)
            for index = indbound'
                if events(index).latency < latin
                    latout = latout + events(index).duration;
                end;
            end;
        elseif ~isempty(indbound) % boundaries but with no associated duration
            latout = NaN;
        end;
    end;
    