% findboundaries() - Find boundaries (data discontinuities) in event
%                    structure of continuous EEG dataset
%
% Usage:
%   >> boundaries = findboundaries(EEG.event);
%
% Inputs:
%   EEG.event     - EEGLAB EEG event structure
%
% Outputs:
%   boundaries    - scalar or vector of boundary event latencies
%
% Author: Andreas Widmann, University of Leipzig, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

function boundaries = findboundaries(event)

if isfield(event, 'type') & isfield(event, 'latency') & cellfun('isclass', {event.type}, 'char')

    % Boundary event indices
    boundaries = strmatch('boundary', {event.type});

    % Boundary event latencies
    boundaries = [event(boundaries).latency];

    % Shift boundary events to epoch onset
    boundaries = fix(boundaries + 0.5);

    % Remove duplicate boundary events
    boundaries = unique(boundaries);

    % Epoch onset at first sample?
    if isempty(boundaries) || boundaries(1) ~= 1
        boundaries = [1 boundaries];
    end

else

    boundaries = 1;

end
