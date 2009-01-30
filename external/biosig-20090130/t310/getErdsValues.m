function values = extractErdsValues(r, time, freq)
% Extracts ERDS values of ERDS maps.
%
% This function returns the ERDS values within a specified time and frequency
% window from an already calculated map.
%
% Usage:
%   values = extractErdsValues(r, time, freq);
%
% Input parameters:
%   r    ... ERDS map structure
%   time ... Time window (s)
%   freq ... Frequency window (Hz)
%
% Output parameter:
%   values ... Structure containing information about the ERDS map

% Copyright by Clemens Brunner, based on timef of EEGLAB
% $Revision: 1.1 $ $Date: 2009-01-30 06:04:44 $
% E-Mail: clemens.brunner@tugraz.at

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

for k = 1:length(r)  % Channels
    values{k} = r{k}.PP(find(r{1}.freqs >= freq(1) & r{1}.freqs <= freq(2)), ...
                        find(r{1}.times >= time(1)*1000 & r{1}.times <= time(2) * 1000));
end;