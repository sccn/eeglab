% std_findsameica() - find groups of datasets with identical ICA decomposiotions
%
% Usage: 
%        >> clusters = std_findsameica(ALLEEG);
% Inputs:
%   ALLEEG  - a vector of loaded EEG dataset structures of all sets in the STUDY set.
%
% Outputs:
%   cluster - cell array of groups of datasets
%
% Authors:  Arnaud Delorme, SCCN, INC, UCSD, July 2009-

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, arno@sccn.ucsd.edu
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

% Coding notes: Useful information on functions and global variables used.

function cluster = std_findsameica(ALLEEG);

cluster = { [1] };
for index = 2:length(ALLEEG)
    
    found = 0;
    for c = 1:length(cluster)
        if all(size(ALLEEG(cluster{c}(1)).icaweights) == size(ALLEEG(index).icaweights))
            %if isequal(ALLEEG(cluster{c}(1)).icaweights, ALLEEG(index).icaweights) 
            if sum(sum(abs(ALLEEG(cluster{c}(1)).icaweights-ALLEEG(index).icaweights))) < 2e-5
                cluster{c}(end+1) = index;
                found = 1;
                break;
            end;
        end;
    end;
    if ~found
        cluster{end+1} = index;
    end;
end;
