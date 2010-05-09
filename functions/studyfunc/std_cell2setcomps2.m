% std_cell2setcomps - convert .sets and .comps to cell array. The .sets and
%                     .comps format is useful for GUI but the cell array
%                     format is used for plotting and statistics.
%            
% Usage:
%   [ struct sets comps ] = std_cell2setcomps(STUDY, clustind);
%
% Author: Arnaud Delorme, CERCO/CNRS, UCSD, 2009-

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

% $Log: std_cell2setcomps.m,v $
% Revision 1.1  2009/10/20 02:28:35  arno
% Updated conversion between sets and indices formats
%

function [ tmpstruct setlist complist ] = std_cell2setcomps(STUDY, ALLEEG, setinds, allinds)

if nargin < 4
    tmpstruct = STUDY.cluster(setinds);
    sets = [ STUDY.cluster(setinds).setinds{:} ];
    inds = [ STUDY.cluster(setinds).allinds{:} ];
else
    tmpstruct  = [];
    sets       = [ setinds{:} ];
    inds       = [ allinds{:} ];
end;

% find datasets with common ICA decompositions
clusters = std_findsameica(ALLEEG);

for ind = 1:length(clusters)
    for c = 2:length(clusters{ind})
        sets(find(sets == clusters{ind}(c))) = clusters{ind}(1);
    end;
end;

setlist  = [];
complist = [];
count    = 1;
for i = 1:size(inds,1)
    for j = 1:size(inds,2)
        for ind = 1:length(inds{i,j})
            if ~flag{i,j}(ind)
                
                % found one good component
                complist(count) = inds{i,j}(ind);
                %if complist(count) == 12, dfds; end;
                
                % search for the same component in other datasets
                for c = 1:length(clusters)
                    if any(clusters{c} == sets{i,j}(ind))
                        
                        setlist(:,count)  = clusters{c}';
                        
                        % flag all of these datasets
                        for i2 = 1:size(inds,1)
                            for j2 = 1:size(inds,2)
                                for ind2 = 1:length(sets{i2,j2})
                                    if any(sets{i2,j2}(ind2) == clusters{c}) && complist(count) == inds{i2, j2}(ind2)
                                        flag{i2,j2}(ind2) = 1;
                                    end;
                                end;
                            end;
                        end;
                    end;
                end;
                
                count = count+1;
                
            end;
        end;
    end;
end;
tmpstruct.sets  = setlist;
tmpstruct.comps = complist;
