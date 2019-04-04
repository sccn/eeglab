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

function [ tmpstruct setlist complist ] = std_cell2setcomps(STUDY, ALLEEG, setinds, allinds)

if nargin < 4
    tmpstruct = STUDY.cluster(setinds);
    sets = STUDY.cluster(setinds).setinds;
    inds = STUDY.cluster(setinds).allinds;
else
    tmpstruct  = [];
    sets       = setinds;
    inds       = allinds;
end

% initialize flag array
% ---------------------
flag = cell(size(inds));
for i = 1:size(inds,1)
    for j = 1:size(inds,2)
        flag{i,j} = zeros(size(inds{i,j}));
    end
end

% find datasets with common ICA decompositions
clusters = std_findsameica(ALLEEG);

setlist  = [];
complist = [];
count    = 1;
for i = 1:size(inds,1)
    for j = 1:size(inds,2)
        for ind = 1:length(inds{i,j})
            if ~flag{i,j}(ind)
                
                % found one good component
                complist(count) = inds{i,j}(ind);
                %if complist(count) == 12, dfds; end
                
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
                                    end
                                end
                            end
                        end
                    end
                end
                
                count = count+1;
                
            end
        end
    end
end
tmpstruct.sets  = setlist;
tmpstruct.comps = complist;
