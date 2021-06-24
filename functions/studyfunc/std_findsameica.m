% std_findsameica() - find groups of datasets with identical ICA decomposiotions
%                     (search identical weight*sphere matrices)
%
% Usage: 
%        >> clusters = std_findsameica(ALLEEG);
%        >> clusters = std_findsameica(ALLEEG,icathreshold);
% Inputs:
%   ALLEEG           - a vector of loaded EEG dataset structures of all sets 
%                      in the STUDY set.
%   icathreshold     - Threshold to compare icaweights. Default 2e-4.
%
% Outputs:
%   cluster - cell array of groups of datasets
%   indices - cluster index for each dataset
%
% Authors:  Arnaud Delorme, SCCN, INC, UCSD, July 2009-
% 2016 change: as of May 2016, the function now compares the product of the
%              weight and the sphere matrices instead of just the weight
%              matrices.

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, arno@sccn.ucsd.edu
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

% Coding notes: Useful information on functions and global variables used.

function [cluster, inds] = std_findsameica(ALLEEG, varargin)

% 6/2/2014 Ramon : Allow ica threshold as input.
if nargin == 1
    icathreshold = 2e-4;
elseif nargin == 2
    icathreshold = varargin{1};
end
    
cluster = { [1] };
inds = [1];
for index = 2:length(ALLEEG)
    
    found = 0;
    for c = 1:length(cluster)
        w1 = ALLEEG(cluster{c}(1)).icaweights*ALLEEG(cluster{c}(1)).icasphere;
        w2 = ALLEEG(index).icaweights*ALLEEG(index).icasphere;
        if all(size(w1) == size(w2))
            %if isequal(ALLEEG(cluster{c}(1)).icaweights, ALLEEG(index).icaweights) 
            if sum(sum(abs(w1-w2))) < icathreshold
                cluster{c}(end+1) = index;
                inds(index) = c;
                found = 1;
                break;
            end
        end
    end
    if ~found
        cluster{end+1} = index;
        inds(index) = index;
    end
end
