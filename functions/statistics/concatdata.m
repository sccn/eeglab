% concatdata - concatenate data stored into a cell array into a single
%              array. only concatenate along the last dimension
% Usage:
%    [dataarray len dims] = concatata(cellarraydata);
%
% Input:
%    cellarraydata - cell array containing data
%
% Output:
%    dataarray - single array containing all data
%    len       - limits of each array
%    dim       - dimension of the original array
%
% Example:
%   a = rand(3, 4, 3, 10);
%   b = rand(3, 4, 3, 4);
%   c = rand(3, 4, 3, 3);
%   [ alldata len ] = concatdata({ a b c});
%   % alldata is size [ 3 4 3 17 ]
%   % len contains [ 0 10 14 17 ]
%   % to access array number i, type "alldata(len(i)+1:len(i+1))
%
% Author: Arnaud Delorme, CERCO/CNRS & SCCN/INC/UCSD, 2009-

% Copyright (C) Arnaud Delorme
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

function [ datac, alllen, dims ] = concatdata(data);
    
    alllen    = cellfun('size', data, myndims(data{1}) ); % by chance, pick up the last dimension
    dims      = size(data);
    alllen    = [ 0 alllen(:)' ];
    switch myndims(data{1})
        case 1,
            datac = zeros(sum(alllen),1, 'single');
            
            for i = 1:prod(dims)
                alllen(i+1) = alllen(i+1) + alllen(i);
                datac(alllen(i)+1:alllen(i+1)) = data{i};
            end
        case 2,
            datac = zeros(size(data{1},1), sum(alllen), 'single');
            
            for i = 1:prod(dims)
                alllen(i+1) = alllen(i+1) + alllen(i);
                datac(:,alllen(i)+1:alllen(i+1)) = data{i};
            end
        case 3,
            datac = zeros(size(data{1},1), size(data{1},2), sum(alllen), 'single');
            
            for i = 1:prod(dims)
                alllen(i+1) = alllen(i+1) + alllen(i);
                datac(:,:,alllen(i)+1:alllen(i+1)) = data{i};
            end
        case 4,
            datac = zeros(size(data{1},1), size(data{1},2), size(data{1},3), sum(alllen), 'single');
            
            for i = 1:prod(dims)
                alllen(i+1) = alllen(i+1) + alllen(i);
                datac(:,:,:,alllen(i)+1:alllen(i+1)) = data{i};
            end
        case 5,
            datac = zeros(size(data{1},1), size(data{1},2), size(data{1},3), size(data{1},4), sum(alllen), 'single');
            
            for i = 1:prod(dims)
                alllen(i+1) = alllen(i+1) + alllen(i);
                datac(:,:,:,:,alllen(i)+1:alllen(i+1)) = data{i};
            end
        case 6,
            datac = zeros(size(data{1},1), size(data{1},2), size(data{1},3), size(data{1},4),size(data{1},5), sum(alllen), 'single');
            
            for i = 1:prod(dims)
                alllen(i+1) = alllen(i+1) + alllen(i);
                datac(:,:,:,:,:,alllen(i)+1:alllen(i+1)) = data{i};
            end
        case 7,
            datac = zeros(size(data{1},1), size(data{1},2), size(data{1},3), size(data{1},4), size(data{1},5), size(data{1},6), sum(alllen), 'single');
            
            for i = 1:prod(dims)
                alllen(i+1) = alllen(i+1) + alllen(i);
                datac(:,:,:,:,:,:,alllen(i)+1:alllen(i+1)) = data{i};
            end
    end

    
function val = myndims(a)
    if ndims(a) > 2
        val = ndims(a);
    else
        if size(a,1) == 1,
            val = 2;
        elseif size(a,2) == 1,
            val = 1;
        else
            val = 2;
        end
    end; 
            
