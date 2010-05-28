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
%    dim       - dimension of the orginal array
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
            end;
        case 2,
            datac = zeros(size(data{1},1), sum(alllen), 'single');
            
            for i = 1:prod(dims)
                alllen(i+1) = alllen(i+1) + alllen(i);
                datac(:,alllen(i)+1:alllen(i+1)) = data{i};
            end;
        case 3,
            datac = zeros(size(data{1},1), size(data{1},2), sum(alllen), 'single');
            
            for i = 1:prod(dims)
                alllen(i+1) = alllen(i+1) + alllen(i);
                datac(:,:,alllen(i)+1:alllen(i+1)) = data{i};
            end;
        case 4,
            datac = zeros(size(data{1},1), size(data{1},2), size(data{1},3), sum(alllen), 'single');
            
            for i = 1:prod(dims)
                alllen(i+1) = alllen(i+1) + alllen(i);
                datac(:,:,:,alllen(i)+1:alllen(i+1)) = data{i};
            end;
        case 5,
            datac = zeros(size(data{1},1), size(data{1},2), size(data{1},3), size(data{1},4), sum(alllen), 'single');
            
            for i = 1:prod(dims)
                alllen(i+1) = alllen(i+1) + alllen(i);
                datac(:,:,:,:,alllen(i)+1:alllen(i+1)) = data{i};
            end;
    end;

    
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
        end;
    end; 
            
