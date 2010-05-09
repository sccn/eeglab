% sum() - sum of memory mapped underlying array
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, Nov. 2008

% Copyright (C) 2008 Arnaud Delorme, SCCN, INC, UCSD
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

function [s] = sum(a,dim)
    
    if nargin < 2
        dim = 1;
    end;
    
    %b = (:,:,:);
    if ~strcmpi(a.fileformat, 'transposed')
        s = sum(a.data.data.x, dim);
    else
        alldim = [3 1 2];
        if length(size(a)) == 3
            dim = alldim(dim);
            s = sum(a.data.data.x, dim);
            s = permute(s, [3 1 2]);
        else
            if dim == 1
                 dim = 2;
            else dim = 1;
            end;
            s = sum(a.data.data.x, dim)';
        end;
    end;
    return;
    
    % do pnts by pnts if dim = 1
%     if dim == 1 & length(
%         
%         s = zeros(size(a,2), size(a,3));
%         for i=1:size(a,2)
%             s(i,:) = mean(a.data.data.x(:,i,:));
%         end;
%     elseif dim == 1
%          s = zeros(size(a,1), size(a,1));
%         for i=1:size(a,1)
%             s(i,:) = mean(a.data.data.x(:,:,:));
%         end;
%        
%         
%     s = builtin('sum', rand(10,10), dim);
    
    %if length(size(a)) > 2
    %else s = sum(a(:,:,:), dim);
    %end;
