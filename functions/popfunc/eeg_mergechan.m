% eeg_mergechan() - merge channel structure while preserving channel
%                   order
%
% >> mergelocs = eeg_mergechan(locs1, locs2);
%
% Inputs: 
%     locs1     - EEGLAB channel location structure
%     locs2     - second EEGLAB channel location structure
%
% Output: 
%     mergelocs - merged channel location structure
%
% Author: Arnaud Delorme, August 2006

% Copyright (C) Arnaud Delorme, CERCO, 2006, arno@salk.edu
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

% union of two channel location structure
% without loosing the order information
% ---------------------------------------
function alllocs = myunion(locs1, locs2)

   labs1 = { locs1.labels };
   labs2 = { locs2.labels };
   
   count1 = 1;
   count2 = 1;
   count3 = 1;
   alllocs = locs1; alllocs(:) = [];
   while count1 <= length(locs1) | count2 <= length(locs2)
       
       if count1 > length(locs1)
           alllocs(count3) = locs2(count2);
           count2 = count2 + 1;
           count3 = count3 + 1;
       elseif count2 > length(locs2)
           alllocs(count3) = locs1(count1);
           count1 = count1 + 1;
           count3 = count3 + 1;
       elseif strcmpi(labs1{count1}, labs2{count2})
           alllocs(count3) = locs1(count1);
           count1 = count1 + 1;
           count2 = count2 + 1;
           count3 = count3 + 1;
       elseif isempty(strmatch(labs1{count1}, labs2))
           alllocs(count3) = locs1(count1);
           count1 = count1 + 1;
           count3 = count3 + 1;
       else
           alllocs(count3) = locs2(count2);
           count2 = count2 + 1;
           count3 = count3 + 1;
       end;
       
   end;
