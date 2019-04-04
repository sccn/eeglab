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
       end
       
   end
