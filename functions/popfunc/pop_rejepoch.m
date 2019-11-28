% pop_rejepoch() - Reject pre-labeled trials in a EEG dataset. 
%                   Ask for confirmation and accept the rejection
%
% Usage:
%         >> OUTEEG = pop_rejepoch( INEEG, trialrej, confirm)
%
% Inputs:
%   INEEG      - Input dataset
%   trialrej   - Array of 0s and 1s (depicting rejected trials) (size is 
%                number of trials)
%   confirm    - Display rejections and ask for confirmation. (0=no. 1=yes;
%                default is 1).
% Outputs:
%   OUTEEG     - output dataset
%
% Example:
%          >> data2 = pop_rejepoch( EEG, [1 1 1 0 0 0] );
%             % reject the 3 first trials of a six-trial dataset
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab(), eegplot()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% 01-25-02 reformated help & license -ad 

function [EEG, com] = pop_rejepoch( EEG, tmprej, confirm);

com = '';
if nargin < 1
   help pop_rejepoch;
   return;
end
if nargin < 2
    tmprej =  find(EEG.reject.rejglobal);
end
if nargin < 3
   confirm = 1;
end
if islogical(tmprej), tmprej = tmprej+0; end

uniquerej = double(sort(unique(tmprej)));
if length(tmprej) > 0 && length(uniquerej) <= 2 && ...
    ismember(uniquerej(1), [0 1]) && ismember(uniquerej(end), [0 1]) && any(~tmprej)
    format0_1 = 1;
    fprintf('%d/%d trials rejected\n', sum(tmprej), EEG.trials);
else 
    format0_1 = 0;
    fprintf('%d/%d trials rejected\n', length(tmprej), EEG.trials);
end

if confirm ~= 0
    ButtonName=questdlg2('Are you sure, you want to reject the labeled trials ?', ...
                         'Reject pre-labelled epochs -- pop_rejepoch()', 'NO', 'YES', 'YES');
    switch ButtonName,
        case 'NO', 
        	disp('Operation cancelled');
 			return;
       case 'YES',
       		disp('Compute new dataset');
    end % switch

end

% create a new set if set_out is non nul 
% --------------------------------------
if format0_1 || length(tmprej) == EEG.trials
    tmprej = find(tmprej);
end
EEG = pop_select( EEG, 'notrial', tmprej);

com = sprintf( 'EEG = pop_rejepoch( EEG, %s);', vararg2str({ tmprej 0 }));		
return;
