% std_checkdesign() - Check if STUDY design is compatible with EEGLAB
%                     2-categorical factor statistics
% Usage:
%   >> [STUDY] = std_checkdesign(STUDY, designind);
%
% Inputs:
%   STUDY      - EEGLAB STUDY set
%   designind  - [integer] design index
%
% Authors: Arnaud Delorme

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2017
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

function retVal = std_checkdesign(STUDY, designind)

if nargin < 2
    help std_checkdesign;
    return;
end

retVal = 1;
if length(STUDY.design(designind).variable) > 2 && any(cellfun(@(x)length(x)>1, {STUDY.design(designind).variable(3:end).value} )) ... % it is OK to have more than 2 var if single value
    || ( ~isempty(STUDY.design(designind).variable) && any(strmatch('continuous', {STUDY.design(designind).variable.vartype })))   % not OK to have continuous variables
    warndlg2( [ 'These plotting functions cannot process this design.' 10 ...
        'This design either has more than 2 categorical variables' 10 ...
        'or it uses one continuous variable. To process this' 10 ...
        'design, use the new LIMO EEG tools.' ]);
    retVal = 0;
end

