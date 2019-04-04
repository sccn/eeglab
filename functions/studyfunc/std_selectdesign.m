% std_selectdesign() - select an existing STUDY design. 
%                      Use std_makedesign() to add a new STUDY.design.
%
% Usage:
%   >> [STUDY] = std_selectdesign(STUDY, ALLEEG, designind);
%
% Inputs:
%   STUDY     - EEGLAB STUDY structure
%   STUDY     - EEGLAB ALLEEG structure
%   designind - desired (existing) design index
%
% Outputs:
%   STUDY     - EEGLAB STUDY structure with currentdesign set to the input design index. 
%
% Author: Arnaud Delorme, Institute for Neural Computation, UCSD, 2010-

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

function STUDY = std_selectdesign(STUDY, ALLEEG, designind);

if nargin < 3
    help std_selectdesign;
    return;
end

if designind < 1 || designind > length(STUDY.design) || isempty(STUDY.design(designind).name)
    disp('Cannot select an empty STUDY.design');
    return;
end
STUDY.currentdesign = designind;
STUDY = std_rmalldatafields( STUDY );
