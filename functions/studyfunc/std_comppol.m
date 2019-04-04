% std_comppol() - inverse component polarity in a component cluster
%
% Usage: [compout pol] = std_comppol(compin);
%
% Inputs:
%    compin  - component scalp maps, one per column.
%
% Outputs:
%    compout - component scalp maps some of them with inverted
%              polarities, one per column.
%    pol     - logical vector of component with inverted 
%              polarities (same length as the number of rows in 
%              compin)
%
% Author: Arnaud Delorme & Hilit Serby, SCCN, INC, UCSD, 2004

% Copyright (C) 2004 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [compin, pol] = std_comppol(compin);

if nargin < 1
    help std_comppol;
    return;
end

% remove the NaN
% --------------
for index = 1:size(compin,2)
    compin(isnan(compin(:,index)),:) =[];
end

% run several iterations
% ----------------------
pol     = ones(1,size(compin,2));
for repeat=1:3
    compave = mean(compin,2);
    for index = 1:size(compin,2)
        
        % remove diagonal and put 0 and 1
        % -------------------------------
        if ~all(compin(:,index) == 0)
             r = corrcoef(compave, compin(:,index) );
        else r = zeros(2,2);
        end
        
        % invert component polarities
        % ---------------------------
        if r(2) < 0
            compin(:,index) = -compin(:,index);
            pol(index)      = -pol(index);
        end
    end
end
