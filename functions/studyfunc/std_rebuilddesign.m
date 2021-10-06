% std_rebuilddesign - reduild design structure when datasets have been 
%                     removed or added.
% Usage:
%    STUDY = std_rebuilddesign(STUDY, ALLEEG);
%    STUDY = std_rebuilddesign(STUDY, ALLEEG, designind);
%
% Inputs:
%   STUDY     - EEGLAB STUDY set
%   ALLEEG    - vector of the EEG datasets included in the STUDY structure 
%
% Optional inputs:
%   designind - [integer>0] indices (number) of the design to rebuild.
%               Default is all.
% Output:
%   STUDY     - updated EEGLAB STUDY set
%
% Author: Arnaud Delorme, Institute for Neural Computation UCSD, 2010-

% Copyright (C) Arnaud Delorme, arno@ucsd.edu
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

function STUDY = std_rebuilddesign(STUDY, ALLEEG, designind)

if nargin < 2
    help std_rebuilddesign;
    return;
end

if nargin < 3
    designind = 1:length(STUDY.design);
end

disp('Rebuilding designs...');
[indvars, indvarvals] = std_getindvar(STUDY);
for indDesign = designind
    
    % find out if some independent variables or independent 
    % variable values have been removed
    tmpdesign = STUDY.design(indDesign);
    if isempty(tmpdesign.cases.value) || isempty(tmpdesign.cases.value{1})
        tmpdesign.cases.value = STUDY.subject;
    else
        tmpdesign.cases.value = intersect(tmpdesign.cases.value, STUDY.subject);
    end
    for iVar = length(tmpdesign.variable):-1:1
        indVar = strmatch(tmpdesign.variable(iVar).label, indvars, 'exact');
        if isempty(indVar)
            fprintf('Design %d: removing independent variable "%s"\n', indDesign, tmpdesign.variable(iVar).label);
            tmpdesign.variable(iVar) = []; % remove independent var
        else
            tmpvals = myintersect( tmpdesign.variable(iVar).value, indvarvals{indVar} );
            if isempty(tmpvals)
                 fprintf('Design %d: removing independent variable "%s"\n', indDesign, tmpdesign.variable(iVar).label);
                 tmpdesign.variable(iVar) = []; % remove independent var
            else tmpdesign.variable(iVar).value = tmpvals;
            end
        end
    end
    
    STUDY = std_makedesign(STUDY, ALLEEG, indDesign, tmpdesign);
end

STUDY = std_selectdesign(STUDY, ALLEEG, STUDY.currentdesign);

% take the intersection for independent variables
% -----------------------------------------------
function a = myintersect(a,b)

    if isempty(b) || isempty(a), a = {}; return; end
    
    for index = length(a):-1:1
        if ischar(a{index})
            res = intersect_bc(a(index), b);
            if isempty(res)
                 a{index} = [];
            else a(index) = res;
            end
        elseif iscell(a{index})
            a{index} = intersect_bc(a{index}, b);
        elseif isnumeric(a{index})
            a{index} = intersect_bc(a{index}, [ b{:} ]);
        end
        if isempty(a{index}), a(index) = []; end
    end
        
