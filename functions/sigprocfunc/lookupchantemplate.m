% lookupchantemplate - look up channel template.
%
% Usage:
%  [ found transform ] = lookupchantemplate( filename, template_struct);
%
% Inputs:
%   filename        - channel location file name
%   template_struct - template structure as defined in dipfitdefs
%
% Outputs:
%   found     - [0|1] 1 if a transformation was found for this template
%   transform - [array] tranditional tailairach transformation
%
% Author: Arnaud Delorme, SCCN, INC, 2007

% Copyright (C) 2007 Arnaud Delorme
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

function [allkeywordstrue, transform] = lookupchantemplate(chanfile, tmpl);

if nargin < 2
    help lookupchantemplate;
    return;
end

allkeywordstrue = 0;
transform = [];
for ind = 1:length(tmpl)
    allkeywordstrue = 1;
    if isempty(tmpl(ind).keywords), allkeywordstrue = 0; end
    for k = 1:length(tmpl(ind).keywords)
        if isempty(findstr(lower(chanfile), lower(tmpl(ind).keywords{k}))), allkeywordstrue = 0; end
    end
    if allkeywordstrue,
        transform = tmpl(ind).transform;
        break;
    end
end


