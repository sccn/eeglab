% lookupchantemplate - look up channel template.
%
% Usage:
%  [ found transform ] = lookupchantemplate( filename, template_struct);
%
% Inputs:
%   filename        - channel location file name
%   template_struct - template strcuture as defined in dipfitdefs
%
% Outputs:
%   found     - [0|1] 1 if a transformation was found for this template
%   transform - [array] tranditional tailairach transformation
%
% Author: Arnaud Delorme, SCCN, INC, 2007

% Copyright (C) 2007 Arnaud Delorme
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


