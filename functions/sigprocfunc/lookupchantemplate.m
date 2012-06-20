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

function [allkeywordstrue, transform] = lookupchantemplate(chanfile, tmpl);

if nargin < 2
    help lookupchantemplate;
    return;
end;

allkeywordstrue = 0;
transform = [];
for ind = 1:length(tmpl)
    allkeywordstrue = 1;
    if isempty(tmpl(ind).keywords), allkeywordstrue = 0; end;
    for k = 1:length(tmpl(ind).keywords)
        if isempty(findstr(lower(chanfile), lower(tmpl(ind).keywords{k}))), allkeywordstrue = 0; end;
    end;
    if allkeywordstrue,
        transform = tmpl(ind).transform;
        break;
    end;
end;


