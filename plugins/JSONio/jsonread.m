function json = jsonread(filename, opts)
% JSON (JavaScript Object Notation) parser - a compiled routine
% FORMAT json = jsonread(filename, opts)
% filename - name of a JSON file or JSON string
% json     - JSON structure
% opts     - structure of optional parameters:
%              replacementStyle: string to control how non-alphanumeric
%                characters are replaced {'underscore','hex','delete','nop'}
%                [Default: 'underscore']
% 
% References:
%   JSON Standard: http://www.json.org/
%   JSMN C parser: http://zserge.com/jsmn.html
%   jsondecode: http://www.mathworks.com/help/matlab/ref/jsondecode.html

% Guillaume Flandin
% $Id: jsonread.m 7045 2017-03-17 10:41:12Z guillaume $


%-This is merely the help file for the compiled routine
error('jsonread.c not compiled.');

% mex jsonread.c jsmn.c -DJSMN_PARENT_LINKS
