% loadtxt() - load ascii text file into numeric or cell arrays
%
% Usage:
%   >> array = loadtxt( filename, 'key', 'val' ...);
%
% Inputs:
%    filename - name of the input file
%
% Optional inputs
%   'skipline' - number of lines to skip {default:0}. If this number is
%                negative the program will only skip non-empty lines 
%                (can be usefull for files transmitted from one platform
%                to an other, as CR may be inserted at every lines).
%   'convert'  - 'on' standard text conversion, see note 1
%                'off' no conversion, considers text only
%                'force' force conversion, NaN are returned 
%                for non-numeric inputs {default:'on'}
%   'delim'    - ascii character for delimiters. {default:[9 32]
%                i.e space and tab}. It is also possible to enter 
%                strings, Ex: [9 ' ' ','].
%   'blankcell' - ['on'|'off'] extract blank cells {default:'on'}
%   'verbose'  - ['on'|'off'] {default:'on'}
%   'convertmethod' - ['str2double'|'str2num'] default is 'str2double'
%   'nlines'   - [integer] number of lines to read {default: all file}
%
% Outputs:
%    array - cell array. If the option 'force' is given, the function
%            retrun a numeric array.
%
% Notes: 1) Since it uses cell arrays, the function can handle text input.
%        The function reads each token and then try to convert it to a 
%        number. If the conversion is unsucessfull, the string itself
%        is included in the array.
%        2) The function adds empty entries for rows that contains
%        fewer columns than others.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 29 March 2002

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 29 March 2002
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

function array = loadtxt( filename, varargin );

if nargin < 1
	help loadtxt;
	return;
end;	
if ~isempty(varargin)
   try, g = struct(varargin{:});
   catch, disp('Wrong syntax in function arguments'); return; end
else
    g = [];
end

g = finputcheck( varargin, { 'convert'   'string'   { 'on';'off';'force' }   'on';
                             'skipline'  'integer'  [0 Inf]          0;
                             'verbose'   'string'   { 'on';'off' }   'on';
                             'uniformdelim' 'string'   { 'on';'off' }   'off';                             
                             'blankcell' 'string'   { 'on';'off' }   'on';
                             'convertmethod' 'string'   { 'str2double';'str2num' }   'str2double';
                             'delim'     { 'integer';'string' } []               [9 32];
                             'nlines'    'integer'  []               Inf });
if ischar(g), error(g); end
if strcmpi(g.blankcell, 'off'), g.uniformdelim = 'on'; end
g.convert = lower(g.convert);
g.verbose = lower(g.verbose);
g.delim = char(g.delim);

% open the file
% -------------
if exist(filename) ~=2, error( ['file ' filename ' not found'] ); end;  
fid=fopen(filename,'r','ieee-le');
if fid<0, error( ['file ' filename ' found but error while opening file'] ); end;  

index = 0;
while index < abs(g.skipline)
    tmpline = fgetl(fid); 
    if g.skipline > 0 || ~isempty(tmpline)
        index = index + 1;
    end;    
end; % skip lines ---------

inputline = fgetl(fid);
linenb = 1;
if strcmp(g.verbose, 'on'), fprintf('Reading file (lines): '); end
while isempty(inputline) || inputline(1)~=-1
     colnb = 1;
     if ~isempty(inputline)
         tabFirstpos = 1;
         
         % convert all delimiter to the first one
         if strcmpi(g.uniformdelim, 'on')
             for index = 2:length(g.delim)
                 inputline(find(inputline == g.delim(index))) = g.delim(1);
             end
         end
         
         while ~isempty(deblank(inputline))
             if strcmpi(g.blankcell,'off'), inputline = strtrim(inputline); end
             if tabFirstpos && length(inputline) > 1 && all(inputline(1) ~= g.delim), tabFirstpos = 0; end
             [tmp inputline tabFirstpos] = mystrtok(inputline, g.delim, tabFirstpos);
             switch g.convert
                case 'off', array{linenb, colnb} = tmp;
                case 'on',  
                     if strcmpi(g.convertmethod, 'str2double')
                         tmp2 = str2double(tmp);
                         if isnan( tmp2 )  , array{linenb, colnb} = tmp;
                         else                array{linenb, colnb} = tmp2;
                         end
                     else
                         tmp2 = str2num(tmp);
                         if isempty( tmp2 )  , array{linenb, colnb} = tmp;
                         else                  array{linenb, colnb} = tmp2;
                         end
                     end
                case 'force', array{linenb, colnb} = str2double(tmp);
             end
             colnb = colnb+1;
         end
	     linenb = linenb +1;
     end
     inputline = fgetl(fid);
     if linenb > g.nlines
         inputline = -1;
     end
     if ~mod(linenb,10) && strcmp(g.verbose, 'on'), fprintf('%d ', linenb); end
end;        
if strcmp(g.verbose, 'on'),  fprintf('%d\n', linenb-1); end
if strcmp(g.convert, 'force'), array = [ array{:} ]; end
fclose(fid); 

% problem strtok do not consider tabulation
% -----------------------------------------
function [str, strout, tabFirstpos] = mystrtok(strin, delim, tabFirstpos);
    % remove extra spaces at the beginning
    while any(strin(1) == delim) && strin(1) ~= 9 && strin(1) ~= ','
         strin = strin(2:end);
    end
    % for tab and coma, consider empty cells
    if length(strin) > 1 && any(strin(1) == delim)
        if tabFirstpos || any(strin(2) == delim)
            str = '';
            strout = strin(2:end);
            if strin(2) ~= 9 && strin(2) ~= ','
                tabFirstpos = 0;
                strout = strtrim(strout);
            end
        else
            [str, strout] = strtok(strin, delim);
        end
    else
        [str, strout] = strtok(strin, delim);
    end
