% strmultiline() - format a long string as a multi-line string. 
%
% Usage:
%   >> strout = strmultiline( strin, maxlen, delimiter);
%
% Inputs:
%   strin     - one-line or several-line string
%   maxlen    - maximum line length
%   delimiter - enter 10 here to separate lines with character 10. Default is
%               empty: lines are on different rows of the returned array.
% Outputs:
%   strout    - string with multiple line
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2005
%
% See also: eegfilt(), eegfiltfft(), eeglab()

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

function strout = strmultiline( strinori, maxlen, delimiter);

if nargin < 1
	help strmultiline;
	return;
end;	
if nargin < 3
    delimiter = [];
end

% first format input string
% -------------------------
if ~isempty(find(strinori == 10))
    tmpinds = [1 find(strinori == 10)+1 length(strinori)];
    strinori = [ ' ' strinori ' ' ];
    for ind = 2:length(tmpinds)
        allstrs{ind} = strinori(tmpinds(ind-1)+1:tmpinds(ind)-1);
        if isempty(allstrs{ind}), allstrs{ind} = ' '; end
    end
    strinori = strvcat(allstrs{:});
end
if size(strinori,2) < maxlen, strout = strinori; return; end

strout = [];
for index = 1:size(strinori,1) % scan lines
    strin = deblank(strinori(index, :));
    lines = {};
    count = 1;
    curline = '';
    
    while ~isempty( strin )
        [tok strin] = strtok(strin);
        if length(curline) + length(tok) +1 > maxlen
            lines{count} = curline;
            curline      = '';
            count = count + 1;
        end
        if isempty(curline) curline = tok;
        else                curline = [ curline ' ' tok ];
        end
    end
    if ~isempty(curline), lines{count} = curline; end
    
    % type of delimiter
    % -----------------
    if isempty(delimiter)
        if ~isempty(lines)
             strouttmp = strvcat(lines{:});
        else strouttmp = '';
        end
        if isempty(strouttmp)
            strouttmp = ones(1,maxlen)*' ';
        elseif size(strouttmp, 2) < maxlen
            strouttmp(:,end+1:maxlen) = ' ';
        end
        strout = strvcat(strout, strouttmp);
    else
        strouttmp = lines{1};
        for index = 2:length(lines)
            strouttmp = [ strouttmp 10 lines{index} ];
        end
        if index == 1, strout = strouttmp;
        else strout = [ strout 10 strouttmp ];
        end
    end
end
	 
