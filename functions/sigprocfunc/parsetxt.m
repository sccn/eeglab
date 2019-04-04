% parsetxt() - parse text input into cell array
%
% Usage: >> cellarray = parsetxt( txt, delims );
%
% Inputs: 
%    txt    - input text
%    delims - optional char array of delimiters (default: [' ' ',' 9]);
%
% Note: commas, and simple quotes are ignored
%
% Author: Arnaud Delorme, CNL / Salk Institute, 18 April 2002

% Copyright (C) 18 April 2002 Arnaud Delorme, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

function cellarray = parsetxt(txt, delims);

if nargin < 1
	help parsetxt;
	return;
end;	

if nargin < 2
    delims = [' ' ',' 9 '"' '''' ];
end
    
cellarray = {};
tmptxt = '';
for index =1:length(txt)
    if ~isempty(findstr(txt(index), delims))
        if ~isempty(tmptxt), cellarray = { cellarray{:}, tmptxt }; end
        tmptxt = '';
    else
        tmptxt = [ tmptxt txt(index) ];
    end
end;        
if ~isempty(tmptxt)
    cellarray = { cellarray{:}, tmptxt };
end
return;
