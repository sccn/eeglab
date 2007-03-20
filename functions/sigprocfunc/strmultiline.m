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

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% $Log: not supported by cvs2svn $
% Revision 1.4  2006/09/20 12:22:37  arno
% multiple line fix
%
% Revision 1.3  2005/09/27 22:02:04  arno
% allow processing multiple lines
%
% Revision 1.2  2005/03/20 18:39:14  scott
% help msg
%
% Revision 1.1  2005/03/09 19:15:22  arno
% Initial revision
%

function strout = strmultiline( strinori, maxlen, delimiter);

if nargin < 1
	help strmultiline;
	return;
end;	
if nargin < 3
    delimiter = [];
end;

% first format input string
% -------------------------
if ~isempty(find(strinori == 10))
    tmpinds = [1 find(strinori == 10)+1 length(strinori)];
    strinori = [ ' ' strinori ' ' ];
    for ind = 2:length(tmpinds)
        allstrs{ind} = strinori(tmpinds(ind-1)+1:tmpinds(ind)-1);
        if isempty(allstrs{ind}), allstrs{ind} = ' '; end;
    end;
    strinori = strvcat(allstrs{:});
end;
if size(strinori,2) < maxlen, strout = strinori; return; end;

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
        end;
        if isempty(curline) curline = tok;
        else                curline = [ curline ' ' tok ];
        end;
    end;
    if ~isempty(curline), lines{count} = curline; end;
    
    % type of delimiter
    % -----------------
    if isempty(delimiter)
        strouttmp = strvcat(lines{:});
        if isempty(strouttmp)
            strouttmp = ones(1,maxlen)*' ';
        elseif size(strouttmp, 2) < maxlen
            strouttmp(:,end+1:maxlen) = ' ';
        end;
        strout = strvcat(strout, strouttmp);
    else
        strouttmp = lines{1};
        for index = 2:length(lines)
            strouttmp = [ strouttmp 10 lines{index} ];
        end;
        if index == 1, strout = strouttmp;
        else strout = [ strout 10 strouttmp ];
        end;
    end;
end;
	 
