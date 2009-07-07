%% Copyright (C) 2004 by Alois Schloegl
%%
%% This program is free software; you can redistribute it and/or
%% modify it under the terms of the GNU General Public License
%% as published by the Free Software Foundation; either version 2
%% of the License, or (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with Octave; see the file COPYING.  If not, write to the Free
%% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%% 02110-1301, USA.

%% -*- texinfo -*-
%% @deftypefn {Function File} {@var{idx} =} strfind (@var{str}, @var{pattern})
%% @deftypefnx {Function File} {@var{idx} =} strfind (@var{cellstr}, @var{pattern})
%% Search for @var{pattern} in the string @var{str} and return the
%% starting index of every such occurrence in the vector @var{idx}.
%% If there is no such occurrence, or if @var{pattern} is longer
%% than @var{str}, then @var{idx} is the empty array @code{[]}.
%%
%% If the cell array of strings @var{cellstr} is specified instead of the
%% string @var{str}, then @var{idx} is a cell array of vectors, as specified
%% above.
%% @seealso{findstr, strmatch, strcmp, strncmp, strcmpi, strncmpi}
%% @end deftypefn

%% Author: alois schloegl <a.schloegl@ieee.org>
%% Created: 1 November 2004
%% Adapted-By: William Poetra Yoga Hadisoeseno <williampoetra@gmail.com>

function idx = strfind (text, pattern)

	if (nargin ~= 2)
		usage ('idx = strfind (text, pattern)');
  	elseif (~ischar (pattern))
    		error ('strfind: pattern must be a string value');
  	end

  	lp = length (pattern);
  	if (ischar (text))
    		idx = 1:(length (text) - lp + 1);
  		k = 0;
  		while ((k < lp) && ~isempty (idx))
    	    		idx = idx(text(idx + k) == pattern(k+1));
    	    		k = k+1;
  		end
  	elseif (iscellstr (text))
    		idx = cell (size (text));
    		for i = 1:(numel (text))
			idx{i} = strfind(text{i}, pattern);
    		end
  	else
    		error ('strfind: text must be a string or cell array of strings');
  	end

%% end
