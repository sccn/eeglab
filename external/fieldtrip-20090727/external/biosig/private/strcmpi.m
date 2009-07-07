%% Copyright (C) 2000  Bill Lash
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

%% usage: strcmpi (s1, s2)
%%
%% Compare two strings, ignoring case, returning 1 if
%% they are the same, and 0 otherwise.
%%
%% Note: For compatibility with Matlab, Octave's strcmpi function
%% returns 1 if the strings are equal, and 0 otherwise.  This is
%% just the opposite of the corresponding C library function.

%% Author: Bill Lash <lash@tellabs.com>

function status = strcmpi(s1, s2)

  if (nargin ~= 2)
    usage ('strcmpi (s, t)');
  end

  status = 0;		%% Assume strings are different
  if (isstr (s1) && isstr(s2))
    status = strcmp(upper(s1),upper(s2));
  end

end

