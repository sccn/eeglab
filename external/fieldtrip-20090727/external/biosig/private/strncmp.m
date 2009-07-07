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

%% usage: strncmp (s1, s2, n)
%%
%% Compare the first n characters of two strings, returning 1 if
%% they are the same, and 0 otherwise.
%%
%% Note: For compatibility with Matlab, Octave's strncmp function
%% returns 1 if the strings are equal, and 0 otherwise.  This is
%% just the opposite of the corresponding C library function.

%% Author: Bill Lash <lash@tellabs.com>

function status = strncmp(s1, s2, n)

  status = 0;		% Assume strings are different
  if (isstr(s1) & isstr(s2))
    c1 = columns (s1);
    c2 = columns (s2);
    if (n < 1)		
			% Comparing less than 1 character of the string
      status = 0;	% will always say they are not equal (for Matlab
		        % compatibility
    elseif ((n <= c1) & (n <= c2))
       status = strcmp(s1(:,1:n),s2(:,1:n));
    else 
       status = 0;
    end
  end
end

