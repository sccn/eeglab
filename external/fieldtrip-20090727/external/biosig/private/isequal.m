%% Copyright (C) 2005, 2006, 2007 William Poetra Yoga Hadisoeseno
%%
%% This file is part of Octave.
%%
%% Octave is free software; you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3 of the License, or (at
%% your option) any later version.
%%
%% Octave is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with Octave; see the file COPYING.  If not, see
%% <http://www.gnu.org/licenses/>.

%% -*- texinfo -*-
%% @deftypefn {Function File} {} isequal (@var{x1}, @var{x2}, @dots{})
%% Return true if all of @var{x1}, @var{x2}, @dots{} are equal.
%% @seealso{isequalwithequalnans}
%% @end deftypefn

%% modified by Alois Schloegl for the use with Matlab

function retval = isequal (x, varargin)

  if (nargin > 1)
    retval = isequal3 (0, x, varargin{:});
  else
    print('usage:  isequal(x,y)');
  end

end

